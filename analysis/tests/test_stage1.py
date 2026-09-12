import ast
import copy
import inspect
from pathlib import Path
import statistics
import tempfile
import unittest
from src import cohort, config, expression, signatures
from src.utils import tcga_ids, open_text, output_path

def clinical(sample='TCGA-AA-0001-01', er='Negative', pr='Negative', her2='Negative', ish=''):
    return dict(sampleID=sample, sample_type='Primary Tumor', sample_type_id='01',
                **{config.ER_COLUMN:er, config.PR_COLUMN:pr, config.HER2_COLUMN:her2, config.HER2_ISH_COLUMN:ish})

class Stage1Tests(unittest.TestCase):
    def test_expression_cannot_select_tnbc(self):
        rows=[clinical(),clinical('TCGA-AA-0002-01',er='Positive')]
        expected=cohort.build_cohort(rows)
        for r in rows:r.update(ESR1=-999999,PGR=999999,ERBB2=0,ESR1_median=0,PAM50='Basal')
        self.assertEqual(expected,cohort.build_cohort(rows))
        for r in rows:r.update(ESR1=999999,PGR=-999999,ERBB2=999999)
        self.assertEqual(expected,cohort.build_cohort(rows))
        self.assertNotIn('expression', inspect.signature(cohort.build_cohort).parameters)
    def test_unknown_equivocal_are_not_negative(self):
        for value in ['', 'Unknown', 'Equivocal', 'Indeterminate', 'NA']:
            for key in ['er','pr','her2']:
                self.assertEqual(cohort.build_cohort([clinical(**{key:value})])[0],[])
    def test_ish_conflict_and_no_automatic_rescue(self):
        self.assertFalse(cohort.build_cohort([clinical(ish='Positive')])[0])
        self.assertFalse(cohort.build_cohort([clinical(her2='Equivocal',ish='Negative')])[0])
        self.assertTrue(cohort.build_cohort([clinical(ish='')])[0])
    def test_primary_only(self):
        normal=clinical('TCGA-AA-0001-11');normal.update(sample_type='Solid Tissue Normal',sample_type_id='11')
        self.assertFalse(cohort.build_cohort([normal])[0])
        self.assertFalse(cohort.build_cohort([clinical('TCGA-AA-0001-06')])[0])
        unknown=clinical('TCGA-AA-0001');unknown.update(sample_type='',sample_type_id='')
        self.assertFalse(cohort.build_cohort([unknown])[0])
    def test_deterministic_ids_and_duplicates(self):
        rows=[clinical('TCGA-AA-0001-01B'),clinical('TCGA-AA-0001-01A'),clinical('TCGA-AA-0002-01')]
        self.assertEqual(cohort.build_cohort(rows),cohort.build_cohort(list(reversed(rows))))
        selected,_=cohort.build_cohort(rows)
        self.assertEqual(selected[0]['sample_id'],'TCGA-AA-0001-01A')
        self.assertEqual(tcga_ids('tcga-aa-0001-01a-01r-0000-01')['patient_id'],'TCGA-AA-0001')
    def test_discordant_patient_excluded(self):
        self.assertFalse(cohort.build_cohort([clinical(),clinical('TCGA-AA-0001-01B',er='Positive')])[0])
    def test_matching_unique_and_ambiguous(self):
        c,_=cohort.build_cohort([clinical()])
        a=['TCGA-AA-0001-01A-01R-0000-01']
        self.assertEqual(expression.match_samples(c,a),{'TCGA-AA-0001-01':a[0]})
        a+=['TCGA-AA-0001-01B-01R-0000-01']
        self.assertEqual(expression.match_samples(c,a),{})
        self.assertEqual(expression.match_samples(c,list(reversed(a))),{})
        self.assertEqual(expression.match_samples(c,['TCGA-AA-0001-11']),{})
    def test_predictor_math_and_frozen_reference(self):
        x={g:{'a':0.,'b':2.,'heldout':3.} for g in config.CANDIDATE_GENES}
        p=signatures.fit_standardization(x,['b','a'])
        scores=signatures.score_predictors(x,['heldout'],p)[0]
        for c in signatures.PREDICTOR_COLUMNS:self.assertEqual(scores[c],2.)
        self.assertEqual(p['HORMAD1']['mean'],1.)
        self.assertEqual(p['HORMAD1']['sd'],1.)
    def test_incomplete_signature_is_not_zero(self):
        x={g:{'a':0.,'b':2.} for g in config.CANDIDATE_GENES if g!='REC8'}
        p=signatures.fit_standardization(x,['a','b'])
        r=signatures.score_predictors(x,['a'],p)[0]
        self.assertIsNone(r['candidate_mean_z']);self.assertIsNone(r['candidate_median_z'])
        self.assertIsNone(r['REC8_z']);self.assertEqual(r['HORMAD1_z'],-1.)
    def test_constant_gene_is_unscorable(self):
        x={g:{'a':1.,'b':1.} for g in config.CANDIDATE_GENES}
        p=signatures.fit_standardization(x,['a','b'])
        self.assertIsNone(signatures.score_predictors(x,['a'],p)[0]['candidate_mean_z'])
    def test_no_retired_framework_in_executable_source(self):
        # Executable token inspection
        for path in (config.PROJECT_ROOT/'src').glob('*.py'):
            tree=ast.parse(path.read_text(encoding='utf-8'))
            for node in ast.walk(tree):
                if isinstance(node,(ast.Name,ast.FunctionDef,ast.ClassDef,ast.Attribute)):
                    identifier=getattr(node,'id',getattr(node,'name',getattr(node,'attr','')))
                    self.assertNotIn('gips',identifier.lower())
        self.assertNotIn('median', inspect.getsource(cohort))
    def test_candidates_are_predictors_only(self):
        self.assertFalse(set(config.CANDIDATE_GENES)&set(config.COMPARISON_GENES))
        for func in [signatures.fit_standardization,signatures.score_predictors]:
            self.assertFalse({'outcome','hrd','cin','sv','y'} & set(inspect.signature(func).parameters))
        x={g:{'a':0.,'b':1.} for g in config.CANDIDATE_GENES}
        r=signatures.score_predictors(x,['a'],signatures.fit_standardization(x,['a','b']))[0]
        self.assertEqual(r['predictor_role'],'expression_predictor')
        self.assertFalse(any(k.lower() in {'hrd','cin','sv','outcome','genomic_instability'} for k in r))
    def test_plain_gz_and_gene_loader(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/'verified_gene.tsv.gz'
            p.write_text('gene\tTCGA-AA-0001-01\nHORMAD1\t3\n')
            ids,values=expression.load_gene_expression(p,['HORMAD1'],provenance='synthetic test',transform='log2p1')
            self.assertEqual(values['HORMAD1'][ids[0]],2.)
            with self.assertRaises(ValueError):expression.load_gene_expression(p,['HORMAD1'],provenance=None,transform='identity')
    def test_exon_loader_rejected(self):
        with self.assertRaises(ValueError):expression.load_gene_expression('exon.tsv',[],provenance='x',transform='identity')
    def test_cannot_write_frozen_source(self):
        with self.assertRaises(ValueError):output_path('../../dataset/forbidden.csv')

if __name__=='__main__':unittest.main()
