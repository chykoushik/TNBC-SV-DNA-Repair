import copy
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch
import numpy as np
from src.star import parse_star, resolve_metadata, select_replicate
from src import stage1
from src.signatures import fit_standardization, score_predictors
from src.config import CANDIDATE_GENES

HEADER='gene_id\tgene_name\tgene_type\tunstranded\tstranded_first\tstranded_second\ttpm_unstranded\tfpkm_unstranded\tfpkm_uq_unstranded\n'

class StarTests(unittest.TestCase):
    def parse(self,body):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/'star.tsv';p.write_text('# gene-model: GENCODE v36\n'+HEADER+body)
            return parse_star(p)
    def test_augmented_tpm_not_stranded_counts(self):
        genes,tpm,counts,counters,comments=self.parse('N_unmapped\t\t\t10\t10\t10\t\t\t\nENSG000001.2\tHORMAD1\tprotein_coding\t100\t90\t10\t3\t99\t999\n')
        self.assertEqual(genes,[('ENSG000001.2','HORMAD1','protein_coding')])
        self.assertEqual(tpm.tolist(),[3.]);self.assertEqual(counts.tolist(),[100.])
        self.assertEqual(counters,{'N_unmapped':10})
        self.assertEqual(np.log2(tpm+1).tolist(),[2.])
    def test_negative_or_nan_tpm_rejected(self):
        for v in ['-1','nan','inf']:
            with self.assertRaises(ValueError):self.parse(f'ENSG000001.2\tHORMAD1\tprotein_coding\t1\t1\t1\t{v}\t1\t1\n')
    def test_duplicate_ensembl_rejected(self):
        row='ENSG000001.2\tHORMAD1\tprotein_coding\t1\t1\t1\t1\t1\t1\n'
        with self.assertRaises(ValueError):self.parse(row+row)
    def test_missing_field_rejected(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/'bad.tsv';p.write_text('gene_id\tunstranded\nENSG1\t1\n')
            with self.assertRaises(ValueError):parse_star(p)
    def fixture(self):
        c={'patient_id':'TCGA-AA-0001','sample_id':'TCGA-AA-0001-01'}
        m={'file_id':'f','md5sum':'m','cases':[{'submitter_id':c['patient_id'],'samples':[{'submitter_id':'TCGA-AA-0001-01A','sample_type':'Primary Tumor'}]}]}
        return c,m
    def test_match_restricted_to_clinical_sample(self):
        c,m=self.fixture();self.assertEqual(resolve_metadata(m,c),('TCGA-AA-0001-01A',''))
        c['patient_id']='TCGA-AA-0002'
        with self.assertRaises(ValueError):resolve_metadata(m,c)
    def test_nonprimary_and_wrong_vial_rejected(self):
        c,m=self.fixture();c['sample_id']='TCGA-AA-0001-01B'
        with self.assertRaises(ValueError):resolve_metadata(m,c)
        c,m=self.fixture();m['cases'][0]['samples'][0]['sample_type']='Solid Tissue Normal'
        with self.assertRaises(ValueError):resolve_metadata(m,c)
    def test_multiple_sample_metadata_rejected(self):
        c,m=self.fixture();m['cases'][0]['samples']*=2
        with self.assertRaises(ValueError):resolve_metadata(m,c)
    def test_full_aliquot_preserved(self):
        c,m=self.fixture();d={'data':copy.deepcopy(m)};a='TCGA-AA-0001-01A-01R-A000-01'
        d['data']['cases'][0]['samples'][0]['portions']=[{'analytes':[{'aliquots':[{'submitter_id':a}]}]}]
        self.assertEqual(resolve_metadata(m,c,d)[1],a)
    def test_duplicate_choice_is_order_independent(self):
        rows=[dict(file_id='b',assigned_reads=100,gdc_sample_id='s',aliquot_id='a'),dict(file_id='a',assigned_reads=100,gdc_sample_id='s',aliquot_id='a')]
        self.assertEqual(select_replicate(rows)['file_id'],'a')
        self.assertEqual(select_replicate(rows),select_replicate(rows[::-1]))
        rows[0]['assigned_reads']=101;self.assertEqual(select_replicate(rows)['file_id'],'b')
    def test_only_discovery_reference_used(self):
        x={g:{'tnbc1':0.,'tnbc2':2.,'non_tnbc':10000.} for g in CANDIDATE_GENES}
        p=fit_standardization(x,['tnbc1','tnbc2'])
        self.assertEqual(score_predictors(x,['tnbc2'],p)[0]['candidate_mean_z'],1.)
        self.assertEqual(p['HORMAD1']['n'],2)
    def test_continuation_does_not_rebuild_clinical(self):
        with patch('src.star.run',return_value='continued') as run,patch.object(stage1,'build_initial_cohort') as build:
            self.assertEqual(stage1.run(),'continued');run.assert_called_once();build.assert_not_called()

if __name__=='__main__':unittest.main()
