import gzip
import tempfile
import unittest
from pathlib import Path
import numpy as np
import pandas as pd
from src.config import CANDIDATE_GENES
from src.stage3_patient import signature,gene_mean,cox,geo_metadata,gse25066
from src.stage3_single_cell import accumulate_triplets,paired_result
from src.stage3 import evidence_row

class Stage3Tests(unittest.TestCase):
    def test_same_nine_genes_equal_weight(self):
        x=pd.DataFrame({f's{i}':np.arange(9)+i for i in range(4)},index=CANDIDATE_GENES)
        score,parameters,missing=signature(x,CANDIDATE_GENES)
        expected=(np.arange(4)-1.5)/np.std(np.arange(4))
        np.testing.assert_allclose(score,expected);self.assertEqual(missing,[]);self.assertEqual(parameters['genes'],list(CANDIDATE_GENES))
    def test_missing_gene_is_not_zero_or_dropped(self):
        x=pd.DataFrame({'a':np.arange(8),'b':np.arange(8)+1},index=CANDIDATE_GENES[:-1])
        score,_,missing=signature(x,CANDIDATE_GENES)
        self.assertTrue(score.empty);self.assertEqual(missing,[CANDIDATE_GENES[-1]])
    def test_duplicate_probes_averaged_not_selected(self):
        frame=pd.DataFrame({'Hugo_Symbol':['STAG3','STAG3','REC8'],'s1':[2.,10.,5.],'s2':[4.,12.,7.]})
        self.assertEqual(gene_mean(frame).loc['STAG3','s1'],6.)
        pd.testing.assert_frame_equal(gene_mean(frame),gene_mean(frame.iloc[::-1]))
    def test_sparse_pseudobulk_matches_dense_sum(self):
        triples=np.array([[1,1,2],[2,1,3],[3,2,7],[1,3,4],[2,3,5]])
        groups=np.array([0,1,0]);lookup=np.array([0,-1,1]);totals=np.zeros(2);counts=np.zeros((2,2))
        accumulate_triplets(triples,groups,lookup,totals,counts)
        np.testing.assert_equal(totals,[14,7]);np.testing.assert_equal(counts,[[6,0],[0,7]])
    def test_sparse_bad_coordinate_rejected(self):
        with self.assertRaises(ValueError):accumulate_triplets(np.array([[1,0,2]]),np.array([0]),np.array([0]),np.zeros(1),np.zeros((1,1)))
    def test_paired_unit_n_and_determinism(self):
        a=paired_result([1,2,3,4,5,6,7,8],'test');b=paired_result([1,2,3,4,5,6,7,8],'test')
        self.assertEqual(a,b);self.assertEqual(a['N'],8);self.assertEqual(a['effect'],4.5);self.assertAlmostEqual(a['p'],.0078125)
    def test_all_zero_paired_difference_is_negative_result(self):
        r=paired_result([0]*8,'test');self.assertEqual(r['p'],1);self.assertEqual(r['effect'],0)
    def test_cox_known_direction(self):
        rng=np.random.default_rng(22);x=rng.normal(size=200);event_time=rng.exponential(size=200)/np.exp(x);censor=np.ones(200)*2
        d=pd.DataFrame({'time':np.minimum(event_time,censor),'event':(event_time<censor).astype(int),'candidate_perSD':x})
        r=cox(d,'time','event');self.assertEqual(r['status'],'estimated');self.assertGreater(r['beta'],0);self.assertGreater(r['hr_ci_low'],1)
    def test_platform_570_never_substituted(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/'geo.gz'
            with gzip.open(p,'wt') as f:f.write('!Sample_geo_accession\t"GSM1"\n!Sample_platform_id\t"GPL570"\n!series_matrix_table_begin\n')
            with self.assertRaises(ValueError):geo_metadata(p)
    def test_clinical_geo_not_expression_rescued(self):
        with tempfile.TemporaryDirectory() as d:
            root=Path(d)/'scaffold'/'V2';(root/'outputs/v2/stage3').mkdir(parents=True)
            source=root.parent.parent/'dataset';source.mkdir()
            fields={'er_status_ihc':['I','N'],'pr_status_ihc':['N','N'],'her2_status':['N','N'],
                    'er_status_ihc_esr1_for indeterminate':['N','N'],'pathologic_response_pcr_rd':['pCR','RD'],
                    'drfs_even_time_years':['1','2'],'drfs_1_event_0_censored':['0','1']}
            with gzip.open(source/'GSE25066_series_matrix.txt.gz','wt') as f:
                f.write('!Sample_geo_accession\t"GSM1"\t"GSM2"\n!Sample_platform_id\t"GPL96"\t"GPL96"\n')
                for k,v in fields.items():f.write('!Sample_characteristics_ch1\t'+ '\t'.join('"'+k+': '+x+'"' for x in v)+'\n')
                f.write('!series_matrix_table_begin\n')
            r,c,s=gse25066(root);self.assertEqual(s['clinical_tnbc_n'],1);self.assertEqual(c.GSM.tolist(),['GSM2']);self.assertEqual(r[0]['status'],'unavailable')
    def test_unavailable_is_not_failed_replication(self):
        r=evidence_row('HORMAD1','DepMap','dependency');self.assertEqual(r['support'],'unavailable');self.assertTrue(np.isnan(r['p']))
    def test_evidence_uses_fdr_and_direction(self):
        r=dict(status='estimated',beta=1,p=.001,fdr=.2,N=10)
        self.assertEqual(evidence_row('HORMAD1','TCGA_HRD','HRD_total',r)['support'],'does_not_support')
        r['fdr']=.01;r['beta']=-1
        self.assertEqual(evidence_row('HORMAD1','TCGA_HRD','HRD_total',r)['support'],'does_not_support')

if __name__=='__main__':unittest.main()
