import inspect
import tempfile
import unittest
from pathlib import Path
import numpy as np
import pandas as pd
from src import outcomes
from src.config import CANDIDATE_GENES,COMPARISON_GENES
from src.covariates import PROLIFERATION_GENES
from src.statistics import bh,fit_model
from src.null_models import matched_random,permutation,permute_labels,bootstrap,leave_one_out,association

class Stage2Tests(unittest.TestCase):
    def test_bh_known_values_and_missing(self):
        np.testing.assert_allclose(bh([.01,.04,.03,np.nan]),[.03,.04,.04,np.nan],equal_nan=True)
        with self.assertRaises(ValueError):bh([-1])
    def test_exact_patient_matching_no_row_order(self):
        f=pd.DataFrame({'HRD_total':[7.,9.],'_patient':['TCGA-AA-0001','TCGA-AA-0002'],
                        '_sample_key':['TCGA-AA-0001-01','TCGA-AA-0002-01'],'_source_id':['TCGA-AA-0001-01A','TCGA-AA-0002-01A']})
        a=outcomes.match_record(f,'TCGA-AA-0002-01','')[0]
        b=outcomes.match_record(f.iloc[::-1],'TCGA-AA-0002-01','')[0]
        self.assertEqual(a.HRD_total,9.);self.assertEqual(a.HRD_total,b.HRD_total)
        self.assertIsNone(outcomes.match_record(f,'TCGA-AA-0003-01','')[0])
    def test_conflicting_outcome_duplicates_rejected(self):
        f=pd.DataFrame({'HRD_total':[7.,9.],'_patient':['TCGA-AA-0001']*2,'_sample_key':['TCGA-AA-0001-01']*2,
                        '_source_id':['TCGA-AA-0001-01A','TCGA-AA-0001-01B']})
        self.assertEqual(outcomes.match_record(f,'TCGA-AA-0001-01','')[1],'conflicting_genomic_records')
    def test_outcomes_independent_of_predictor(self):
        with tempfile.TemporaryDirectory() as d:
            root=Path(d)/'outer'/'V2';p=root/'data/external/tcga_hrd';p.mkdir(parents=True)
            (p/'TCGA.HRD_withSampleID.txt').write_text('sampleID\tHRD\thrd-loh\tlst1\tai1\nTCGA-AA-0001-01\t30\t10\t11\t9\n')
            c=pd.DataFrame({'patient_id':['TCGA-AA-0001'],'clinical_sample_id':['TCGA-AA-0001-01'],'rna_sample_id':['TCGA-AA-0001-01A'],'candidate_mean_z':[100.]})
            a=outcomes.load_outcomes(root,c)[0];c['candidate_mean_z']=-999.
            pd.testing.assert_frame_equal(a,outcomes.load_outcomes(root,c)[0]);self.assertEqual(a.HRD_total.iloc[0],30)
            self.assertTrue(a.FGA.isna().all())
        self.assertNotIn('expression',inspect.getsource(outcomes.load_outcomes))
    def test_no_candidates_in_proliferation(self):
        self.assertEqual(len(PROLIFERATION_GENES),11);self.assertFalse(set(CANDIDATE_GENES)&set(PROLIFERATION_GENES))
    def test_leave_one_out_exact_eight_genes(self):
        z=pd.DataFrame(np.arange(27).reshape(9,3),index=CANDIDATE_GENES)
        result=leave_one_out(z)
        self.assertEqual(len(result),9)
        for g in CANDIDATE_GENES:np.testing.assert_allclose(result[g],z.drop(index=g).mean(axis=0))
    def test_ols_and_hc3_against_matrix_algebra(self):
        rng=np.random.default_rng(4);x=rng.normal(size=80);y=1+2*x+rng.normal(size=80)*(1+abs(x))
        result,fit,_=fit_model(pd.DataFrame({'candidate_mean_z':x,'HRD_total':y}))
        X=np.column_stack([np.ones(len(x)),x]);inv=np.linalg.inv(X.T@X);b=inv@X.T@y;r=y-X@b;h=np.sum((X@inv)*X,axis=1)
        cov=inv@X.T@np.diag(r*r/(1-h)**2)@X@inv
        self.assertAlmostEqual(result['beta'],b[1]);self.assertAlmostEqual(result['hc3_se'],np.sqrt(cov[1,1]))
    def test_insufficient_n_is_not_negative_result(self):
        result,_,_=fit_model(pd.DataFrame({'candidate_mean_z':[1,2],'HRD_total':[np.nan,np.nan]}))
        self.assertEqual(result['status'],'unavailable');self.assertTrue(np.isnan(result['p']))
    def test_permutation_preserves_labels_and_seed(self):
        y=np.arange(20.);rng=np.random.default_rng(4);permuted=permute_labels(y,rng)
        np.testing.assert_equal(np.sort(permuted),y);self.assertFalse(np.array_equal(permuted,y))
        x=np.sin(y);a,ts=permutation(x,y,n=50,seed=9);b,ts2=permutation(x,y,n=50,seed=9)
        self.assertEqual(a,b);np.testing.assert_equal(ts,ts2)
        self.assertAlmostEqual(a['empirical_p'],(1+sum(abs(ts)>=a['observed_abs_t']))/51)
    def test_bootstrap_deterministic(self):
        x=np.arange(30.);y=2*x+np.sin(x)
        a,dist=bootstrap(x,y,n=50,seed=2);b,dist2=bootstrap(x,y,n=50,seed=2)
        self.assertEqual(a,b);self.assertEqual(dist,dist2);self.assertLess(a['ci_low'],2);self.assertGreater(a['ci_high'],2)
    def test_random_sets_fixed_seed_no_forbidden_genes(self):
        rng=np.random.default_rng(3);names=[f'G{i}' for i in range(250)]+list(CANDIDATE_GENES)+list(COMPARISON_GENES)
        e=pd.DataFrame(rng.uniform(.5,5,(len(names),25)),index=names,columns=[f'p{i}' for i in range(25)])
        z=e.sub(e.mean(axis=1),axis=0).div(e.std(axis=1,ddof=0),axis=0);x=z.loc[list(CANDIDATE_GENES)].mean(axis=0);y=rng.normal(size=25)
        a,s=matched_random(e,z,e.columns.tolist(),y,x,n=20,seed=3);b,t=matched_random(e,z,e.columns.tolist(),y,x,n=20,seed=3)
        pd.testing.assert_frame_equal(a,b);self.assertEqual(s,t)
        for genes in a.genes:
            selected=genes.split(';');self.assertEqual(len(set(selected)),9);self.assertFalse(set(selected)&(set(CANDIDATE_GENES)|set(COMPARISON_GENES)))
    def test_association_matches_primary_beta(self):
        x=np.arange(15.);y=3*x+np.cos(x);b,t=association(x,y)
        self.assertAlmostEqual(b,np.polyfit(x,y,1)[0]);self.assertTrue(np.isfinite(t))

if __name__=='__main__':unittest.main()
