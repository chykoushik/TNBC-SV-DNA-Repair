import inspect
import unittest
from pathlib import Path
import numpy as np
import pandas as pd
from src.config import CANDIDATE_GENES
from src.stage4_scanb import biological_key, map_samples, require_unambiguous, construct_predictors, sample_null_sets
from src.null_models import permutation, bootstrap


class Stage4Tests(unittest.TestCase):
    def record(self,title='F2632',m='m2'):
        return dict(sample='GSM'+title,title=title,characteristics=[
            'scan-b external id: Q009039.C009084.S000003.l.r.'+m+'.c.lib.g.k.a.t',
            'age at diagnosis: 50'])

    def hrd(self):
        return pd.DataFrame({'PD_ID':['PD1'],'RBA':['S000003.l.r.m2.c.lib.g.k2.a.t']})

    def test_mapping_k_normalization_and_determinism(self):
        records=[self.record(),self.record('F2632repl')]
        a=map_samples(self.hrd(),records,['F2632','F2632repl'])
        b=map_samples(self.hrd(),records[::-1],['F2632repl','F2632'])
        pd.testing.assert_frame_equal(a,b)
        self.assertEqual(a.loc[a.mapping_status=='mapped','F_expression_column'].tolist(),['F2632'])
        self.assertEqual(a.mapping_status.tolist().count('excluded_replicate'),1)
        self.assertEqual(biological_key(self.hrd().RBA.iloc[0]),('S000003','m2'))

    def test_m_variants_not_rescued(self):
        a=map_samples(self.hrd(),[self.record(m='m')],['F2632'])
        self.assertEqual(a.mapping_status.tolist(),['unmatched'])

    def test_ambiguity_rejected(self):
        a=map_samples(self.hrd(),[self.record(),self.record('F99')],['F2632','F99'])
        with self.assertRaisesRegex(ValueError,'Ambiguous'):require_unambiguous(a)
        self.assertFalse((a.mapping_status=='mapped').any())

    def test_duplicate_patient_rejected(self):
        with self.assertRaises(ValueError):map_samples(pd.concat([self.hrd()]*2),[self.record()],['F2632'])

    def expression(self):
        return pd.DataFrame(np.random.default_rng(7).normal(size=(9,20)),index=CANDIDATE_GENES,columns=['F'+str(i) for i in range(20)])

    def test_frozen_genes_and_standardization(self):
        self.assertEqual(tuple(CANDIDATE_GENES),('HORMAD1','HORMAD2','STAG3','REC8','SMC1B','SYCP2','SYCP3','MSH4','MSH5'))
        p,z,_=construct_predictors(self.expression())
        np.testing.assert_allclose(z.mean(axis=1),0,atol=1e-14)
        np.testing.assert_allclose(z.std(axis=1,ddof=0),1)
        np.testing.assert_allclose(p.candidate_mean_z,z.mean(axis=0))
        with self.assertRaises(ValueError):construct_predictors(self.expression(),genes=CANDIDATE_GENES[:-1])
        p,_,_=construct_predictors(self.expression().drop('REC8'))
        self.assertTrue(p.candidate_mean_z.isna().all())

    def test_outcomes_and_receptors_cannot_define_predictors_or_selection(self):
        e=self.expression();p,_,_=construct_predictors(e)
        for g in ['HRDetect_probability','scarHRD_score','ESR1','PGR','ERBB2']:e.loc[g]=np.arange(20)*100
        q,_,_=construct_predictors(e)
        pd.testing.assert_frame_equal(p,q)
        self.assertEqual(list(inspect.signature(construct_predictors).parameters),['expression','genes'])
        selection=inspect.getsource(map_samples)
        for gene in ['ESR1','PGR','ERBB2']:self.assertNotIn(gene,selection)

    def test_retired_score_absent(self):
        for name in ['stage4.py','stage4_scanb.py']:
            self.assertNotIn('gips',(Path(__file__).parents[1]/'src'/name).read_text().lower())

    def test_fixed_seed_controls(self):
        neighbors={g:['R'+str(i) for i in range(30)] for g in CANDIDATE_GENES}
        a=sample_null_sets(neighbors,20,71)
        self.assertEqual(a,sample_null_sets(neighbors,20,71))
        self.assertTrue(all(len(set(s))==9 for s in a))
        rng=np.random.default_rng(1);x=rng.normal(size=50);y=.4*x+rng.normal(size=50)
        for method in [permutation,bootstrap]:
            r,a=method(x,y,n=50,seed=21);s,b=method(x,y,n=50,seed=21)
            self.assertEqual(r,s);np.testing.assert_array_equal(a,b)


if __name__=='__main__':unittest.main()
