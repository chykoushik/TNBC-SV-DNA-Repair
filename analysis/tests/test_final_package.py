import unittest,sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
from validate_final_package import forbidden_fields,FROZEN
from build_final_package import FIGURES,TABLES,value
from src.config import CANDIDATE_GENES

class FinalPackageTests(unittest.TestCase):
    def test_frozen_signature(self):self.assertEqual(tuple(CANDIDATE_GENES),FROZEN)
    def test_retired_field_recursive(self):
        self.assertTrue(forbidden_fields({'results':[{'old_gips_score':1}]}))
        self.assertFalse(forbidden_fields({'methods':'GIPS was not used'}))
    def test_expected_scope(self):
        self.assertEqual(len(FIGURES),5);self.assertEqual(len(TABLES),9)
    def test_missing_not_zero(self):
        self.assertEqual(value({'a':float('nan'),'b':0},'a','b'),0)

if __name__=='__main__':unittest.main()
