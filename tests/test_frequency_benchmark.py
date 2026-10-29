"""Small independent checks of weighting, leakage and assessment semantics."""
import sys
import unittest
from pathlib import Path
import numpy as np
import scipy.sparse as sp
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from frequency_benchmark import (curated_counts, excluded_reference, frequency_value,
                                 normalize, rank_interval, resolved_weight)
from rank_engine import Engine


class FrequencyChecks(unittest.TestCase):
    def test_pooled_reference_excluded_if_any_component_is_held_out(self):
        self.assertTrue(excluded_reference('PMID:11;PMID:22',{'22'},{},set()))
        self.assertTrue(excluded_reference('PMID:11',set(),{'11':'author'},{'author'}))
        self.assertFalse(excluded_reference('PMID:11',{'22'},{},set()))

    def test_missing_and_conflicting_values_are_not_zero(self):
        self.assertEqual(resolved_weight({''}),1.)
        self.assertEqual(resolved_weight({'0/4'}),.1)
        self.assertEqual(resolved_weight({'1/4','3/4'}),1.)
        self.assertEqual(resolved_weight({'1/4',''}),.3)
        self.assertEqual(frequency_value('HP:0040282','low'),.3)
        self.assertEqual(frequency_value('HP:0040282','high'),.79)

    def test_family_weights_preserve_unknown(self):
        annotations=[dict(cohort_id='c',disease_id='d')]
        meta=[dict(cohort_id='c',patient_label=str(i),family_id=f) for i,f in [(1,'a'),(2,'a'),(3,'b')]]
        assessments=[dict(cohort_id='c',patient_label=str(i),hpo_id='h',pmid='p',status=s)
                     for i,s in [(1,'present'),(2,'absent'),(3,'unknown')]]
        self.assertEqual(curated_counts(annotations,assessments,meta),{('d','h'):[1.,2.]})
        self.assertEqual(curated_counts(annotations,assessments,meta,family=True),{('d','h'):[.5,1.]})

    def test_weighted_bma_against_hand_calculation(self):
        # Tree: root -> a -> b, root -> c. Query a; profile [b,c].
        # Resnik matches [2,0]; forward=2; reverse=(3*2+1*0)/4=1.5.
        e=Engine.__new__(Engine)
        e.names=['root','a','b','c'];e.tid={t:i for i,t in enumerate(e.names)}
        e.anc_sets={'root':{'root'},'a':{'root','a'},'b':{'root','a','b'},'c':{'root','c'}}
        e.A=sp.csr_matrix([[1,0,0,0],[1,1,0,0],[1,1,1,0],[1,0,0,1]],dtype=float)
        e.ic=np.array([0.,2.,3.,4.]);e.genes=['d'];e.gid={'d':0}
        e.build_gene_matrices({'d':['b','c']})
        self.assertAlmostEqual(e.bma_all_genes(['a'])[0],1.5)
        uniform=normalize(sp.csr_matrix([[0,0,1,1]],dtype=float))
        np.testing.assert_allclose(e.bma_all_genes(['a'],GTnorm=uniform),e.bma_all_genes(['a']))
        weighted=normalize(sp.csr_matrix([[0,0,3,1]],dtype=float))
        self.assertAlmostEqual(e.bma_all_genes(['a'],GTnorm=weighted)[0],1.75)

    def test_ties_use_midrank(self):
        self.assertEqual(rank_interval(np.array([5.,3.,3.,1.]),1),(2,3,2.5))


if __name__=='__main__':unittest.main()
