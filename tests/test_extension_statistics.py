"""Independent checks of inferential units and source/frequency semantics."""
import sys,unittest
from pathlib import Path
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from source_robustness import cluster_boot
from gene_regional_analysis import permutation,eligible
from hpoa_provenance_audit import frequency
from genetic_context import protein_class
from paired_source_analysis import metrics, paired_intervals
from curate_frequency_tables import binary_counts

class StatisticalChecks(unittest.TestCase):
    def test_binary_table_uses_only_explicit_assessments(self):
        self.assertEqual(binary_counts(['+','-','NA','','?']),(1,2,3))
        self.assertEqual(binary_counts(['-','−']),(0,2,0))
        with self.assertRaises(ValueError):binary_counts(['+','not described'])

    def test_protein_consequence_preserves_unknown_and_syntax_scope(self):
        for value,expected in [('Arg41Trp','amino_acid_substitution'),('Gln44Ter','stop_gained'),
                               ('Glu45AspfsTer8','frameshift'),('Arg12_Gly14delinsSer','protein_delins'),
                               ('Met1? ',None),('Gly12=', 'synonymous'),('Ter177GlnextTer20','protein_extension')]:
            self.assertEqual(protein_class('hgvs.p:NP_002428.1:p.'+value.strip()),expected)
        self.assertIsNone(protein_class('hgvs.c:NM_002437.4:c.121C>T'))

    def test_source_share_excludes_mixed_and_unknown_from_denominator(self):
        m=metrics({'A','B','C','D'},{'A':'laboratory','B':'clinical_exam','C':'mixed'})
        self.assertEqual(m['n_terms'],4)
        self.assertEqual(m['test_share'],.5)
        self.assertEqual(m['n_unclassified'],1)
        self.assertTrue(np.isnan(metrics({'D'}, {})['test_share']))

    def test_paired_intervals_preserve_fixed_difference_and_missingness(self):
        n,m,lo,hi=paired_intervals([2,2,np.nan,2],np.random.default_rng(0),100)
        self.assertEqual((n,m,lo,hi),(3,2,2,2))

    def test_cluster_bootstrap_uses_case_weighted_cluster_totals(self):
        # Two equally sampled publications, with unequal numbers of cases.
        values=np.array([[0.],[0.],[1.]])
        ci,n=cluster_boot(values,['A','A','B'],np.random.default_rng(1),500)
        self.assertEqual(n,2)
        np.testing.assert_equal(ci,[[0.],[1.]])
        ci,n=cluster_boot(values,['A','A','A'],np.random.default_rng(1),100)
        self.assertTrue(np.isnan(ci).all())

    def test_stratification_removes_constructed_source_confounding(self):
        cs=[];Y=[]
        for source,nr,ne,y in [('abstract',8,2,1),('fulltext',2,8,0)]:
            for region,n in [('MENA',nr),('EUR',ne)]:
                for i in range(n):
                    cs.append(dict(region=region,corpus=source,group=f'{source}:{region}:{i}'))
                    Y.append([y])
        Y=np.array(Y,float)
        unadjusted=Y[[c['region']=='MENA' for c in cs]].mean()-Y[[c['region']=='EUR' for c in cs]].mean()
        self.assertAlmostEqual(unadjusted,.6)
        obs,*_=permutation(cs,'MENA',lambda c:(c['corpus'],),Y,np.random.default_rng(3),99)
        self.assertAlmostEqual(obs[0],0)

    def test_care_confounding_is_not_declared_estimable(self):
        cs=[dict(region='MENA',care='SA',pmid=str(i)) for i in range(20)]
        cs += [dict(region='EUR',care='US',pmid=str(i+20)) for i in range(20)]
        kept,_,_,status=eligible(cs,'MENA',lambda c:(c['care'],),10,3)
        self.assertEqual(kept,[])
        self.assertNotEqual(status,'eligible')

    def test_frequency_preserves_zero_and_unknown(self):
        self.assertEqual(frequency('0/3'),('count',0,3))
        self.assertEqual(frequency('0/0')[0],'invalid_count')
        self.assertEqual(frequency('4/3')[0],'invalid_count')
        self.assertEqual(frequency(''),('missing',None,None))
        self.assertEqual(frequency('HP:0040282'),('category',None,None))
        self.assertEqual(frequency('42%'),('percentage',None,None))
        self.assertEqual(frequency('142%')[0],'invalid_percentage')

if __name__=='__main__':unittest.main()
