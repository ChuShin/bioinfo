import unittest
from filter_variants_tassel import FilterVariantTassel


class MyTestCase(unittest.TestCase):

    def test_get_observed_type_biallelic(self):
        filter_vars = FilterVariantTassel("/Users/kkoh/GIFS/projects/bjun_varcal/test", 3, 23, 0.3)
        num_nulls, num_refs, num_alts, num_typed = filter_vars.get_observed_types('A', 'T',['A/A','T/T','N/N'])
        self.assertEqual(num_nulls, 1)
        self.assertEqual(num_refs, 1)
        self.assertEqual(num_alts, 1)
        self.assertEqual(num_typed, 2)

    def test_get_observed_type_triallelic(self):
        filter_vars = FilterVariantTassel("/Users/kkoh/GIFS/projects/bjun_varcal/test", 3, 23, 0.3)
        num_nulls, num_refs, num_alts, num_typed = filter_vars.get_observed_types('A', 'T,C',['A/A', 'T/T', 'C/C', 'T/C', 'T/A',
                                                                                              'C/A', 'N/N'])
        self.assertEqual(num_nulls, 1)
        self.assertEqual(num_refs, 1)
        self.assertEqual(num_alts, 2)
        self.assertEqual(num_typed, 6)


    def test_get_observed_type_triallelic_indel(self):
        filter_vars = FilterVariantTassel("/Users/kkoh/GIFS/projects/bjun_varcal/test", 3, 23, 0.3)
        num_nulls, num_refs, num_alts, num_typed = filter_vars.get_observed_types('A', 'TT,CC',
                                                                                  ['A/TT', 'A/CC', 'TT/CC', 'CC/CC', 'TT/TT',
                                                                                   'A/A', 'N/N'])
        self.assertEqual(num_nulls, 1)
        self.assertEqual(num_refs, 1)
        self.assertEqual(num_alts, 2)
        self.assertEqual(num_typed, 6)

    def test_get_mutation_type(self):
        filter_vars = FilterVariantTassel("/Users/kkoh/GIFS/projects/bjun_varcal/test", 3, 23, 0.3)
        self.assertEquals(filter_vars.get_mutation_type('A', 'T'), 'simple_snp')
        self.assertEquals(filter_vars.get_mutation_type('A', 'T,G'), 'simple_snp')
        self.assertEquals(filter_vars.get_mutation_type('A', 'T,GG'), 'unknown')
        self.assertEquals(filter_vars.get_mutation_type('A', 'TGG'), 'insertion')
        self.assertEquals(filter_vars.get_mutation_type('ATT', 'TT'), 'deletion')
        self.assertEquals(filter_vars.get_mutation_type('AT', 'TGG'), 'insertion')


    def test_get_consensus_call(self):
        filter_vars = FilterVariantTassel("/Users/kkoh/GIFS/projects/bjun_varcal/test", 3, 23, 0.3)
        self.assertEqual(filter_vars.get_consensus_call(['N/N:0,0', 'N/N:0,0', 'N/N:0,0']), "N/N")
        self.assertEqual(filter_vars.get_consensus_call(['N/N:0,0', 'N/N:0,0', 'a/a:2,0']), "N/N")
        self.assertEqual(filter_vars.get_consensus_call(['N/N:0,0', 't/t:0,2', 'a/a:2,0']), "N/N")
        self.assertEqual(filter_vars.get_consensus_call(['G/G:0,0,5', 't/t:0,2,0', 'a/a:2,0,0']), "N/N")
        self.assertEqual(filter_vars.get_consensus_call(['A/A:5,0', 'C/C:0,10', 'c/c:0,1']), 'C/C')
        self.assertEqual(filter_vars.get_consensus_call(['c/c:0,2', 'C/C:0,10', 'c/c:0,1']), 'C/C')
        #triallelic
        self.assertEqual(filter_vars.get_consensus_call(['cc/cc:0,2,0', 'CC/CC:0,10,0', 't/t:0,0,1']), 'CC/CC')
        self.assertEqual(filter_vars.get_consensus_call(['tt/cc:0,2,2', 'CC/CC:0,0,10', 'tt/cc:0,1,1']), 'TT/CC')



if __name__ == '__main__':
    unittest.main()
