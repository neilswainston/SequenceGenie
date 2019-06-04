'''
SequenceGenie (c) University of Manchester 2018

All rights reserved.

@author: neilswainston
'''
# pylint: disable=invalid-name
# pylint: disable=protected-access
# pylint: disable=wrong-import-order
import os
import unittest

from sbc_ngs import vcf_utils


class Test(unittest.TestCase):
    '''Class to test utils module.'''

    def test_vcf_to_df(self):
        '''Test vcf_to_df method.'''
        directory = os.path.dirname(os.path.realpath(__file__))
        vcf_filename = os.path.join(directory, 'variants_indel.vcf')
        df, _ = vcf_utils.vcf_to_df(vcf_filename)

        self.assertFalse(df.empty)

    def test_analyse_vcf_indel(self):
        '''Test analyse_vcf method.'''
        directory = os.path.dirname(os.path.realpath(__file__))
        vcf_filename = os.path.join(directory, 'variants_indel.vcf')
        _, _, _, indels, _, _, _, _ = vcf_utils.analyse_vcf(vcf_filename)

        self.assertAlmostEqual(indels[0][1], 0.563636, 3)

    def test_analyse_vcf_mut(self):
        '''Test analyse_vcf method.'''
        directory = os.path.dirname(os.path.realpath(__file__))
        vcf_filename = os.path.join(directory, 'variants.vcf')
        _, mutations, _, _, _, _, _, _ = vcf_utils.analyse_vcf(vcf_filename)

        self.assertFalse(mutations)
        # self.assertAlmostEqual(mutations[0][1], 0.510134, 3)

    def test_analyse_vcf_mut_thresh(self):
        '''Test analyse_vcf method.'''
        directory = os.path.dirname(os.path.realpath(__file__))
        vcf_filename = os.path.join(directory, 'variants.vcf')
        _, mutations, _, _, _, _, _, _ = vcf_utils.analyse_vcf(vcf_filename)

        self.assertFalse(mutations)

    def test_get_probs(self):
        '''Test _get_probs method.'''
        num_forw = 5
        num_rev = 8
        reverse = [0.2, 0.4, 0.3, 0.1]
        forward = [0.1, 0.5, 0.2, 0.2]
        prior_nucl = 'A'
        probs = (vcf_utils._get_probs(
            forward, reverse, num_forw, num_rev, prior_nucl))

        for pair in zip(probs, [0.748583, 0.1555631, 0.0648895, 0.0309640]):
            self.assertAlmostEqual(pair[0], pair[1], 5)


if __name__ == "__main__":
    unittest.main()
