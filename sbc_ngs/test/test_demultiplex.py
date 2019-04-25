'''
SequenceGenie (c) University of Manchester 2018

All rights reserved.

@author: neilswainston
'''
# pylint: disable=wrong-import-order
import os
from shutil import copyfile
import tempfile
import unittest
from Bio import SeqIO
import pandas as pd
from sbc_ngs import demultiplex


class Test(unittest.TestCase):
    '''Class to test utils module.'''

    def test_demuliplex_simple_multi(self):
        '''Test demuliplex method.'''
        self.__test_demuliplex_simple(2)

    def test_demuliplex_complex_multi(self):
        '''Test demuliplex method.'''
        self.__test_demuliplex_complex(2)

    def __test_demuliplex_simple(self, num_threads):
        '''Test demuliplex method.'''
        barcodes = [('AAAAAAGGGGGG', 'AAAAAAGGGGGG')]

        test_dir = os.path.dirname(os.path.realpath(__file__))
        in_dir = tempfile.mkdtemp()

        copyfile(os.path.join(test_dir, 'simple_seqs.fasta'),
                 os.path.join(in_dir, 'simple_seqs.fasta'))

        barcode_seqs = demultiplex.demultiplex(barcodes,
                                               in_dir,
                                               min_length=0,
                                               max_read_files=1,
                                               out_dir=tempfile.mkdtemp(),
                                               tolerance=1,
                                               num_threads=num_threads)

        fasta_file = barcode_seqs[barcodes[0]]
        records = list(SeqIO.parse(fasta_file, 'fasta'))

        self.assertEqual(len(records), 2)

    def __test_demuliplex_complex(self, num_threads):
        '''Test demuliplex method.'''
        directory = os.path.dirname(os.path.realpath(__file__))

        barcodes_df = \
            pd.read_csv(os.path.join(directory, 'barcodes.csv'))

        barcodes = [tuple(pair)
                    for pair in
                    barcodes_df[['forward', 'reverse']].values.tolist()]

        test_dir = os.path.dirname(os.path.realpath(__file__))
        in_dir = tempfile.mkdtemp()

        copyfile(os.path.join(test_dir, 'reads.fasta'),
                 os.path.join(in_dir, 'reads.fasta'))

        barcode_seqs = demultiplex.demultiplex(barcodes,
                                               in_dir,
                                               min_length=0,
                                               max_read_files=1,
                                               out_dir=tempfile.mkdtemp(),
                                               tolerance=8,
                                               num_threads=num_threads)

        bc_pair = ('GAGTCTTGTGTCCCAGTTACCAGG', 'CGGGCCCTTCATCTCTCAGCCGAT')
        fasta_file = barcode_seqs[bc_pair]
        records = list(SeqIO.parse(fasta_file, 'fasta'))

        self.assertEqual(len(records), 261)


if __name__ == "__main__":
    unittest.main()
