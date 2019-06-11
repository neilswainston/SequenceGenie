'''
sbc-ngs (c) University of Manchester 2019

All rights reserved.

@author: neilswainston
'''
# pylint: disable=no-member
# pylint: disable=too-few-public-methods
# pylint: disable=too-many-arguments
# pylint: disable=too-many-instance-attributes
# pylint: disable=unused-argument
# pylint: disable=wrong-import-order
from __future__ import division

import os
import subprocess
import sys
import uuid

import multiprocessing as mp
import pandas as pd
from sbc_ngs import demultiplex, results, utils, vcf_utils


class PathwayAligner():
    '''Class to align NGS data to pathways.'''

    def __init__(self, out_dir, in_dir, seq_files, min_length, max_read_files):
        # Initialise project directory:
        self.__out_dir = out_dir

        if not os.path.exists(self.__out_dir):
            os.makedirs(self.__out_dir)

        self.__in_dir = in_dir
        self.__seq_files = seq_files
        self.__min_length = min_length
        self.__max_read_files = max_read_files

        self.__barcodes, self.__barcodes_df = \
            demultiplex.get_barcodes(os.path.join(in_dir, 'barcodes.csv'))

        # Backwards compatibility:
        self.__barcodes_df.rename(columns={'actual_ice_id': 'known_seq_id'},
                                  inplace=True)
        self.__barcodes_df['known_seq_id'] = \
            self.__barcodes_df['known_seq_id'].astype(str)

        # Index sequence / template files:
        for templ_filename in self.__seq_files.values():
            subprocess.call(['bwa', 'index', templ_filename])

    def score_alignments(self, tolerance, num_threads):
        '''Score alignments.'''
        num_threads = num_threads if num_threads > 0 else mp.cpu_count()
        print('Running pathway with %d threads' % num_threads)

        barcode_reads = demultiplex.demultiplex(self.__barcodes,
                                                self.__in_dir,
                                                self.__min_length,
                                                self.__max_read_files,
                                                self.__out_dir,
                                                tolerance=tolerance,
                                                num_threads=num_threads)

        write_queue = mp.Manager().Queue()
        results_thread = results.ResultsThread(sorted(self.__seq_files.keys()),
                                               self.__barcodes_df,
                                               write_queue)
        results_thread.start()

        for barcodes, reads_filename in barcode_reads.items():
            _score_alignment(self.__out_dir,
                             barcodes,
                             reads_filename,
                             self.__get_seq_files(barcodes),
                             num_threads,
                             write_queue)

        # Update summary:
        results_thread.close()
        results_thread.write(self.__out_dir)

    def __get_seq_files(self, barcodes):
        '''Get appropriate sequence files.'''
        try:
            seq_id = self.__barcodes_df.loc[barcodes, 'known_seq_id']

            if seq_id:
                return {seq_id: self.__seq_files[seq_id]}
        except KeyError:
            print('Unexpected barcodes: ' + str(barcodes))
            return {}

        return self.__seq_files


def _get_barcode_seq(barcode_seq_filename):
    '''Get barcode seq dict.'''
    barcode_seq = pd.read_csv(barcode_seq_filename,
                              dtype={'barcode': str, 'seq_id': str}) \
        if barcode_seq_filename else None

    return barcode_seq.set_index('barcode')['seq_id'].to_dict()


def _score_alignment(dir_name, barcodes, reads_filename, seq_files,
                     num_threads, write_queue):
    '''Score an alignment.'''
    for seq_id, seq_filename in seq_files.items():
        barcode_dir_name = utils.get_dir(dir_name, barcodes, seq_id)
        bam_filename = os.path.join(barcode_dir_name, '%s.bam' % barcodes[2])
        vcf_filename = bam_filename.replace('.bam', '.vcf')

        prc = subprocess.Popen(('bwa', 'mem',
                                '-x', 'ont2d',
                                '-O', '6',
                                '-t', str(num_threads),
                                seq_filename, reads_filename),
                               stdout=subprocess.PIPE)

        subprocess.check_output(('samtools', 'sort',
                                 '-@%i' % num_threads,
                                 '-o', bam_filename, '-'),
                                stdin=prc.stdout)
        prc.wait()

        # Generate and analyse variants file:
        prc = subprocess.Popen(['samtools',
                                'mpileup',
                                '-uvf',
                                seq_filename,
                                '-t', 'DP',
                                '-o', vcf_filename,
                                bam_filename])

        prc.communicate()

        vcf_utils.analyse(vcf_filename, seq_id, barcodes, write_queue)

        print('Scored: %s against %s' % (reads_filename, seq_id))


def _get_seq_files(filename):
    '''Get seq files.'''
    seq_files = {}

    if os.path.isdir(filename):
        for fle in os.listdir(filename):
            name, ext = os.path.splitext(os.path.basename(fle))

            if ext == '.fasta':
                seq_files[name] = os.path.join(filename, fle)
    else:
        seq_files[os.path.splitext(os.path.basename(filename))[0]] = filename

    return seq_files


def main(args):
    '''main method.'''
    seq_files = {}

    for seq_file in args[6:]:
        seq_files.update(_get_seq_files(seq_file))

    aligner = PathwayAligner(out_dir=os.path.join(args[0], str(uuid.uuid4())),
                             in_dir=args[1],
                             seq_files=seq_files,
                             min_length=int(args[2]),
                             max_read_files=int(args[3]))

    aligner.score_alignments(int(args[4]), num_threads=int(args[5]))


if __name__ == '__main__':
    main(sys.argv[1:])
