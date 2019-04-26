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
import sys
import uuid

import pysam

import multiprocessing as mp
import pandas as pd
from sbc_ngs import demultiplex, ice, results, utils, vcf_utils


class PathwayAligner():
    '''Class to align NGS data to pathways.'''

    def __init__(self, out_dir, in_dir,
                 ice_url, ice_username, ice_password,
                 for_primer, rev_primer, min_length, max_read_files,
                 dp_filter=0.25):
        # Initialise project directory:
        self.__out_dir = os.path.join(out_dir, str(uuid.uuid4()))
        os.makedirs(self.__out_dir)

        # Get pathway sequences from ICE:
        self.__ice_files, _, _ = \
            ice.get_ice_files(ice_url, ice_username, ice_password,
                              os.path.join(in_dir, 'ice_ids.txt'),
                              for_primer, rev_primer,
                              self.__out_dir)

        self.__in_dir = in_dir
        self.__min_length = min_length
        self.__max_read_files = max_read_files
        self.__dp_filter = dp_filter

        self.__barcodes, self.__barcodes_df = \
            demultiplex.get_barcodes(os.path.join(in_dir, 'barcodes.csv'))

    def score_alignments(self, tolerance, num_threads):
        '''Score alignments.'''
        for templ_filename, _ in self.__ice_files.values():
            utils.index(templ_filename)

        barcode_reads = demultiplex.demultiplex(self.__barcodes,
                                                self.__in_dir,
                                                self.__min_length,
                                                self.__max_read_files,
                                                self.__out_dir,
                                                tolerance=tolerance,
                                                num_threads=num_threads)

        pool = mp.Pool(processes=num_threads)
        write_queue = mp.Manager().Queue()
        results_thread = results.ResultsThread(sorted(self.__ice_files.keys()),
                                               self.__barcodes_df,
                                               write_queue)
        results_thread.start()

        rslts = [pool.apply_async(_score_alignment,
                                  args=(self.__out_dir,
                                        barcodes,
                                        reads_filename,
                                        self.__get_ice_files(barcodes),
                                        self.__dp_filter,
                                        write_queue))
                 for barcodes, reads_filename in barcode_reads.items()]

        for res in rslts:
            res.get()

        # Update summary:
        results_thread.close()
        results_thread.write(self.__out_dir)

    def __get_ice_files(self, barcodes):
        '''Get appropriate ICE files.'''
        try:
            ice_id = self.__barcodes_df.loc[barcodes, 'actual_ice_id']

            if ice_id:
                return {ice_id: self.__ice_files[ice_id]}
        except KeyError:
            print('Unexpected barcodes: ' + str(barcodes))
            return {}

        return self.__ice_files


def _get_barcode_ice(barcode_ice_filename):
    '''Get barcode ice dict.'''
    barcode_ice = pd.read_csv(barcode_ice_filename,
                              dtype={'barcode': str, 'ice_id': str}) \
        if barcode_ice_filename else None

    return barcode_ice.set_index('barcode')['ice_id'].to_dict()


def _score_alignment(dir_name, barcodes, reads_filename, ice_files,
                     dp_filter, write_queue):
    '''Score an alignment.'''
    for ice_id, (templ_filename, _) in ice_files.items():
        _score_barcodes_ice(templ_filename, dir_name, barcodes,
                            ice_id, reads_filename,
                            dp_filter, write_queue)

        print('Scored: %s against %s' % (reads_filename, ice_id))


def _score_barcodes_ice(templ_pcr_filename, dir_name, barcodes,
                        ice_id, reads_filename,
                        dp_filter, write_queue):
    '''Score barcodes ice pair.'''
    barcode_dir_name = utils.get_dir(dir_name, barcodes, ice_id)
    sam_filename = os.path.join(barcode_dir_name, 'alignment.sam')
    bam_filename = os.path.join(barcode_dir_name, 'alignment.bam')

    # Align:
    utils.mem(templ_pcr_filename, reads_filename, sam_filename)

    # Convert sam to bam and sort:
    pysam.view(sam_filename, '-o', bam_filename, catch_stdout=False)
    pysam.sort('-o', bam_filename, bam_filename)
    os.remove(sam_filename)

    # Generate and analyse variants file:
    vcf_filename = vcf_utils.get_vcf(bam_filename, templ_pcr_filename,
                                     pcr_offset=0)

    vcf_utils.analyse(vcf_filename, ice_id, barcodes, dp_filter, write_queue)


def main(args):
    '''main method.'''
    try:
        num_threads = int(args[-1])
    except ValueError:
        num_threads = mp.cpu_count()

    print('Running pathway with %d threads' % num_threads)

    aligner = PathwayAligner(*args[:-4],
                             min_length=int(args[-4]),
                             max_read_files=int(args[-3]))

    aligner.score_alignments(int(args[-2]), num_threads)


if __name__ == '__main__':
    main(sys.argv[1:])
