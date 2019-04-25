'''
sbc-ngs (c) University of Manchester 2019

All rights reserved.

@author: neilswainston
'''
# pylint: disable=invalid-name
# pylint: disable=too-many-arguments
# pylint: disable=wrong-import-order
from collections import OrderedDict
import os.path
from threading import Thread

from Bio import Seq, SeqIO

import multiprocessing as mp
import pandas as pd
from sbc_ngs import reads


class ReadThread(Thread):
    '''Thread-safe class to write demultiplexed reads to Fasta.'''

    def __init__(self, queue, parent_dir):
        self.__queue = queue
        self.__parent_dir = parent_dir
        self.__files = {}
        Thread.__init__(self)

    def run(self):
        '''Run.'''
        while True:
            task = self.__queue.get()

            if task is None:
                break

            self.__write(task)
            self.__queue.task_done()

        for fle in self.__files.values():
            fle.close()

    def get_filenames(self):
        '''Get filenames.'''
        return {barcodes: fle.name
                for barcodes, fle in self.__files.items()}

    def close(self):
        '''Close.'''
        self.__queue.put(None)

    def __write(self, task):
        barcodes = task[0]

        if barcodes not in self.__files:
            dir_name = os.path.join(self.__parent_dir, '_'.join(barcodes))
            filename = os.path.join(dir_name, 'reads.fasta')

            if not os.path.exists(dir_name):
                os.makedirs(dir_name)

            self.__files[barcodes] = open(filename, 'w')

        SeqIO.write(task[1], self.__files[barcodes], 'fasta')


def get_barcodes(filename):
    '''Get barcodes.'''
    barcodes_df = pd.read_csv(os.path.join(filename))
    barcodes_df.fillna('', inplace=True)

    barcodes = \
        [tuple(pair) for pair in barcodes_df[['forward', 'reverse']].values]

    barcodes_dfs = []

    for barcode_type in ['all', 'forward', 'reverse']:
        barcode_type_df = barcodes_df.copy()
        barcode_type_df['barcode_type'] = barcode_type
        barcodes_dfs.append(barcode_type_df)

    return barcodes, pd.concat(barcodes_dfs)


def demultiplex(barcodes, in_dir, min_length, max_read_files, out_dir,
                tolerance, num_threads, search_len=48):
    '''Bin sequences according to barcodes.'''
    max_barcode_len = max([len(barcode)
                           for pair in barcodes
                           for barcode in pair])

    pool = mp.Pool(processes=num_threads)
    write_queue = mp.Manager().Queue()
    read_thread = ReadThread(write_queue, out_dir)
    read_thread.start()

    filenames = reads.get_filenames(in_dir, max_read_files)

    results = [pool.apply_async(_bin_seqs, args=(fle,
                                                 min_length,
                                                 max_barcode_len,
                                                 search_len,
                                                 _format_barcodes(barcodes),
                                                 tolerance,
                                                 idx,
                                                 len(filenames),
                                                 write_queue))
               for idx, fle in enumerate(filenames)]

    for res in results:
        res.get()

    write_queue.join()
    read_thread.close()
    return read_thread.get_filenames()


def _format_barcodes(barcodes):
    '''Format barcodes to reduce number of get_rev_complement calls.'''
    for_brcds = OrderedDict()
    rev_brcds = OrderedDict()

    for pair in barcodes:
        for_brcds[pair] = \
            [list(pair[0]), list(Seq.Seq(pair[1]).reverse_complement())]
        rev_brcds[pair] = \
            [list(pair[1]), list(Seq.Seq(pair[0]).reverse_complement())]

    return for_brcds, rev_brcds


def _bin_seqs(reads_filename, min_length, max_barcode_len, search_len,
              barcodes, tolerance, idx, num_read_files, write_queue):
    '''Bin a batch of sequences.'''
    barcode_seqs = 0

    seqs = reads.get_reads(reads_filename, min_length)

    for seq in seqs:
        if seq:
            for pairs in barcodes:
                if _check_seq(seq, max_barcode_len, search_len, pairs,
                              tolerance, write_queue):
                    barcode_seqs += 1
                    break

    _report_barcodes(idx, num_read_files, len(seqs), barcode_seqs)


def _check_seq(seq, max_barcode_len, search_len, pairs,
               tolerance, write_queue, indiv_strand=False):
    '''Check sequence against barcode sequences.'''
    seq_len = min(max_barcode_len + search_len, len(seq))
    seq_start = list(seq.seq[:seq_len])
    seq_end = list(seq.seq[-(seq_len):])

    # Check all barcodes:
    for orig, bc_pair in pairs.items():
        selected_barcodes = \
            _check_pair(orig, bc_pair, [seq_start,
                                        seq_end], seq_len, tolerance)

        if selected_barcodes[0] and selected_barcodes[1]:
            write_queue.put([tuple(selected_barcodes + ['all']), seq])

            if indiv_strand:
                if orig[0] == ''.join(bc_pair[0]):
                    write_queue.put([tuple(selected_barcodes + ['forward']),
                                     seq])
                else:
                    write_queue.put([tuple(selected_barcodes + ['reverse']),
                                     seq])

            return True

    return False


def _check_pair(orig, pair, seqs, seq_len, tolerance):
    '''Check similarity scores.'''
    selected_barcodes = [None, None]

    for idx in range(2):
        resp = _check_barcode(orig[idx], pair[idx], seqs[idx], seq_len,
                              tolerance)

        if resp:
            selected_barcodes[idx] = resp

    return selected_barcodes


def _check_barcode(orig, barcode, seq, seq_len, tolerance):
    '''Check barcode.'''
    bc_len = len(barcode)

    for substr in [seq[i:i + bc_len] for i in range(seq_len - bc_len + 1)]:
        diff = 0

        for idx, s in enumerate(substr):
            if s != barcode[idx]:
                diff += 1

                if diff > tolerance:
                    break

        if diff <= tolerance:
            return orig

    return None


def _report_barcodes(idx, num_read_files, num_seqs, barcode_seqs):
    '''Report barcodes.'''
    if barcode_seqs:
        print('Files: %d/%d\tMatched: %d/%d' % ((idx + 1),
                                                num_read_files,
                                                barcode_seqs,
                                                num_seqs))
    else:
        print('Files: %d/%d' % ((idx + 1), num_read_files))
