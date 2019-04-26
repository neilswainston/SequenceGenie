'''
sbc-ngs (c) University of Manchester 2019

All rights reserved.

@author: neilswainston
'''
# pylint: disable=no-name-in-module
import os
import subprocess

from Bio import Seq, SeqIO, SeqRecord
from synbiochem.utils import io_utils


def index(filename):
    '''Index file.'''
    subprocess.call(['bwa', 'index', filename])


def mem(templ_filename, reads_filename, out_filename=None,
        readtype='ont2d', gap_open=6):
    '''Runs BWA MEM.'''
    out_file = io_utils.get_filename(out_filename)

    with open(out_file, 'w') as out:
        subprocess.call(['bwa', 'mem',
                         '-x', readtype,
                         '-O', str(gap_open),
                         templ_filename, reads_filename],
                        stdout=out)

    return out_file


def reject_indels(sam_filename_in, templ_filename, sam_filename_out):
    '''Rejects indels.'''
    sam_file = Samfile(sam_filename_in, 'r')
    out_file = Samfile(sam_filename_out, 'wh',
                       template=sam_file,
                       header=sam_file.header)
    templ_seq = get_seq(templ_filename)

    all_reads = 0
    passed_reads = 0

    for read in sam_file:
        all_reads += 1

        if read.cigarstring and str(len(templ_seq)) + 'M' in read.cigarstring:
            out_file.write(read)
            passed_reads += 1

    print('%s: %i/%i passed reject_indels filter' % (sam_filename_in,
                                                     passed_reads,
                                                     all_reads))

    out_file.close()


def replace_indels(sam_filename_in, templ_filename, sam_filename_out):
    '''Replace indels, replacing them with wildtype.'''
    sam_filename_out = io_utils.get_filename(sam_filename_out)
    templ_seq = get_seq(templ_filename)
    records = []

    all_reads = 0

    for read in Samfile(sam_filename_in, 'r'):
        # Perform mapping of nucl indices to remove spurious indels:
        all_reads += 1

        seq = ''.join([read.seq[pair[0]]
                       if pair[0]
                       else templ_seq[pair[1]]
                       for pair in read.aligned_pairs
                       if pair[1] is not None])

        if seq:
            records.append(SeqRecord.SeqRecord(Seq.Seq(seq), read.qname,
                                               '', ''))

    reads_filename = io_utils.get_filename(None)

    with open(reads_filename, 'w') as fle:
        SeqIO.write(records, fle, 'fasta')

    mem(templ_filename, reads_filename,
        out_filename=sam_filename_out,
        gap_open=12)

    print('%s: %i/%i passed replace_indels filter' % (sam_filename_in,
                                                      len(records),
                                                      all_reads))

    return sam_filename_out


def get_dir(parent_dir, barcodes, ice_id=None):
    '''Get directory from barcodes.'''
    dir_name = os.path.join(parent_dir, '_'.join(barcodes))

    if ice_id:
        dir_name = os.path.join(dir_name, ice_id)

    if not os.path.exists(dir_name):
        os.makedirs(dir_name)

    return dir_name


def get_seq(filename):
    '''Get sequence from Fasta file.'''
    for record in SeqIO.parse(filename, 'fasta'):
        return record.seq

    return None
