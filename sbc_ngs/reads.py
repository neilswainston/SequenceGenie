'''
sbc-ngs (c) University of Manchester 2019

All rights reserved.

@author: neilswainston
'''
import os
from os.path import splitext

from Bio import SeqIO


def get_filenames(reads_filename, max_read_files=1e16):
    '''Get reads filename.'''
    filenames = []

    if os.path.isdir(reads_filename):
        for filename in os.listdir(os.path.abspath(reads_filename)):
            filenames.append(os.path.join(reads_filename, filename))

        return filenames[:max_read_files]

    # else:
    return [reads_filename]


def get_reads(filename, min_length):
    '''Gets reads.'''
    _, ext = splitext(filename)

    try:
        with open(filename, 'rU') as fle:
            print('Reading: %s' % filename)

            all_reads = [record for record in SeqIO.parse(fle, ext[1:])]

            return [record for record in all_reads
                    if len(record.seq) > min_length]
    except (IOError, ValueError) as err:
        print(err)
        return []
