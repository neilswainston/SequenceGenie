'''
sbc-ngs (c) University of Manchester 2019

All rights reserved.

@author: neilswainston
'''
from os.path import splitext

from Bio import SeqIO


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
