'''
sbc-ngs (c) University of Manchester 2019

All rights reserved.

@author: neilswainston
'''
# pylint: disable=no-name-in-module
import os
import subprocess

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


def get_dir(parent_dir, barcodes, ice_id=None):
    '''Get directory from barcodes.'''
    dir_name = os.path.join(parent_dir, '_'.join(barcodes))

    if ice_id:
        dir_name = os.path.join(dir_name, ice_id)

    if not os.path.exists(dir_name):
        os.makedirs(dir_name)

    return dir_name
