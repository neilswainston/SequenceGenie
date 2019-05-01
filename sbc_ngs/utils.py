'''
sbc-ngs (c) University of Manchester 2019

All rights reserved.

@author: neilswainston
'''
# pylint: disable=no-name-in-module
import os


def get_dir(parent_dir, barcodes, ice_id=None):
    '''Get directory from barcodes.'''
    dir_name = os.path.join(parent_dir, '_'.join(barcodes[:2]))

    if ice_id:
        dir_name = os.path.join(dir_name, ice_id)

    if not os.path.exists(dir_name):
        try:
            os.makedirs(dir_name)
        except FileExistsError:
            pass

    return dir_name
