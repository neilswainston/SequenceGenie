'''
sbc-ngs (c) University of Manchester 2019

All rights reserved.

@author: neilswainston
'''
# pylint: disable=too-many-arguments
# pylint: disable=too-many-locals
import os

from synbiochem.utils import ice_utils, seq_utils


def get_ice_files(url, username, password, ice_ids_filename,
                  for_primer, rev_primer, dir_name):
    '''Get ICE sequences.'''
    ice_client = ice_utils.ICEClient(url, username, password)

    with open(ice_ids_filename, 'rU') as ice_ids_file:
        ice_ids = [line.strip() for line in ice_ids_file]

    seqs_offsets = [seq_utils.pcr(ice_client.get_ice_entry(ice_id).get_seq(),
                                  for_primer, rev_primer)
                    for ice_id in ice_ids]

    seqs, offsets = zip(*seqs_offsets)

    ice_files = {ice_id:
                 (seq_utils.write_fasta({ice_id: seq},
                                        os.path.join(dir_name,
                                                     ice_id + '.fasta')),
                  len(seq))
                 for ice_id, seq in zip(ice_ids, seqs)}

    pcr_offsets = {ice_id: offset for ice_id, offset in zip(ice_ids, offsets)}

    # Get Genbank files for subsequent data analysis:
    for ice_id in ice_ids:
        gb_filename = os.path.join(dir_name, ice_id + '.gb')
        ice_client.get_genbank(ice_id, gb_filename)

    return ice_files, pcr_offsets, [len(seq) for seq in seqs]
