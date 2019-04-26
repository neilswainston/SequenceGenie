'''
sbc-ngs (c) University of Manchester 2019

All rights reserved.

@author: neilswainston
'''
# pylint: disable=too-many-arguments
# pylint: disable=too-many-locals
import os
import sys
from synbiochem.utils import ice_utils, seq_utils

from sbc_ngs.pathway import PathwayAligner


def get_ice_files(url, username, password, ice_ids_filename,
                  for_primer, rev_primer, out_dir):
    '''Get ICE sequences.'''
    ice_client = ice_utils.ICEClient(url, username, password)

    with open(ice_ids_filename, 'r') as ice_ids_file:
        ice_ids = [line.strip() for line in ice_ids_file]

    seqs_offsets = [seq_utils.pcr(ice_client.get_ice_entry(ice_id).get_seq(),
                                  for_primer, rev_primer)
                    for ice_id in ice_ids]

    seqs, offsets = zip(*seqs_offsets)

    ice_files = {ice_id:
                 (seq_utils.write_fasta({ice_id: seq},
                                        os.path.join(out_dir,
                                                     ice_id + '.fasta')),
                  len(seq))
                 for ice_id, seq in zip(ice_ids, seqs)}

    pcr_offsets = {ice_id: offset for ice_id, offset in zip(ice_ids, offsets)}

    # Get Genbank files for subsequent data analysis:
    for ice_id in ice_ids:
        gb_filename = os.path.join(out_dir, ice_id + '.gb')
        ice_client.get_genbank(ice_id, gb_filename)

    return ice_files, pcr_offsets, [len(seq) for seq in seqs]


def main(args):
    '''main method.'''
    out_dir = args[0]
    in_dir = args[1]

    # Get pathway sequences from ICE:
    seq_files, _, _ = \
        get_ice_files(url=args[2],
                      username=args[3],
                      password=args[4],
                      ice_ids_filename=os.path.join(in_dir, 'ice_ids.txt'),
                      for_primer=args[5],
                      rev_primer=args[6],
                      out_dir=out_dir)

    aligner = PathwayAligner(out_dir,
                             in_dir,
                             seq_files,
                             min_length=int(args[7]),
                             max_read_files=int(args[8]))

    aligner.score_alignments(int(args[9]), num_threads=int(args[10]))


if __name__ == '__main__':
    main(sys.argv[1:])
