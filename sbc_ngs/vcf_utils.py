'''
sbc-ngs (c) University of Manchester 2019

All rights reserved.

@author: neilswainston
'''
# pylint: disable=chained-comparison
# pylint: disable=consider-using-dict-comprehension
# pylint: disable=invalid-name
# pylint: disable=too-many-arguments
# pylint: disable=too-many-locals
import operator
import os
import re
import sys

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

import numpy as np
import pandas as pd


def analyse(vcf_filename, target_id, src_id, write_queue):
    '''Analyse a given vcf file.'''
    num_matches, mutations, indels, deletions, templ_len, \
        consensus_seq, depths = analyse_vcf(vcf_filename)

    consensus_filename = vcf_filename.replace('.vcf', '_consensus.fasta')
    record = SeqRecord(Seq(consensus_seq), id=vcf_filename)
    SeqIO.write([record], consensus_filename, 'fasta')

    write_queue.put(
        ['identity', num_matches / float(templ_len), target_id, src_id])
    write_queue.put(
        ['mutations', mutations, target_id, src_id])
    write_queue.put(
        ['indels', indels, target_id, src_id])
    write_queue.put(
        ['deletions', deletions, target_id, src_id])
    write_queue.put(
        ['depths', max(depths) if depths else 0, target_id, src_id])


def vcf_to_df(vcf_filename):
    '''Convert vcf to Pandas dataframe.'''
    if not os.path.exists(vcf_filename):
        return None, 0

    data = []
    templ_len = float('NaN')

    with open(vcf_filename) as vcf:
        for line in vcf:
            if line.startswith('##'):
                mtch = re.match(r'(?:.*),length=(.*)>', line)

                if mtch:
                    templ_len = int(mtch.group(1))
                else:
                    pass
            elif line.startswith('#'):
                columns = line[1:].split()[:-1] + ['DATA']
            else:
                data.append(line.split())

    df = _expand_info(pd.DataFrame(columns=columns, data=data))

    df['POS'] = df['POS'].astype(int)

    if 'DP' in df.columns:
        df['DP'] = df['DP'].astype(int)
        df['DP_PROP'] = df['DP'] / df['DP'].max()

    if 'INDEL' not in df.columns:
        df['INDEL'] = False

    df['INDEL'] = df['INDEL'].fillna(value=False)

    return df, templ_len


def analyse_vcf(vcf_filename, dp_filter=0.0, qs_threshold=0.0):
    '''Analyse vcf file, returning number of matches, mutations and
    indels.'''
    num_matches = 0
    mutations = []
    indels = []
    deletions = []
    consensus_seq = []
    depths = []

    df, templ_len = vcf_to_df(vcf_filename)

    for _, row in df.iterrows():
        if (dp_filter > 1 and row['DP'] > dp_filter) \
                or row['DP_PROP'] > dp_filter:
            # Extract QS values and order to find most-likely base:
            # Compare most-likely term to reference:
            max_gs = max(_get_qs(row).items(), key=operator.itemgetter(1))

            if row.get('INDEL', False):
                if row['REF'] != max_gs[0]:
                    indels.append((row['REF'] + str(row['POS']) + max_gs[0],
                                   max_gs[1]))
            else:
                consensus_seq.append(max_gs[0])

                if row['REF'] != max_gs[0] and max_gs[1] > qs_threshold:
                    consensus_seq.append(max_gs[0])

                    mutations.append((row['REF'] + str(row['POS']) + max_gs[0],
                                      max_gs[1]))
                else:
                    consensus_seq.append(row['REF'])
                    num_matches += 1

            depths.append(row['DP'])
        else:
            deletions.append(row['POS'])

    return num_matches, mutations, indels, _get_ranges_str(deletions), \
        templ_len, ''.join(consensus_seq), depths


def analyse_dir(parent_dir):
    '''Analyse directory.'''
    for fle in os.listdir(parent_dir):
        if os.path.isdir(fle):
            for subfle in fle:
                forward_df, _ = vcf_to_df(os.path.join(subfle, 'forward.vcf'))
                reverse_df, _ = vcf_to_df(os.path.join(subfle, 'reverse.vcf'))
                forward_qs = _get_qs_by_pos(forward_df) if forward_df else None
                reverse_qs = _get_qs_by_pos(reverse_df) if reverse_df else None


def _expand_info(df):
    '''Expand out INFO column from vcf file.'''
    infos = []

    for row in df.itertuples():
        info = [term.split('=') for term in row.INFO.split(';')]

        infos.append({term[0]: (term[1] if len(term) == 2 else True)
                      for term in info})

    return df.join(pd.DataFrame(infos, index=df.index))


def _get_ranges_str(vals):
    '''Convert list of integers to range strings.'''
    return ['-'.join([str(r) for r in rnge])
            if rnge[0] != rnge[1]
            else rnge[0]
            for rnge in _get_ranges(vals)]


def _get_ranges(vals):
    '''Convert list of integer to ranges.'''
    ranges = []

    if vals:
        sorted_vals = sorted(vals)

        i = 0

        for j in range(1, len(sorted_vals)):
            if sorted_vals[j] > 1 + sorted_vals[j - 1]:
                ranges.append((sorted_vals[i], sorted_vals[j - 1]))
                i = j

        ranges.append((sorted_vals[i], sorted_vals[-1]))

    return ranges


def _get_qs_by_pos(df):
    '''Get QS values per position.'''
    return {row['POS']: _get_qs(row) for _, row in df.iterrows()}


def _get_qs(row):
    '''Get QS for a given positional row.'''
    nucl_qs = {n: 0.0 for n in 'ACGT'}

    nucl = [row['REF']] + [nucl for nucl in row['ALT'].split(',')
                           if nucl != '<*>']

    qs = [float(val) for val in dict([term.split('=')
                                      for term in row['INFO'].split(';')
                                      if '=' in term])
          ['QS'].split(',')
          if float(val) > 0]

    nucl_qs.update(dict(zip(nucl, qs)))

    return nucl_qs


def _get_probs(forward, reverse, num_forw, num_rev, prior_nucl,
               prior_prob=0.97):
    '''Get probabilities.'''
    prior = {nucl: (1 - prior_prob) / 3 for nucl in 'ACGT'}
    prior[prior_nucl] = prior_prob

    consistency = _get_consistency(forward, reverse)

    return (consistency / (num_forw + num_rev) *
            (num_forw * forward + num_rev * reverse)) + \
        (1 - consistency) * np.fromiter(prior.values(), dtype=float)


def _get_consistency(forward, reverse):
    '''Get consistency.'''
    minima = np.minimum(forward, reverse)
    return np.true_divide(np.sum(minima), np.sum(reverse))


def main(args):
    '''main method.'''
    analyse_dir(args[0])


if __name__ == '__main__':
    main(sys.argv[1:])
