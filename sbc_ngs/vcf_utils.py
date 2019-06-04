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
from collections import defaultdict
import operator
import os
import re
import sys

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from scipy.stats import binom

import numpy as np
import pandas as pd


def analyse(vcf_filename, target_id, src_id, write_queue):
    '''Analyse a given vcf file.'''
    num_matches, mutations, nucleotides, indels, deletions, templ_len, \
        consensus_seq, depths = analyse_vcf(vcf_filename)

    consensus_filename = vcf_filename.replace('.vcf', '_consensus.fasta')
    record = SeqRecord(Seq(consensus_seq), id=vcf_filename)
    SeqIO.write([record], consensus_filename, 'fasta')

    write_queue.put(
        ['identity', num_matches / float(templ_len), target_id, src_id])
    write_queue.put(
        ['mutations', mutations, target_id, src_id])
    write_queue.put(
        ['nucleotides', nucleotides, target_id, src_id])
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


def analyse_vcf(vcf_filename, dp_filter=0.0):
    '''Analyse vcf file, returning number of matches, mutations and indels.'''
    print(vcf_filename)
    num_matches = 0
    muts = []
    nucls = []
    indels = []
    deletions = []
    consensus_seq = []
    depths = []

    df, templ_len = vcf_to_df(vcf_filename)

    dir_name, name = os.path.split(vcf_filename)
    str_specs = {}

    if name == 'sum.vcf':
        str_specs = get_strand_specific(dir_name).get(dir_name, {})

    for _, row in df.iterrows():
        if (dp_filter > 1 and row['DP'] > dp_filter) \
                or row['DP_PROP'] > dp_filter:
            # Extract QS values and order to find most-likely base:
            # Compare most-likely term to reference:
            qs, _, _ = _get_qs(row)
            max_gs = max(qs.items(), key=operator.itemgetter(1))

            if row.get('INDEL', False):
                if row['REF'] != max_gs[0]:
                    indels.append((row['REF'] + str(row['POS']) + max_gs[0],
                                   max_gs[1]))
            else:
                consensus_res = str_specs.get(row['POS'], max_gs)
                consensus_seq.append(consensus_res[0])
                match = True

                if row['REF'] != max_gs[0]:
                    match = False
                    nucls.append((row['REF'] + str(row['POS']) + max_gs[0],
                                  max_gs[1]))

                if name == 'sum.vcf':
                    str_spec = str_specs.get(row['POS'], None)

                    if str_spec and row['REF'] != str_spec[0]:
                        match = False
                        muts.append((row['REF'] + str(row['POS']) +
                                     str_spec[0], str_spec[1]))
                    else:
                        match = True

                if match:
                    num_matches += 1

            depths.append(row['DP'])
        else:
            deletions.append(row['POS'])

    return num_matches, muts, nucls, indels, _get_ranges_str(deletions), \
        templ_len, ''.join(consensus_seq), depths


def get_strand_specific(parent_dir):
    '''Get strand specific scores.'''
    strand_specific = defaultdict(lambda: defaultdict(list))
    vcf_files = defaultdict(list)

    for root, _, files in os.walk(parent_dir, topdown=False):
        for name in files:
            if re.search(r'^(forward|reverse)\.vcf$', name):
                vcf_files[root].append(os.path.join(root, name))

    for root, files in vcf_files.items():
        qs_df = _get_qs_df(files)

        for pos, row in qs_df.iterrows():
            missing_0 = row.iloc[qs_df.columns.get_level_values(0) == 0].empty
            missing_1 = row.iloc[qs_df.columns.get_level_values(0) == 1].empty
            ref = row[1, 'REF'] if missing_0 else row[0, 'REF']

            probs = _get_probs(np.array([0.25] * 4)
                               if missing_0
                               else np.array(list(row[0, 'nucls'].values())),
                               np.array([0.25] * 4)
                               if missing_1
                               else np.array(list(row[1, 'nucls'].values())),
                               0 if missing_0 else row[0, 'DP'],
                               0 if missing_1 else row[1, 'DP'],
                               ref)

            if ref != 'ACGT'[np.argmax(probs)]:
                strand_specific[root][pos] = (['ACGT'[np.argmax(probs)],
                                               np.max(probs),
                                               probs,
                                               ref])

    return {key: dict(val) for key, val in strand_specific.items()}


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


def _get_qs_df(filenames):
    '''Get QS values per position.'''
    qs_dfs = []

    for filename in filenames:
        # name = os.path.splitext(os.path.basename(filename))[0]
        df, _ = vcf_to_df(filename)

        if df is not None and not df.empty:
            df = df.set_index('POS')
            df = df[~df['INFO'].str.contains('INDEL')]
            # df = df[~df['INFO'].str.contains('DP=0')]
            qs = df.apply(_get_qs, axis=1)
            qs_df = pd.DataFrame(list(qs),
                                 index=qs.index,
                                 columns=['nucls', 'DP', 'REF']).dropna(axis=1)
            qs_df.columns = pd.MultiIndex.from_tuples(
                [(len(qs_dfs), col) for col in qs_df.columns])
            qs_dfs.append(qs_df)

    if len(qs_dfs) > 1:
        return pd.merge(qs_dfs[0], qs_dfs[1],
                        left_index=True, right_index=True)

    return qs_dfs[0] if qs_dfs else pd.DataFrame()


def _get_qs(row):
    '''Get QS for a given positional row.'''
    nucl_qs = {n: 0.0 for n in 'ACGT'}

    nucl = [row['REF']] + [nucl for nucl in row['ALT'].split(',')
                           if nucl != '<*>']

    info = dict([term.split('=')
                 for term in row['INFO'].split(';')
                 if '=' in term])

    qs = [float(val) for val in info['QS'].split(',')
          if float(val) > 0]

    nucl_qs.update(dict(zip(nucl, qs)))

    return (nucl_qs if sum(nucl_qs.values()) else {n: 0.25 for n in 'ACGT'},
            int(info['DP']), row['REF'])


def _get_probs(forward, reverse, num_forw, num_rev, prior_nucl):
    '''Get probabilities.'''
    mutprob = 0.0005

    prior = [1 - mutprob if nucl == prior_nucl else mutprob / 3
             for nucl in 'ACGT']

    prior_nucl2 = np.array([0.25, 0.25, 0.25, 0.25])

    nf = [round(p * num_forw) for p in forward]
    nr = [round(p * num_rev) for p in reverse]
    misprob_f = (sum(nf) - max(nf)) / \
        max(sum(nf), 1e-16)  # misread probability
    misprob_r = (sum(nr) - max(nr)) / \
        max(sum(nr), 1e-16)  # misread probability

    like_f = binom.pmf(nf, num_forw, 1 - misprob_f)
    like_r = binom.pmf(nr, num_rev, 1 - misprob_r)

    prob_nucl1 = prior_nucl2 * like_r * like_f  # Probs reads are informative
    prob_nucl1 = prob_nucl1 / sum(prob_nucl1)

    # Estimate of discordance between forward and reverse:
    concord = sum(like_f * like_r) / (sum(like_f) * sum(like_r))

    return concord * prob_nucl1 + (1 - concord) * np.array(prior)


def _get_probs_old(forward, reverse, num_forw, num_rev, prior_nucl,
                   prior_prob=0.95):
    '''Get probabilities.'''
    prior = {nucl: (1 - prior_prob) / 3 for nucl in 'ACGT'}
    prior[prior_nucl] = prior_prob

    intersection = _get_intersection(forward, reverse)

    if num_forw + num_rev:
        return (intersection / (num_forw + num_rev) *
                (num_forw * forward + num_rev * reverse)) + \
            (1 - intersection) * np.fromiter(prior.values(), dtype=float)

    return np.array(list(prior.values()))


def _get_intersection(forward, reverse):
    '''Get intersection.'''
    minima = np.minimum(forward, reverse)
    return np.true_divide(np.sum(minima), np.sum(reverse))


def main(args):
    '''main method.'''
    import multiprocessing as mp
    from sbc_ngs import demultiplex, results

    in_dir = args[0]
    results_dir = args[1]
    out_dir = args[2]

    _, barcodes_df = \
        demultiplex.get_barcodes(os.path.join(in_dir, 'barcodes.csv'))

    barcodes_df.rename(columns={'actual_ice_id': 'known_seq_id'},
                       inplace=True)
    seq_ids = barcodes_df['known_seq_id'].unique()

    write_queue = mp.Manager().Queue()
    results_thread = results.ResultsThread(sorted(seq_ids),
                                           barcodes_df,
                                           write_queue)
    results_thread.start()

    for dirpath, _, filenames in os.walk(os.path.abspath(results_dir)):
        for filename in filenames:
            filepath, filename = os.path.split(os.path.join(dirpath, filename))

            if filename.endswith('.vcf'):
                dirs = filepath.split(os.sep)

                analyse(os.path.join(dirpath, filename),
                        dirs[-1],
                        dirs[-2].split('_') + [filename.replace('.vcf', '')],
                        write_queue)

    # Update summary:
    results_thread.close()
    results_thread.write(out_dir)


if __name__ == '__main__':
    np.seterr(all='raise')
    main(sys.argv[1:])
