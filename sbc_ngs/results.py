'''
sbc-ngs (c) University of Manchester 2019

All rights reserved.

@author: neilswainston
'''
# pylint: disable=invalid-name
import os
from threading import Thread

import numpy as np
import pandas as pd


class ResultsThread(Thread):
    '''Thread-safe class to write results to DataFrames.'''

    def __init__(self, columns, barcodes_df, queue):
        self.__dfs = {}

        self.__dfs['raw_summary'] = barcodes_df
        self.__dfs['raw_summary'].set_index(
            ['forward', 'reverse', 'barcode_type'], inplace=True)

        # Initialise specific Dataframes:
        for name in ['mutations', 'nucleotides', 'indels', 'deletions',
                     'identity', 'depths']:
            self.__dfs[name] = self.__init_df(columns)

        self.__queue = queue
        self.__closed = False
        Thread.__init__(self)

    def run(self):
        '''Run.'''
        while True:
            task = self.__queue.get()

            if task is None:
                break

            self.__update_df(task)
            self.__queue.task_done()

    def write(self, dir_name):
        '''Write result output files.'''
        self.__update_raw_summary()

        self.__add_summary()

        for name, df in self.__dfs.items():
            df.to_csv(os.path.join(dir_name, name + '.csv'))

    def close(self):
        '''Close.'''
        self.__queue.put(None)

    def __init_df(self, columns):
        '''Initialise a results dataframe.'''
        df = pd.concat([pd.DataFrame(columns=columns),
                        self.__dfs['raw_summary'].copy()],
                       sort=False)
        df.index = self.__dfs['raw_summary'].index
        return df

    def __update_df(self, task):
        '''Update DataFrame.'''
        name, val, col_id, row_ids = task
        self.__dfs[name][col_id].loc[row_ids[0], row_ids[1], row_ids[2]] = val

    def __update_raw_summary(self):
        '''Update raw summary.'''
        self.__dfs['identity'].fillna(0, inplace=True)

        numerical_df = self.__dfs['identity'].copy()

        if 'plate_idx' in self.__dfs['identity'].columns:
            numerical_df = numerical_df.drop(['plate_idx'], axis=1)

        numerical_df = numerical_df.select_dtypes(include=[np.float])

        if not numerical_df.empty:
            self.__dfs['raw_summary']['matched_seq_id'] = \
                numerical_df.idxmax(axis=1)
            self.__dfs['raw_summary']['identity'] = numerical_df.max(axis=1)

            for name, df in self.__dfs.items():
                if name != 'raw_summary':
                    self.__dfs['raw_summary'][name] = \
                        df.lookup(df.index,
                                  self.__dfs['raw_summary']['matched_seq_id'])

            # Remove spurious unidentified entries:
            self.__dfs['raw_summary'] = \
                self.__dfs['raw_summary'][self.__dfs['raw_summary']
                                          ['identity'] != 0]

            # Sort:
            self.__dfs['raw_summary'] = \
                self.__dfs['raw_summary'].sort_values(
                    'identity', ascending=False)

    def __add_summary(self):
        '''Add summary.'''
        summary_df = self.__dfs['raw_summary'][
            ['well', 'mutations', 'deletions', 'depths']]

        summary_df = summary_df.iloc[
            summary_df.index.get_level_values('barcode_type') == 'sum']

        summary_df.index = summary_df.index.droplevel('barcode_type')

        summary_df['mutations'] = summary_df['mutations'].apply(
            lambda x: x if x else np.nan)

        summary_df['deletions'] = summary_df['deletions'].apply(
            lambda x: x if x else np.nan)

        self.__dfs['summary'] = summary_df
