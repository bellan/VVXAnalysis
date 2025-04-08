#!/usr/bin/env python3

################################################################
# Utilities to print LaTeX tables with the yield of processes  #
# of interest in a given region ("bin" in Combine parlance)    #
#                                                              #
# Author: A. Mecca (alberto.mecca@cern.ch) - 2024              #
################################################################

import logging
from math import sqrt, isnan
import sys
import pandas as pd

try:
    from VVXAnalysis.TreeAnalysis.utils23 import lumi_dict
except ImportError:
    sys.path.append('../TreeAnalysis/python')
    from utils23 import lumi_dict

pandas_version = tuple(int(i) for i in pd.__version__.split('.'))
pandas_brackets = pandas_version[0] >= 2

class EventYield:
    '''
    A class that handles correctly the statistical error on the event yield (sum in quadrature)
    '''
    def __init__(self, val=0, err=0):
        self.val = val
        self.err = err

    def __add__(self, other):
        new = EventYield(self.val, self.err)
        new += other
        return new

    def __iadd__(self, other):
        self.val += other.val
        self.err = sqrt(self.err**2 + other.err**2)
        return self

    def __repr__(self):
        return '(' + str(self.val) + '+-' + str(self.err) + ')'

    def to_string(self, fmt='%.4g'):
        return (fmt+' \\pm '+fmt) %(self.val, self.err)

    def __iter__(self):
        return (i for i in (self.val, self.err))


def print_yield(data, unblind=False, float_format='%.4g', add_col_run2=False, add_row_scaled=False, **kwargs):
    '''
    Prints to stdout a LaTeX table of event yields using Pandas DataFrame's to_latex()
    '''
    # Create list of bins that contains unique elements but mantains the order in which they were in data
    # this could not be archieved with a list comprehension (no uniqueness)
    # nor with sorting a set (the starting order would not be respeted)
    # bin_names = []
    # for _, proc_data in data.items():
    #     for bin_name in proc_data:
    #         if bin_name not in bin_names:
    #             bin_names.append(bin_name)

    # Remove data_obs
    series_data = data.pop('data_obs', None)

    df = pd.DataFrame(data)
    # df.fillna(value=EventYield(0,0), inplace=True)
    # On lxplus, pandas version is 1.2.2, and fillna is bugged (does not accept object values)
    df = df.applymap(lambda x: x if not(isinstance(x, float) and isnan(x)) else EventYield(0,0))

    # Transpose: rows = samples, columns = years
    df = df.transpose()

    # Sort
    def sort_func(row):
        '''sort by sample (signal first), then by yield'''
        index = row['index']
        if('ZZGTo4LG'   in index): k0 = 2
        if('WZGTo3LNuG' in index): k0 = 1
        else:                      k0 = 0

        k1 = (row['2016preVFP']+row['2016postVFP']+row['2017']+row['2018']).val
        return (k0, k1)

    df['index'] = df.index
    df['sort_key'] = df.apply(sort_func, axis=1)
    df.sort_values('sort_key', ascending=False, inplace=True)
    df.drop(columns=['sort_key', 'index'], inplace=True)

    # Total yield of MC for each year
    df.loc['Total'] = df.sum()

    if unblind:
        if(series_data is not None):
            logging.debug('data series: %s', series_data)
            df.loc['Data'] = series_data
        else:
            logging.warning('data series is None')

    # Base columns and their format. Additional columns may be added
    column_format='l >{$}r<{$} >{$}r<{$} >{$}r<{$} >{$}r<{$}'
    yhcell = r'\yhcell{{%s}}' if pandas_brackets else r'\yhcell{%s}'
    header = ['{}'] + [yhcell%(y) for y in ('2016preVFP', '2016postVFP', '2017', '2018')]

    # Compute total yield for Run2
    if(add_col_run2):
        df[r'Run2'] = df.sum(axis=1)
        column_format += ' >{$}r<{$}'
        header.append(yhcell %(r'\Run2'))

    # ARC request (Toni): total yield scaled to Run2
    if(add_row_scaled):
        df.loc[r'Tot. (scaled to \Run2)'] = get_scaled_row(df.loc['Total'])

    # Convert to string using the format supplied by command line args
    formatters = [lambda x:x.to_string(fmt=float_format)]*len(df.columns)

    # Add a column with the same content as the index (sample names)
    df.insert(loc=0, column='process',value=df.index)
    # Column formatters
    formatters = [sample_to_latex] + formatters

    df_string  = df.to_latex(formatters=formatters
                             , escape=False
                             , column_format=column_format
                             , header=header
                             , index=False)

    # Add a small space before the Total row
    out_split = []
    for line in df_string.split('\n'):
        if(line.strip().startswith('Total')):
            out_split.append(r'\noalign{\vspace{.3ex}}\hline\noalign{\vspace{.3ex}}')
        out_split.append(line)
    out_string = '\n'.join(out_split)

    sys.stdout.write(out_string)


def sample_to_latex(sample):
    '''Used to re-index with latex expressions instead of sample names'''
    sample.replace('_Contin_MCFM701', '')
    if('-' in sample):
        base, extra = sample.split('-')
    else:
        base = sample
        extra = None
    if  ('ZZGTo4LG'       in base): base = r'$\PZ\PZ\PGg\to4\Pl\PGg$'
    elif('WZGTo3LNuG'     in base): base = r'$\PW\PZ\PGg\to3\Pl\PGn\PGg$'
    elif('WZTo3LNu'       in base): base = r'$\PW\PZ\to3\Pl\PGn$'
    elif('ZZTo4l-nonpro'== sample): return r'\qqZZnonpro'
    elif('ZZTo4l'         in base): base = r'\qqZZ'
    elif('ggTo4mu'        in base): base = r'\ggtomm'
    elif('ggTo2e2mu'      in base): base = r'\ggtoem'
    elif('ggTo4e'         in base): base = r'\ggtoee'
    elif('ggTo4l'         in base): base = r'\ggtoll'
    elif('ZZZ'            in base): base = r'$\PZ\PZ\PZ$'
    elif('WZZ'            in base): base = r'$\PW\PZ\PZ$'
    elif('WWZ'            in base): base = r'$\PW\PW\PZ$'
    elif('TTZJets'        in base): base = r'$\PQt\PAQt\PZ$+jets'
    elif('ZGToLLG'        in base): base = r'$\PZ\PGg\to\Pl\Pl$'
    elif('TZq'            in base): base = r'$\PQt\PZ\PQq$'
    elif('tW'             in base): base = r'$\PQt\PW$'
    elif('DYJetsToLL_M50' in base): base = r'\DYnonpro'
    elif('fake_leptons'   in base): base = r'Nonprompt leptons'
    elif('fake_photons'   in base): base = r'Nonprompt photons'

    if extra is None:
        return base
    elif('ZZGTo4LG' in sample):
        if(extra == 'nonpro'):
            # Out-of-Acceptance
            return base+' OOA'
        else:
            return base
    else:
        return '-'.join([base, extra])


def get_scaled_row(total_series):
    '''
    Return a dict with the ratio of the yield of each year (scaled to Run2 lumi), to the total Run2 yield
    '''
    out = pd.Series()
    yield_t = total_series['Run2']
    lumi_t = lumi_dict['Run2']['value']
    logging.debug('%-12s: %.4g (%.4g pb-1)', 'Run2', yield_t.val, lumi_t)
    for year, yield_y in total_series.items(): #('2016preVFP', '2016postVFP', '2017', '2018'):
        lumi_y = lumi_dict[year]['value']
        r_y  = yield_y.val/yield_t.val * lumi_t/lumi_y
        e_y = sqrt(yield_y.err**2/yield_t.val**2 + yield_y.val**2/yield_t.val**4 * yield_y.err**2) * lumi_t/lumi_y
        logging.debug('%-12s: %.4g (%.4g pb-1) --> %.4g', year, yield_y.val, lumi_y, r_y)
        out[year] = EventYield(r_y, e_y)
    return out
