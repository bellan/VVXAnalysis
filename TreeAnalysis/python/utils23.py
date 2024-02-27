#!/usr/bin/env python

############################################
# Utilties which do not require ROOT       #
#                                          #
# Author: A. Mecca (alberto.mecca@cern.ch) #
############################################

import sys
if sys.version_info.major <= 2:
    from collections import Mapping
else:
    from collections.abc import Mapping
from argparse import ArgumentParser

# Contains luminosity [pb^-1], error (as a lnN width suitable for Combine datacards)
lumi_dict = {
    '2016':        {'value': 36300, 'error_uncorrelated': 1.010, 'error_correlated':1.006, 'error_1718': 0    },
    '2016preVFP':  {'value': 19500, 'error_uncorrelated': 1.010, 'error_correlated':1.006, 'error_1718': 0    },
    '2016postVFP': {'value': 16800, 'error_uncorrelated': 1.010, 'error_correlated':1.006, 'error_1718': 0    },
    '2017':        {'value': 41480, 'error_uncorrelated': 1.020, 'error_correlated':1.009, 'error_1718': 1.006},
    '2018':        {'value': 59830, 'error_uncorrelated': 1.015, 'error_correlated':1.020, 'error_1718': 1.002},
    'Run2':        {'value':137620, 'error_uncorrelated': 1.0092,'error_correlated':1.013, 'error_1718': 1.0027}
}

def deep_update(orig, new):
    for k, v in new.items():
        if isinstance(v, Mapping):
            orig[k] = deep_update(orig.get(k, {}), v)
        else:
            orig[k] = v
    return orig


def _test_deep_update():
    d1 = {'A': {'a1': 1}, 'B': {'b1': 1}}
    d2 = {'A': {'a2': 2}, 'B': {'b1': 3, 'b2': 2}}
    target = {'A': {'a1': 1, 'a2':2}, 'B': {'b1':3, 'b2': 2}}
    assert deep_update(d1, d2) == target, 'deep_update test failed'

    d3 = {'A': [1,2,3]}
    d4 = {'A': [4,5]}
    target = {'A': [4,5]}
    assert deep_update(d3, d4) == target, 'deep_update failed to replace a list'


def common_parser(**kwargs):
    parser = ArgumentParser(**kwargs)
    parser.add_argument('-y', '--year'    , default='2018'       , choices=list(lumi_dict.keys()), help='Default: %(default)s')
    parser.add_argument('-i', '--inputdir', default='results'    , help='Base directory containing analyzer results (default: %(default)s)')
    parser.add_argument('-r', '--region'  , default='SR4P'       , help='Default: %(default)s')
    parser.add_argument('-A', '--analyzer', default='VVXAnalyzer', help='Name of the analyzer, used to compose the path of the input files (default: %(default)s)')
    parser.add_argument('--log', dest='loglevel', metavar='LEVEL', default='WARNING', help='Level for the python logging module. Can be either a mnemonic string like DEBUG, INFO or WARNING or an integer (lower means more verbose).')
    return parser


def config_logging(loglevel):
    import logging
    if(isinstance(loglevel, int)):
       level = loglevel
    elif(loglevel.isdigit()):
       level = int(loglevel)
    else:
       level = loglevel.upper()
    logging.basicConfig(format='%(levelname)s:%(module)s:%(funcName)s: %(message)s', level=level)


def main():
    '''
    Run the tests
    '''
    _test_deep_update()

    return 0


if __name__ == '__main__':
    exit(main())
