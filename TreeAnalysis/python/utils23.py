#!/usr/bin/env python

############################################
# Utilties which do not require ROOT       #
#                                          #
# Author: A. Mecca (alberto.mecca@cern.ch) #
############################################

import sys
import os
from errno import EEXIST
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

NONPROMPT_LOW = 'nonprompt' # 'non-prompt' # American/British spelling
NONPROMPT_CAP = NONPROMPT_LOW.capitalize()


def deep_update(orig, new):
    for k, v in new.items():
        if isinstance(v, Mapping):
            orig[k] = deep_update(orig.get(k, {}), v)
        else:
            orig[k] = v
    return orig


def get_VVXAnalysis(default=None):
    if('CMSSW_BASE' in os.environ):
        return os.path.join(os.environ['CMSSW_BASE'], 'src', 'VVXAnalysis')
    elif('VVXAnalysis' in os.getcwd()):
        return os.path.join(os.getcwd().split('VVXAnalysis')[0], 'VVXAnalysis')
    else:
        return default


# Emulate os.makekdirs(..., exists_ok=True) for python2
if(sys.version_info.major < 3):
    def makedirs_ok(*args, **kwargs):
        try: os.makedirs(*args, **kwargs)
        except OSError as e:
            if(e.errno != EEXIST): raise e  # Catch only "File esists"
else:
    def makedirs_ok(*args, **kwargs):
        os.makedirs(*args, exist_ok=True, **kwargs)


def byteify(data, ignore_dicts = False):
    if isinstance(data, str):
        return data

    # If this is a list of values, return list of byteified values
    if isinstance(data, list):
        return [ byteify(item, ignore_dicts=True) for item in data ]
    # If this is a dictionary, return dictionary of byteified keys and values
    # but only if we haven't already byteified it
    if isinstance(data, dict) and not ignore_dicts:
        return {
            byteify(key, ignore_dicts=True): byteify(value, ignore_dicts=True)
            for key, value in data.items() # changed to .items() for Python 2.7/3
        }

    # Python 3 compatible duck-typing
    # If this is a Unicode string, return its string representation
    if str(type(data)) == "<type 'unicode'>":
        return data.encode('utf-8')

    # If it's anything else, return it in its original form
    return data


def _test_deep_update():
    d1 = {'A': {'a1': 1}, 'B': {'b1': 1}}
    d2 = {'A': {'a2': 2}, 'B': {'b1': 3, 'b2': 2}}
    target = {'A': {'a1': 1, 'a2':2}, 'B': {'b1':3, 'b2': 2}}
    assert deep_update(d1, d2) == target, 'deep_update test failed'

    d3 = {'A': [1,2,3]}
    d4 = {'A': [4,5]}
    target = {'A': [4,5]}
    assert deep_update(d3, d4) == target, 'deep_update failed to replace a list'


def _test_byteify():
    source = {u'a': [1,2,3], 'b': {u'b1': 2, u'b2': u'2'}, u'c': {u'c1': {u'c11': u'12', u'c12': [u'1', u'2']}}}
    target = { 'a': [1,2,3], 'b': { 'b1': 2,  'b2':  '2'},  'c': { 'c1': { 'c11':  '12',  'c12': [ '1',  '2']}}}
    assert byteify(source) == target, 'byteify failed to convert a dictionary as expected'


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
    _test_byteify()

    return 0


if __name__ == '__main__':
    exit(main())
