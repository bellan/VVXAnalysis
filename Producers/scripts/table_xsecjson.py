#!/usr/bin/env python3

from argparse import ArgumentParser
from os import environ, path
import pandas as pd
from packaging import version

parser = ArgumentParser(description='Visualize in a table the cross sections stored in local json')
parser.add_argument('-y', '--year', choices=['2016preVFP', '2016', '2017', '2018'], default=2017)
parser.add_argument(      '--path', help='Manually specify the path to the json')
args = parser.parse_args()

fname = None
if(args.path is not None):
    fname = args.path
else:
    fname = 'xsections{}.json'.format(args.year)
    if(environ.get('CMSSW_BASE')):
        fname = path.join(environ['CMSSW_BASE'], 'src', 'VVXAnalysis', 'Producers', 'scripts', fname)

df = pd.read_json(fname).transpose()

def my_format(l):
    if(not hasattr(l, '__iter__')):
        return l
    return [float('{:.5g}'.format(x)) for x in l]

df['mcm']    = df['mcm'].apply(my_format)
df['xsecdb'] = df['xsecdb'].apply(my_format)

pd_version = version.parse(pd.__version__)
opt = ''
if(pd_version <= version.parse('1.3')):
    opt = 'display.precision'
elif(pd_version <= version.parse('1.5')):
    opt = 'format.precision'
with pd.option_context(opt, 4):
    print(df[['identifier', 'mcm', 'xsecdb', 'internal', 'external', 'csv']])
