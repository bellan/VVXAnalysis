#!/usr/bin/env python

from __future__ import print_function
import os
import sys
import json
import re
import string
from math import isnan
import pandas as pd
from argparse import ArgumentParser
from samplesByRegion import getSamplesByRegion
from tableSystematics import fillDataFrame

parser = ArgumentParser()
parser.add_argument('template', help='Template for the datacard')
parser.add_argument('-v', '--verbose', action='count', default=0, help='Increase the verbosity level')
parser.add_argument('-r', '--region', default='SR4P')
parser.add_argument('-y', '--year', type=int, default=2016)
parser.add_argument('-C', '--config-file', help='Configuration file. Defaults to one hardcoded in this script')
parser.add_argument('-c', '--config', type=json.loads, help='String convertible to dictionary used to override the config', default={})
parser.add_argument(      '--unblind', action='store_true')
parser.add_argument(      '--path', default='/afs/cern.ch/user/a/amecca/public/histogramsForCombine')

args = parser.parse_args()
if(args.verbose):
    print('INFO: writing card for {args.year}, {args.region}'.format(*globals()))

### Hardcoded configuration ###
config = {
    # Define which samples are signals and which are background
    'signals': ['ZZGTo4LG', 'WZGTo3LNuG'],
    'backgrounds': ['ZZTo4l', 'ggTo4l', 'WZTo3LNu'],
    # Define the observable and the samples in each region
    'regions': {
        'SR4P': {
            'processes': {'ZZGTo4LG': -1, 'ZZTo4l':-1, 'ggTo4l':-1, 'WZTo3LNu':-1},
            'observable': {'name':'mZZGloose', 'observation':-1}  # Combine's name for "observable"
        },
        'CR3P1F':{
            'processes': {'ZZTo4l':-1, 'ggTo4l':-1, 'WZTo3LNu':-1},
            'observable': {'name':'mZZ', 'observation':-1}
        },
        'SR3P': {
            'processes': {'WZGTo3LNuG':-1, 'ZZGTo4LG':-1, 'ZZTo4l':-1, 'ggTo4l':-1, 'WZTo3LNu':-1},
            'observable': {'name':'mWZG', 'observation':-1}
        }
    },

    # General configuration
    'systematics':{
        'lumi': {
            2016: 1.012,
            2017: 1.023,
            2018: 1.025,
            1618: 1.016
        },
        'has_shape': ['QCDscaleF', 'QCD-F']
    }
}

# Update from config file
if(args.config_file is not None):
    with open(args.config_file) as f:
        fconfig = json.load(f)
    config.update(fconfig)

# Update from command line
config.update(args.config)


# Utility functions
def getSystType(syst):
    if(any(re.search(expr, syst) for expr in config['systematics']['has_shape'])):
        return 'shape'
    else:
        return 'lnN'

def formatForCombine(value):
    if   any(isnan(v)      for v in [value['up'], value['dn']]):
        return '-'
    if   all(abs(v) < 1e-4 for v in [value['up'], value['dn']]):
        # print('WARN: negligible syst:', value)
        return '-'
    if   any(abs(v) > 1    for v in [value['up'], value['dn']]):
        print('WARN: very large syst:', value)
        return '-'
    else:
        return '{:f}'.format(1+abs(value['up']-value['dn'])/2)


def getBinName(region, observable):
    return observable  # region+'_'+observable


# Formatter class that ignores missing arguments
try:
    # Python 3
    from _string import formatter_field_name_split
except ImportError:
    formatter_field_name_split = str._formatter_field_name_split

class PartialFormatter(string.Formatter):
    def get_field(self, field_name, args, kwargs):
        try:
            val = super(PartialFormatter, self).get_field(field_name, args, kwargs)
        except (IndexError, KeyError, AttributeError):
            first, _ = formatter_field_name_split(field_name)
            val = '{' + field_name + '}', first
        return val


### Bin section ###
region_config = config['regions'][args.region]

observables=[[
    getBinName(args.region, region_config['observable']['name']),
    region_config['observable']['observation']
]]
df_bin = pd.DataFrame(observables, columns=['bin', 'observation']).transpose()

if(args.verbose > 1): print(df_bin.to_string(header=False))

### Observable section ###
signals     = [proc for proc in region_config['processes'] if proc in config['signals']    ]
backgrounds = [proc for proc in region_config['processes'] if proc in config['backgrounds']]

if(args.verbose):
    print('>>> {signals=}'    .format(*globals))
    print('>>> {backgrounds=}'.format(*globals))
    print('>>> {observables=}'.format(*globals))

minProc = 1 - len(signals)
samples_to_idx      = {s:minProc+i for i,s in enumerate(signals)    }
samples_to_idx.update({b:i+1       for i,b in enumerate(backgrounds)})
if(args.verbose):
    print('>>> {samples_to_idx=}'.format(*globals))

df_rate = pd.DataFrame(
    [[getBinName(args.region, region_config['observable']['name']), k, samples_to_idx[k], v] for k,v in region_config['processes'].items()],
    columns=['bin', 'process', 'process', 'rate']
).transpose()
if(args.verbose > 1):
    print(df_rate.to_string(header=False))
    print()


##### SYSTEMATICS #####
### Read systematics ###
sysFile = 'data/systematics_{year}.json'.format(year=args.year)
systematics = {}
if(os.path.exists(sysFile)):
    with open(sysFile) as f:
        systematics = json.load(f)
else:
    print('ERROR: "{}" does not exist. Try (re-)running systematics.py to generate it'.format(sysFile))
    exit(2)


### Process systematics ###
# MEMO: setdefault(region, {}).setdefault(var, {}).setdefault(sample, {})[syst] = {'up':upVar, 'dn':dnVar}

# Order the systematics so that the samples have the same order of the observable section
data_syst = { sample:systematics[args.region][region_config['observable']['name']][sample] for sample in samples_to_idx }
df_syst = fillDataFrame(data_syst, formatter=formatForCombine)

df_syst = df_syst.rename(lambda x: 'CMS_'+x)

lumi = config['systematics']['lumi'][args.year]
df_syst.loc['CMS_lumi_13TeV'] = pd.Series({ sample: lumi for sample in df_syst.columns })
df_syst.insert(0, 'type', [ getSystType(syst) for syst in df_syst.index ], True)

if(args.verbose > 1):
    print(df_syst)
    print()


### Open template ###
with open(args.template) as ftemplate:
    template = ftemplate.read()
    
# template = template.format(
template = PartialFormatter().format(template,
    imax=1,
    jmax=len(signals)+len(backgrounds)-1,
    kmax=len(df_syst),
    path=os.path.join(args.path, str(args.year), args.region+'.root'),
    bins=df_bin.to_string(header=False),
    processes=df_rate.to_string(header=False),
    systematics='#'+df_syst.to_string()
)

### Write card ###
with open('combine/{}_{}.txt'.format(args.year, args.region), 'w') as fout:
    fout.write(template)
    print('INFO: config written to "{}"'.format(fout.name))