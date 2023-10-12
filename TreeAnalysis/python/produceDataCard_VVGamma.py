#!/usr/bin/env python

from __future__ import print_function
import os
import sys
import json
import re
import string
import copy
from math import isnan
import pandas as pd
from argparse import ArgumentParser
from samplesByRegion import getSamplesByRegion
from tableSystematics import fillDataFrame
import logging

### Hardcoded configuration ###
__builtin_config__ = {
    # Define which samples are signals and which are background
    'signals': ['ZZGTo4LG', 'WZGTo3LNuG'],
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

__builtin_template__ = '''\
imax {imax} number of channels
jmax {jmax} number of backgrounds
kmax {kmax} number of nuisance parameters
------------

shapes * * {path} $CHANNEL/$PROCESS $CHANNEL/$PROCESS_$SYSTEMATIC
------------

{bins}
------------

{processes}
------------

{systematics}
'''


# Utility functions
def getSystType(syst, config):
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


def main():
    parser = ArgumentParser()
    parser.add_argument('config_file', help='Configuration file. Defaults to the one hardcoded in this script', default=None)
    parser.add_argument('-t', '--template', help='Template for the datacard')
    parser.add_argument('-v', '--verbose'  , dest='verbosity', action='count', default=1, help='Increase verbosity')
    parser.add_argument(      '--verbosity', type=int, help='Set verbosity')
    parser.add_argument('-q', '--quiet'    , dest='verbosity', action='store_const', const=0, help='Set verbose to minimum')
    parser.add_argument('-r', '--region', default='SR4P')
    parser.add_argument('-y', '--year', type=int, default=2016)
    parser.add_argument('-c', '--config', type=json.loads, help='String convertible to dictionary used to override the config', default={})
    parser.add_argument(      '--unblind', action='store_true')
    parser.add_argument(      '--path', default='/afs/cern.ch/user/a/amecca/public/histogramsForCombine', help='Path to the histograms')
    parser.add_argument('--log', dest='loglevel', metavar='LEVEL', default='WARNING', help='Level for the python logging module. Can be either a mnemonic string like DEBUG, INFO or WARNING or an integer (lower means more verbose).')

    args = parser.parse_args()
    loglevel = args.loglevel.upper() if not args.loglevel.isdigit() else int(args.loglevel)
    logging.basicConfig(format='%(levelname)s:%(module)s:%(funcName)s: %(message)s', level=loglevel)

    logging.info('writing card for %(year)s, %(region)s', vars(args))


    config = copy.deepcopy(__builtin_config__)

    # Update from config file
    if(args.config_file is not None):
        try:
            with open(args.config_file) as f:
                fconfig = json.load(f)
        except json.decoder.JSONDecodeError as e:
            print('ERROR: Caught', type(e), 'while reading', args.config_file)
            print(e)
            return 1
        config.update(fconfig)

    # Update from command line
    config.update(args.config)

    if(args.verbosity >= 3):
        print('### Configuration ###')
        print(json.dumps(config, indent=2))
        print()


    ### Bin section ###
    region_config = config['regions'][args.region]

    observables=[[
        getBinName(args.region, region_config['observable']['name']),
        region_config['observable']['observation']
    ]]
    df_bin = pd.DataFrame(observables, columns=['bin', 'observation']).transpose()

    if(args.verbosity >= 3):
        print(df_bin.to_string(header=False))
        print()

    ### Observable section ###
    signals     = [proc for proc in region_config['processes'] if proc in config['signals']    ]
    backgrounds = [proc for proc in region_config['processes'] if proc not in config['signals']]

    logging.info('signals:     %s', signals    )
    logging.info('backgrounds: %s', backgrounds)
    logging.info('observables: %s', observables)

    minProc = 1 - len(signals)
    samples_to_idx      = {s:minProc+i for i,s in enumerate(signals)    }
    samples_to_idx.update({b:i+1       for i,b in enumerate(backgrounds)})

    logging.info('samples_to_idx: %s', samples_to_idx)

    df_rate = pd.DataFrame(
        [[getBinName(args.region, region_config['observable']['name']), k, samples_to_idx[k], v] for k,v in region_config['processes'].items()],
        columns=['bin', 'process', 'process_number', 'rate']
    ).transpose()
    df_rate.sort_values('process_number', axis=1, inplace=True)
    df_rate.rename(index={'process_number':'process'}, inplace=True)

    if(args.verbosity >= 3):
        print(df_rate.to_string(header=False))
        print()


    ##### SYSTEMATICS #####
    ### Read systematics ###
    sysFile = 'data/systematics_{year}.json'.format(year=args.year)
    logging.debug('sysFile = %s', sysFile)
    systematics = {}
    if(os.path.exists(sysFile)):
        with open(sysFile) as f:
            systematics = json.load(f)
    else:
        logging.critical('"%s" does not exist. Try (re-)running systematics.py to generate it', sysFile)
        return 1


    ### Process systematics ###
    # MEMO: setdefault(region, {}).setdefault(var, {}).setdefault(sample, {})[syst] = {'up':upVar, 'dn':dnVar}

    # Order the systematics so that the samples have the same order of the observable section
    missing_systematics = False
    data_syst = {}
    for sample in samples_to_idx:
        try:
            data_syst[sample] = systematics[args.region][region_config['observable']['name']][sample]
        except KeyError:
            logging.error('Systematics empty for %s', sample)
            data_syst[sample] = {}
            missing_systematics = True
    df_syst = fillDataFrame(data_syst, formatter=formatForCombine).fillna(0)

    df_syst = df_syst.rename(lambda x: 'CMS_'+x)

    lumi = config['systematics']['lumi'][args.year]
    df_syst.loc['CMS_lumi_13TeV'] = pd.Series({ sample: lumi for sample in df_syst.columns })
    df_syst.insert(0, 'type', [ getSystType(syst, config) for syst in df_syst.index ], True)

    if(args.verbosity >= 4):
        print(df_syst)
        print()


    ### Open template ###
    if(args.template is not None):
        with open(args.template) as ftemplate:
            template = ftemplate.read()
    else:
        template = __builtin_template__

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
    cardname = os.path.join('combine', '{}_{}.txt'.format(args.year,# args.region,
                                                             args.config_file.split('/')[-1].split('.')[0] if args.config_file is not None else 'default'
                                                             ))
    with open(cardname, 'w') as fout:
        if(missing_systematics): fout.write('### WARNING systematics missing for one or more samples ###\n\n')
        fout.write(template)
        logging.info('config written to "%s"', fout.name)

    if(missing_systematics):
        return 2
    return 0

if __name__ == '__main__':
    exit(main())
