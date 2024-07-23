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
from plotUtils23 import TFileContext
from plotUtils import makedirs_ok
from utils23 import lumi_dict, byteify
import logging

### Hardcoded configuration ###
__builtin_config__ = {
    # Define which samples are signals and which are background
    'signals': ['ZZGTo4LG', 'WZGTo3LNuG'],
    # Define the observable and the samples in each region
    'regions': {
        'SR4P': {
            'processes': {'ZZGTo4LG': -1, 'ZZTo4l':-1, 'WZTo3LNu':-1},
            'observable': {'name':'mZZGloose', 'observation':-1}  # Combine's name for "observable"
        },
        'CR3P1F':{
            'processes': {'ZZTo4l':-1, 'WZTo3LNu':-1},
            'observable': {'name':'mZZ', 'observation':-1}
        },
        'SR3P': {
            'processes': {'WZGTo3LNuG':-1, 'ZZGTo4LG':-1, 'ZZTo4l':-1, 'WZTo3LNu':-1},
            'observable': {'name':'mWZG', 'observation':-1}
        }
    },

    # General configuration
    'systematics':{
        'shape': [],
        'correlated'  : ['L1Prefiring', 'PDFVar', 'QCDscale', 'alphas', 'phEffSF', 'phEffMVASF', 'phEScale', 'phESigma', 'muoEffSF', 'eleEffSF', 'puWeight', 'phFakeRate'],
        'uncorrelated': ['electronVeto', 'muoFakeRateSF', 'eleFakeRateSF', 'fake_leptons_norm', 'fake_photons_norm', 'phARstat'],
        'skip-if-signal': ['PDFVar', 'QCDscale', 'alphas'],
        'theory': ['QCDscale', 'alphas', 'PDFVar'],
        'datadriven': ['phFakeRate', 'phARstat', 'muoFakeRateSF', 'eleFakeRateSF'],
        '_end':[]
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

{groups}
'''


# Utility functions
def getSystType(syst, config):
    if  (syst in config['systematics'].get('shape', [])):
        return 'shape'
    elif(syst in config['systematics'].get('gmN'  , [])):
        return 'gmN'
    else:
        return 'lnN'

def format_lnN(value):
    # up = 1 - yield_up/yield --> k_up = 1 + up
    up = value['up']
    dn = value['dn']
    if   any(isnan(v)      for v in [up, dn]):
        return '-'
    if   all(abs(v) < 1e-4 for v in [up, dn]):
        # print('WARN: negligible syst:', value)
        return '-'
    if   any(abs(v) > 1    for v in [up, dn]):
        logging.warning('very large syst: %s', value)
        return '-'
    else:
        symmetric = abs(up-dn)/2
        asymmetry = abs(up+dn)/2  # In case of symmetric effect up and dn have opposite sign
        if  (symmetric == 0):
            return '-'
        elif( up*dn > 0 or asymmetry > 0.01 ):
            return '{:f}/{:f}'.format(1+dn, 1+up)
        else:
            return '{:f}'.format(1 + symmetric)


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

def check_exisiting_histograms(fname, observable, processes):
    '''
    Filters the list "processes" and returns only those that actually have an histogram in the files for the observable
    '''
    logging.debug('fname     : %s', fname     )
    logging.debug('observable: %s', observable)
    logging.debug('processes : %s', processes )
    existing_processes = {}
    with TFileContext(fname) as tf:
        obs_folder = tf.Get(observable)
        if(not obs_folder):
            raise KeyError('Observable "%s" missing from file "%s"' %(observable, fname))
        # logging.debug('obs_folder: %s', obs_folder)
        keys = obs_folder.GetListOfKeys()
        # logging.debug('keys: %s', '\n\t'+'\n\t'.join(sorted([k.GetName() for k in keys if not 'CMS' in k.GetName()])))
        for process in processes:
            if  (not keys.Contains(process)):
                logging.warning('dropping %s, since it is missing for observable "%s" in %s', process, observable, fname)
            else:
                h = obs_folder.Get(process)
                integral = h.Integral(1, h.GetNbinsX())
                if(integral <= 0):
                    logging.warning('dropping %s, since it has norm %.3g for observable "%s" in %s', process, integral, observable, fname)
                else:
                    existing_processes[process] = integral
    return existing_processes

def get_gmN_params(syst, data_syst):
    logging.debug('syst: %s', syst)
    n_affected = 0
    sample_affected = None
    N = 0
    alpha = 0
    for sample, sample_data in data_syst.items():
        syst_data = sample_data[syst]
        if(syst_data['up'] - syst_data['dn'] > 0):
            n_affected += 1
            sample_affected = sample
            N      = syst_data['N']
            integr = syst_data['integr']
            alpha  = integr/N
            logging.debug('sample: %s - N: %d - integr: %g - alpha: %g', sample, N, integr, alpha)

    if  (n_affected >  1):
        raise RuntimeError('The number of samples affected by "%s" is %d, but the specified type is gmN!' %(syst, n_affected))
    elif(n_affected == 0):
        logging.warning('No sample is affected by "%s", but the specified type is "gmN"' %(syst))

    return sample_affected, N, alpha

def get_gmN_params_local(fname, observable, process):
    '''
    Retrieve the gmN parameters N and alpha from the local rootfile (which has an histogramsForCombine layout)
    Parameters:
    :param fname: path to file with histograms
    :param observable: name TDirectory in the file
    :param process: name of the process within the TDirectory
    '''
    with TFileContext(fname) as tf:
        h = tf.Get(observable).Get(process)
        if(not h): raise KeyError('observable "%s", process "%s" in file "%s"' %(observable, process, fname))
        N      = h.GetEntries()
        integr = h.Integral(1, h.GetNbinsX())
        alpha  = integr/N
        logging.debug('process: %s - N: %d - integr: %g - alpha: %g', process, N, integr, alpha)
    return N, alpha

def get_shape_affected(syst, data_syst):
    samples_affected = []
    for sample, sample_data in data_syst.items():
        syst_data = sample_data[syst]
        if(syst_data['up'] - syst_data['dn'] != 0.):
            samples_affected.append(sample)

    if(len(samples_affected) == 0):
        logging.warning('No sample is affected by "%s", but the specified type is "shape"' %(syst))

    logging.debug('syst: %-12s - affected(%d): %s', syst, len(samples_affected), samples_affected)
    return samples_affected

def main():
    parser = ArgumentParser()
    parser.add_argument('config_file', help='Configuration file')
    parser.add_argument('-t', '--template', help='Template for the datacard')
    parser.add_argument('-v', '--verbose'  , dest='verbosity', action='count', default=1, help='Increase verbosity')
    parser.add_argument(      '--verbosity', type=int, help='Set verbosity')
    parser.add_argument('-q', '--quiet'    , dest='verbosity', action='store_const', const=0, help='Set verbose to minimum')
    parser.add_argument('-r', '--region', default=None)
    parser.add_argument('-y', '--year'  , default='2018', choices=lumi_dict.keys())
    parser.add_argument('-c', '--config', type=json.loads, help='String convertible to dictionary used to override the config', default={})
    parser.add_argument(      '--unblind', action='store_true')
    parser.add_argument(      '--path' , default=None                                    , help='Path to the histograms (default: %(default)s)')
    parser.add_argument('-i', '--input', default='histogramsForCombine', dest='localpath', help='Path to the histograms for local checks (default: %(default)s)')
    parser.add_argument('--log', dest='loglevel', metavar='LEVEL', default='WARNING', help='Level for the python logging module. Can be either a mnemonic string like DEBUG, INFO or WARNING or an integer (lower means more verbose).')

    args = parser.parse_args()
    loglevel = args.loglevel.upper() if not args.loglevel.isdigit() else int(args.loglevel)
    logging.basicConfig(format='%(levelname)s:%(module)s:%(funcName)s: %(message)s', level=loglevel)

    logging.info('writing card for %(year)s, %(region)s', vars(args))

    if(args.path is None):
        # This is a temporary hack until I understand why combineCards.py mishandles relative paths
        args.path = os.path.join('/afs/cern.ch/work/a/amecca/Analysis/Combine/CMSSW_11_3_4/src/VVXAnalysis/Combine/test', args.localpath)

    config = copy.deepcopy(__builtin_config__)

    # Update from config file
    try:
        with open(args.config_file) as f:
            fconfig = json.load(f, object_hook=byteify)
    except json.decoder.JSONDecodeError as e:
        print('ERROR: Caught', type(e), 'while reading', args.config_file)
        print(e)
        return 1
    config.update(fconfig)
    # Update 'systematics' with more granularity
    config['systematics'] = copy.deepcopy(__builtin_config__['systematics'])
    config['systematics'].update(fconfig.get('systematics', {}))

    # Update from command line
    config.update(args.config)

    if(args.verbosity >= 3):
        print('### Configuration ###')
        print(json.dumps(config, indent=2))
        print()

    if args.region is None:
        if(len(config['regions']) == 1):
            args.region = list(config['regions'].keys())[0]
            logging.info('Using region %s from the config file', args.region)
        else:
            raise RuntimeError('No region specified in arguments. The config defines: %s' %(config['regions'].keys()) )

    # Path
    path_to_histograms = os.path.join(args.path, str(args.year), args.region+'.root')
    path_to_histograms_local = os.path.join(args.localpath, str(args.year), args.region+'.root')

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

    # Test that all the processes have at least one event - otherwise Combine will crash later
    existing_processes = check_exisiting_histograms(path_to_histograms_local, observables[0][0], region_config['processes'])
    region_config_processes = {}
    for proc_name, proc_yield_config in region_config['processes'].items():
        if proc_name not in existing_processes:
            continue
        # Yield: prefer what is specified in the config, resort to the integral in the (local) rootfile
        y = existing_processes[proc_name] if proc_yield_config == -1 else proc_yield_config
        region_config_processes[proc_name] = y
    region_config['processes'] = region_config_processes

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
            logging.error('Systematics from JSON empty for %s', sample)
            data_syst[sample] = {}
            missing_systematics = True

        # Zero theoretical uncertainty on signal cross section, since we are measuring it
        if(sample in config['signals']):
            for syst, val in data_syst[sample].items():
                if(syst in config['systematics']['skip-if-signal']):
                    val['dn'] = val['up'] = 0
                    logging.debug('Zeroed systematic "%s" for sample "%s"', syst, sample)

    # Set normalization uncertainty (e.g. fake_leptons and fake_photons)
    for sample, val in config['systematics'].get('norm_uncertainty', {}).items():
        if sample in data_syst:
            logging.info('setting norm uncertainty on %s (%s)', sample, val)
            data_syst[sample][sample+'_norm'] = val
            if  (sample == 'fake_leptons'):
                logging.info('Using norm uncertainty instead of lepton fake rate uncertainty')
                data_syst[sample]['eleFakeRateSF'] = {'up':1, 'dn':1}
                data_syst[sample]['muoFakeRateSF'] = {'up':1, 'dn':1}

    # Set normalization for groups of samples
    for group_name, group_info in config['systematics'].get('norm_group_uncertainty', {}).items():
        logging.info('setting group norm uncertainty on %s (%s)', group_name, group_info['value'])
        for sample in group_info['samples']:
            data_syst[sample][group_name+'_norm'] = group_info['value']

    df_syst = fillDataFrame(data_syst, formatter=format_lnN).fillna(0)
    type_column = []
    for syst in df_syst.index:
        # Drop systematics that do not affect any sample
        if(all(df_syst[column].loc[syst] == '-' for column in df_syst.columns)):
            logging.info('dropping systematic "%s"', syst)
            df_syst.drop(syst, inplace=True)
            continue

        # Decide systematic type and append it to column
        syst_type = getSystType(syst, config)
        if(syst_type == 'gmN'):
            sample_affected, N, alpha = get_gmN_params(syst, data_syst)
            if(N > 0):
                df_syst[sample_affected].loc[syst] = alpha
                type_column.append('gmN %d' %(N))
            else:
                type_column.append('lnN')
        elif(syst_type == 'shape'):
            for sample in get_shape_affected(syst, data_syst):
                df_syst[sample].loc[syst] = 1
            type_column.append('shape')
        else:
            type_column.append(syst_type)

    # Additional systematics
    for syst, process in config['systematics'].get('ADD_gmN', {}).items():
        logging.debug('additional syst: %s -> %s', syst, process)
        try:
            N, alpha = get_gmN_params_local(path_to_histograms_local, observables[0][0], process)
        except KeyError as e:
            logging.error("%s --> skipping this gmN in datacard", e)
            continue

        if(N > 0):
            new_row = [ alpha if p == process else '-' for p in df_syst.columns ]
            df_syst.loc[syst] = new_row
            type_column.append('gmN %d' %(N))
        else:
            logging.error('N=%d for gmN syst %s', N, syst)
            return 1

    def rename_syst(syst):
        # Uses implicitly: config, year
        if  (syst in config['systematics'][  'correlated']): suffix = ''
        elif(syst in config['systematics']['uncorrelated']): suffix = '_'+args.year
        elif(syst.endswith('_norm')):
            logging.info('Treating %s as correlated among years', syst)
            suffix = ''
        else: raise RuntimeError('Systematic "%s": unspecified if correlated' %(syst))
        return 'CMS_{}{}'.format(syst, suffix)

    df_syst = df_syst.rename(rename_syst)

    lumi_uncorrelated = lumi_dict[args.year]['error_uncorrelated']
    lumi_correlated   = lumi_dict[args.year]['error_correlated']
    lumi_1718         = lumi_dict[args.year]['error_1718']

    df_syst.loc['CMS_lumi_13TeV_%s'%(args.year)] = pd.Series({ sample: lumi_uncorrelated for sample in df_syst.columns })
    type_column.append('lnN')
    df_syst.loc['CMS_lumi_13TeV_correlated']     = pd.Series({ sample: lumi_correlated   for sample in df_syst.columns })
    type_column.append('lnN')
    if(args.year in ('2017', '2018')):
        df_syst.loc['CMS_lumi_13TeV_1718']       = pd.Series({ sample: lumi_1718         for sample in df_syst.columns })
        type_column.append('lnN')

    df_syst.insert(0, 'type', type_column, False)

    if(args.verbosity >= 4):
        print(df_syst)
        print()

    ### Groups of nusiances ###
    groups = {}
    for syst_fullname in df_syst.index:
        is_lumi   = 'lumi' in syst_fullname
        is_theory = any(syst in syst_fullname for syst in config['systematics']['theory'])
        is_datadr = any(syst in syst_fullname for syst in config['systematics']['datadriven'])
        logging.debug('syst_fullname: %s, is_theory: %d, is_datadr: %d, is_lumi: %d', syst_fullname, is_theory, is_datadr, is_lumi)
        if  (is_lumi):
            groups.setdefault('lumi'  , []).append(syst_fullname)
        elif(is_theory):
            groups.setdefault('theory', []).append(syst_fullname)
        elif(is_datadr):
            groups.setdefault('datadr', []).append(syst_fullname)
    groups_s = '\n'.join( ['%s group = %s'%(k, ' '.join(v)) for k,v in groups.items()] )


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
        path=path_to_histograms,
        bins=df_bin.to_string(header=False),
        processes=df_rate.to_string(header=False),
        systematics='#'+df_syst.to_string(),
        groups=groups_s
    )

    ### Write card ###
    cardname = os.path.join('combine/cards',
                            '{}_{}.txt'.format(args.year,# args.region,
                                                             args.config_file.split('/')[-1].split('.')[0] if args.config_file is not None else 'default'
                                                             ))
    makedirs_ok(os.path.dirname(cardname))  # make directories for card if they do not exist
    with open(cardname, 'w') as fout:
        if(missing_systematics): fout.write('### WARNING systematics missing for one or more samples ###\n\n')
        fout.write(template)
        logging.info('config written to "%s"', fout.name)

    if(missing_systematics):
        return 2
    return 0

if __name__ == '__main__':
    exit(main())
