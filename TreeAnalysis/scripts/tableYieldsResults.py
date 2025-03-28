#!/usr/bin/env python

from argparse import ArgumentParser
from optparse import OptionParser
from math import sqrt
import logging
from ctypes import c_double
import re
import pandas as pd
import os

import ROOT

try:
    from VVXAnalysis.TreeAnalysis.plotUtils23 import InputDir, get_plots_singleyear
    from VVXAnalysis.TreeAnalysis.produceDataCard_VVGamma import get_strategy_config
    from VVXAnalysis.Combine.yieldutils import EventYield, print_yield
except ImportError:
    import sys
    sys.path.extend(['../Combine/python', 'python'])
    from plotUtils23 import InputDir, get_plots_singleyear
    from yieldutils import EventYield, print_yield
    from produceDataCard_VVGamma import get_strategy_config


SAMPLES_INFO = {
    'ggTo4e_Contin_MCFM701'   : {'kfactor': 1.7},
    'ggTo2e2mu_Contin_MCFM701': {'kfactor': 1.7},
    'ggTo4mu_Contin_MCFM701'  : {'kfactor': 1.7},
    'ZZTo4l'                  : {'kfactor': 1.325/1.256}
}

def main(args):
    logging.debug('args = %s', args)
    logging.warning('The k-factors are hardcoded in this scripts, instead of being retrieved from TreeAnalysis/python/samplesByRegion.py')

    # Read the config file to determine the list of processes
    if(args.config_file is not None):
        config = get_strategy_config(args.config_file)

        # Deduce the region from the config name
        if(args.region is None):
            args.region = list(config['regions'].keys())[0]
            logging.info('deduced the region from the config: %s', args.region)

        processes = list( config['regions'][args.region]['processes'].keys() )
    else:
        processes = ["ZZGTo4LG-prompt", "ZZGTo4LG-nonpro", "ZZTo4l-nonpro", "ggTo4e_Contin_MCFM701", "ggTo2e2mu_Contin_MCFM701", "ggTo4mu_Contin_MCFM701", "ZZZ", "WZZ", "WWZ", "TTZJets", "fake_leptons"]
        logging.warning('using the default list of processes: %s', processes)

    if(args.unblind):
        processes.append('data')

    logging.debug('args.region = %s', args.region)
    yields = get_yields(args.resultsdir, processes, args.hist_name, region=args.region, analyzer=args.analyzer, binname=args.binname)

    # logging.debug('yields: %s', yields)

    print_yield(yields, add_col_run2=True, add_row_scaled=True, float_format=args.format, **vars(args))

    return 0


def parse_args():
    parser = ArgumentParser(description='Reads one or more datacards, finds the histograms referenced by their shapeMaps and formats the expected and observed events in a LaTeX table')
    parser.add_argument('resultsdir', metavar='DIR')
    parser.add_argument('hist_name')
    parser.add_argument('-b', '--binname', help='Get the yield from a specific bin (only alphanumeric)')
    parser.add_argument('-c', '--config', dest='config_file', help='Configuration file that is used by produceDatacard_VVGamma (a JSON). It is used to get the processes and deduce the region.')
    parser.add_argument('-r', '--region')
    parser.add_argument('-A', '--analyzer', default='VVGammaAnalyzer', help='Default: %(default)s')
    parser.add_argument('-u', '--unblind', action='store_true', help='Print the "Observation" row')
    parser.add_argument(      '--format', '--float-format', default='%.2f', help='Format string used for floats (default: %(default)s)')
    parser.add_argument('--log', dest='loglevel', metavar='LEVEL', default='WARNING', help='Level for the python logging module. Can be either a mnemonic string like DEBUG, INFO or WARNING or an integer (lower means more verbose).')
    args = parser.parse_args()

    return args


def uniq(orig):
    '''
    Create a list of unique elements respecting the order in the original list
    '''
    out = []
    seen = set()
    for e in orig:
        if(e in seen): continue
        out.append(e)
        seen.update(e)
    return out


def sum_yields(*card_yields, **kwargs):
    '''
    Sums EventYields stored in dictionaries, for each process and bin
    '''
    out = {proc_name: {} for y in card_yields for proc_name in y.keys()}

    for proc_name in out.keys():
        out[proc_name] = {}
        bin_names = uniq(bin_name for card_yield in card_yields for bin_name in card_yield[proc_name].keys())

        for bin_name in bin_names:
            bin_tot_y = sum((card_yield[proc_name][bin_name] for card_yield in card_yields), EventYield())
            out[proc_name][bin_name] = bin_tot_y

    return out


def get_yields(resultsdir, samples, histogram, analyzer='VVGammaAnalyzer', region='SR4P', binname=None):
    '''
    Retrieve the histograms in the files pointed by the shapeMap in the card
    and use their integral and error to construct a dictionary that maps
    {process: {bin: [yield +- error]}}
    '''
    out = {}
    inputdir = InputDir(resultsdir, analyzer=analyzer, year=None, region=region)
    logging.debug('inputdir = %s', repr(inputdir))
    for sample in samples:
        if('-' in sample):
            # The sample is split in prompt/nonpro. We have to modify the plot name
            sample_name, prompt = sample.split('-')
            # and set a special sample name
            # The correct approach would be parsing variablesInfo. We just guess
            if(histogram.startswith('SYS')):
                split = histogram.split('_')
                split[1] += '-'+prompt
                hist_name = '_'.join(split)
            else:
                hist_name = histogram + '_' + prompt
        else:
            sample_name = sample
            hist_name = histogram

        kfactor = SAMPLES_INFO.get(sample_name, {}).get('kfactor', 1.)
        if(kfactor != 1.): logging.debug('kfactor for %15.15s: %f', sample_name, kfactor)
        out[sample] = {}

        for year in ('2016preVFP', '2016postVFP', '2017', '2018'):
            inputdir.year = year

            try:
                h = get_plots_singleyear(inputdir, sample_name, [hist_name])[0]
            except OSError:
                h = None
            if(h):
                if(binname is not None):
                    bin_i = h.GetXaxis().FindFixBin(binname)
                    integral = h.GetBinContent(bin_i)
                    error    = h.GetBinError(bin_i)
                else:
                    double_e = c_double(0.)
                    integral = h.IntegralAndError(0,-1, double_e)
                    error = double_e.value
            else:
                logging.warning('Could not get "%s" from %s', hist_name, os.path.join(inputdir.path(), sample_name+'.root'))
                integral = error = 0
            integral *= kfactor
            error    *= kfactor
            ey = EventYield(integral, error)
            logging.debug('event yield %15.15s (%s): %s', sample_name, year, ey)
            out[sample][year] = ey

    return out


if __name__ == '__main__':
    args = parse_args()
    loglevel = args.loglevel.upper() if not args.loglevel.isdigit() else int(args.loglevel)
    logging.basicConfig(format='%% %(levelname)s:%(module)s:%(funcName)s: %(message)s', level=loglevel)

    exit(main(args))
