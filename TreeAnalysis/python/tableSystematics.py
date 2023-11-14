#!/usr/bin/env python

from __future__ import print_function
import json
import re
import pandas as pd
import logging

# Hierarchy: region | variable | sample | syst

def formatVariation(value):
    if   all(abs(v) < 1e-4 for v in [value['up'], value['dn']]):
        return '-'
    else:
        return '{:+2.2f}/{:+2.2f}'.format(value['up']*100, value['dn']*100)  # For slides or AN


def fillDataFrame(raw_data, formatter=formatVariation):
    # "Unpack" the inner dictionary {'up':x.xx, 'dn':x.xx} into a string
    data = { sample:
             { syst: formatter(value) for syst, value in d.items() }
             for sample, d in raw_data.items() }
    
    # Use pandas for pretty formatting
    return pd.DataFrame.from_dict(data) #, orient='index')


def tableRegion(df, samples=None, phstatus='default', sysdisplay=None, transpose=False, style='plain', **kwargs):

    # df = df.loc[label_order, :]
    # df = df.loc[:, sample_order]

    if(samples is not None):
        # If a list of samples was provided, use it
        df = df.reindex(samples, axis=1)
    else:
        # Otherwise, check if the prompt/nonpro must be kept (default), skipped, or used exclusively
        base_samples = {s.split('-prompt')[0].split('-nonpro')[0] for s in df.columns}
        if  (phstatus == 'skip'):
            df = df.loc[:, list(base_samples)]
        elif(phstatus == 'only'):
            phstat_samples = []
            for sample in base_samples:
                s_prompt = sample+'-prompt'
                s_nonpro = sample+'-nonpro'
                if  (s_prompt in df.columns):
                    phstat_samples.append(s_prompt)
                elif(s_nonpro in df.columns):
                    phstat_samples.append(s_nonpro)
                else:
                    phstat_samples.append(sample)
            df = df.loc[:, phstat_samples]

        df = df.reindex(sorted(df.columns), axis=1)

    # Only display certain systematics, if requested
    if(sysdisplay is not None):
        df = df.loc[sysdisplay,:]  # Reindex fills missing keys, loc raises a KeyError
    else:
        df = df.reindex(sorted(df.index  ), axis=0)

    # Transpose
    if(transpose):
        df = df.transpose()

    # Display
    if  (style == 'plain'):
        print(df.to_string())
    elif(style == 'latex'):
        print(df.to_latex())
        

def main():
    from argparse import ArgumentParser

    parser = ArgumentParser()
    parser.add_argument('-y', '--year', default=2016)
    parser.add_argument('-r', '--region'   , choices=['SR4P','CR3P1F','CR2P2F','SR3P'], default='SR4P')
    parser.add_argument('-i', '--inputfile', help='JSON file with systematics. Overrides year')
    parser.add_argument('-p', '--variable' , dest='variable', default='mZZGloose')
    parser.add_argument('-s', '--samples'  , nargs='+', help='Manually specify which samples to list')
    parser.add_argument('-S', '--systematics'  , dest='sysdisplay', nargs='+', help='Manually specify which systematics to display')
    parser.add_argument('-t', '--transpose'    , action='store_true', help='Transpose the table')
    parser.add_argument(      '--skip-phstatus', dest='phstatus', action='store_const', const='skip', help='Skip samples prompt/nonpro')
    parser.add_argument(      '--only-phstatus', dest='phstatus', action='store_const', const='only', help='Use only prompt-nonpro samples when available')
    parser.add_argument(      '--phstatus', choices=['skip', 'only', 'default'])
    parser.add_argument(      '--style'        , choices=('plain', 'latex'), default='plain')
    parser.add_argument(      '--log'          , dest='loglevel', metavar='LEVEL', default='WARNING', help='Level for the python logging module. Can be either a mnemonic string like DEBUG, INFO or WARNING or an integer (lower means more verbose).')
    args = parser.parse_args()
    loglevel = args.loglevel.upper() if not args.loglevel.isdigit() else int(args.loglevel)
    logging.basicConfig(format='%% %(levelname)s:%(module)s:%(funcName)s: %(message)s', level=loglevel)

    if(args.inputfile):
        sysFileName = args.inputfile
    else:
        sysFileName = 'data/systematics_{year}.json'.format(year=args.year)

    with open(sysFileName) as fin:
        systematics = json.load(fin)

    logging.info('Region: %6s, variable: %s', args.region, args.variable)
    if(not args.region in systematics.keys()):
        logging.error('region "{}" not found'.format(args.region))
        return

    df = fillDataFrame(systematics[args.region][args.variable])

    tableRegion(df, **vars(args))


if __name__ == '__main__':
    exit(main())
