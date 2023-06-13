#!/usr/bin/env python

from __future__ import print_function
import json
import pandas as pd

_padding_length = 80

# Hierarchy: region | variable | sample | syst

def formatVariation(value):
    if   any(abs(v) > 0.7  for v in [value['up'], value['dn']]):
        return '*'
    if   all(abs(v) < 1e-4 for v in [value['up'], value['dn']]):
        return '-'
    else:
        return '{:+2.2f}/{:+2.2f}'.format(value['up']*100, value['dn']*100)  # For slides or AN
        # return '{:f}'.format(1+abs(value['up']-value['dn'])/2)               # For Combine datacard


def fillDataFrame(raw_data, formatter=formatVariation):
    # "Unpack" the inner dictionary {'up':x.xx, 'dn':x.xx} into a string
    data = { sample:
             { syst: formatter(value) for syst, value in d.items() }
             for sample, d in raw_data.items() }
    
    # Use pandas for pretty formatting
    return pd.DataFrame.from_dict(data) #, orient='index')


def tableRegion(systematics, region, variables):
    if(not region in systematics.keys()):
        print('WARN: region "{}" not found'.format(region))
        return
    nPadding_reg = int((_padding_length - len(region))/2) - 1
    print('#'*nPadding_reg, region, '#'*nPadding_reg)
    
    for variable in variables:
        if(not variable in systematics[region].keys()):
            # print('WARN: variable "{}" not found in region "{}"'.format(variable, region))
            continue
        
        nPadding = int((_padding_length - len(variable))/2) - 1
        print('-'*nPadding, variable, '-'*nPadding)

        df = fillDataFrame(systematics[region][variable])
        
        # label_order = ['PDFVar', 'QCD-F', 'QCD-muR', 'alphas', 'L1Prefiring', 'puWeight', 'PhEScale', 'PhESigma', 'phEffSF', 'eleEffSF', 'muoEffSF', 'eleFakeRateSF', 'muoFakeRateSF']
        # df = df.loc[label_order, :]
        
        print(df.to_string())
        
        print('-'*_padding_length, end='\n\n')
    print('#'*_padding_length, end='\n\n')

if __name__ == '__main__':
    from argparse import ArgumentParser

    parser = ArgumentParser()
    parser.add_argument('-y', '--year', default=2016)
    parser.add_argument('-r', '--regions', dest='regions', nargs='+', choices=['SR4P','CR3P1F','CR2P2F','SR3P'], default=['SR4P','SR3P'])
    parser.add_argument('-i', '--inputfile', help='JSON file with systematics. Overrides year')
    parser.add_argument('-p', '--variables', dest='variables', nargs='+', default=['mZZ','mZZG','mWZ','mWZG'])
    args = parser.parse_args()
    # print(args)

    if(args.inputfile):
        sysFileName = args.inputfile
    else:
        sysFileName = 'data/systematics_{year}.json'.format(year=args.year)

    with open(sysFileName) as fin:
        systematics = json.load(fin)
    
    for region in args.regions:
        tableRegion(systematics, region, args.variables)
