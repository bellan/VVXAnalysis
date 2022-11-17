#!/usr/bin/env python2

from __future__ import print_function
import json
import pandas as pd
from argparse import ArgumentParser

parser = ArgumentParser()
parser.add_argument('-r', '--regions', dest='regions', nargs='+', choices=['SR4P','CR3P1F','CR2P2F','SR3P'], default=['SR4P','SR3P'])
parser.add_argument('-p', '--variables', dest='variables', nargs='+', default=['mZZ','mZZG','mWZ','mWZG'])
args = parser.parse_args()
print(args)

_padding_length = 80

with open('data/systematics.json') as fin:
    systematics = json.load(fin)

# Hierarchy: region | variable | sample | syst

def formatVariation(value):
    if   any(abs(v) > 0.7  for v in [value['up'], value['dn']]):
        return '*'
    if   all(abs(v) < 1e-4 for v in [value['up'], value['dn']]):
        return '-'
    else:
        return '{:+2.2f}/{:+2.2f}'.format(value['up']*100, value['dn']*100)  # For slides or AN
        # return '{:f}'.format(1+abs(value['up']-value['dn'])/2)               # For Combine datacard

for region in args.regions:
    if(not region in systematics.keys()):
        print('WARN: region "{}" not found'.format(region))
        continue
    nPadding_reg = (_padding_length - len(region))/2 - 1
    print('#'*nPadding_reg, region, '#'*nPadding_reg)
    
    for variable in args.variables:
        if(not variable in systematics[region].keys()):
            # print('WARN: variable "{}" not found in region "{}"'.format(variable, region))
            continue
        
        nPadding = (_padding_length - len(variable))/2 - 1
        print('-'*nPadding, variable, '-'*nPadding)

        raw_data = systematics[region][variable]
        # "Unpack" the inner dictionary {'up':x.xx, 'dn':x.xx} into a string
        data = { sample:
                 { syst: formatVariation(value) for syst, value in d.items() }
                 for sample, d in raw_data.items() }
        
        # Use pandas for pretty formatting
        df = pd.DataFrame.from_dict(data) #, orient='index')
        label_order = ['PDFVar', 'QCD-F', 'QCD-muR', 'alphas', 'L1Prefiring', 'puWeight', 'PhEScale', 'PhESigma', 'phEffSF', 'eleEffSF', 'muoEffSF', 'eleFakeRateSF', 'muoFakeRateSF']
        # if variable.startswith('mZZ'):
        #     df = df.loc[:, ['ZZGTo4LG', 'ZZTo4l', 'ggTo4l', 'fake_leptons']]
        
        print(df.loc[label_order, :].to_string())
        
        print('-'*_padding_length, end='\n\n')
    print('#'*_padding_length, end='\n\n')
