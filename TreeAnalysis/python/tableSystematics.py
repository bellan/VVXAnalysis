#!/usr/bin/env python

from __future__ import print_function
import json
import pandas as pd

_padding_length = 80

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


def tableRegion(systematics, region, variables, **kwargs):
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

        # df = df.loc[label_order, :]
        # df = df.loc[:, sample_order]

        kwargs_samples = kwargs.get('samples')
        if(kwargs_samples is not None):
            df = df.reindex(kwargs_samples, axis=1)
        else:
            base_samples = {s.split('-prompt')[0].split('-nonpro')[0] for s in df.columns}
            if  (kwargs.get('phstatus') == 'skip'):
                df = df.loc[:, list(base_samples)]
            elif(kwargs.get('phstatus') == 'only'):
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

        df = df.reindex(sorted(df.index  ), axis=0)

        if(kwargs.get('transpose')):
            df = df.transpose()
        print(df.to_string())
        
        print('-'*_padding_length, end='\n\n')
    print('#'*_padding_length, end='\n\n')

if __name__ == '__main__':
    from argparse import ArgumentParser

    parser = ArgumentParser()
    parser.add_argument('-y', '--year', default=2016)
    parser.add_argument('-r', '--region'   , choices=['SR4P','CR3P1F','CR2P2F','SR3P'], default='SR4P')
    parser.add_argument('-i', '--inputfile', help='JSON file with systematics. Overrides year')
    parser.add_argument('-p', '--variables', dest='variables', nargs='+', default=['mZZ','mZZG','mWZ','mWZG'])
    parser.add_argument('-s', '--samples'  , nargs='+', help='Manually specify which samples to list')
    parser.add_argument(      '--transpose', action='store_true', help='Transpose the table')
    parser.add_argument(      '--skip-phstatus', dest='phstatus', action='store_const', const='skip', help='Skip samples prompt/nonpro')
    parser.add_argument(      '--only-phstatus', dest='phstatus', action='store_const', const='only', help='Use only prompt-nonpro samples when available')
    parser.add_argument(      '--phstatus', choices=['skip', 'only', 'default'])
    args = parser.parse_args()

    if(args.inputfile):
        sysFileName = args.inputfile
    else:
        sysFileName = 'data/systematics_{year}.json'.format(year=args.year)

    with open(sysFileName) as fin:
        systematics = json.load(fin)

    tableRegion(systematics, args.region, args.variables, transpose=args.transpose, phstatus=args.phstatus, samples=args.samples)
