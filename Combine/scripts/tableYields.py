#!/usr/bin/env python3

from argparse import ArgumentParser
from optparse import OptionParser
from fnmatch import fnmatch
from string import Template
from math import sqrt, isnan
import logging
from ctypes import c_double
import re
import pandas as pd

import ROOT

from HiggsAnalysis.CombinedLimit import DatacardParser
# from VVXAnalysis.TreeAnalysis.utils23 import deep_update


class EventYield:
    '''
    A class that handles correctly the statistical error on the event yield (sum in quadrature)
    '''
    def __init__(self, val=0, err=0):
        self.val = val
        self.err = err

    def __add__(self, other):
        new = EventYield(self.val, self.err)
        new += other
        return new

    def __iadd__(self, other):
        self.val += other.val
        self.err = sqrt(self.err**2 + other.err**2)
        return self

    def __repr__(self):
        return '(' + str(self.val) + '+-' + str(self.err) + ')'

    def to_string(self, fmt='%.4g'):
        return (fmt+' pm '+fmt) %(self.val, self.err)

    def __iter__(self):
        return (i for i in (self.val, self.err))


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


def parse_args():
    parser = ArgumentParser(description='Reads one or more datacards, finds the histograms referenced by their shapeMaps and formats the expected and observed events in a LaTeX table')
    parser.add_argument('datacards', nargs='+')
    parser.add_argument('-u', '--unblind', action='store_true', help='print the "Observation" row')
    parser.add_argument(      '--float-format', default='%.4g', help='Format string used for floats (default: %(default)s)')
    parser.add_argument('--log', dest='loglevel', metavar='LEVEL', default='WARNING', help='Level for the python logging module. Can be either a mnemonic string like DEBUG, INFO or WARNING or an integer (lower means more verbose).')
    args, unknown = parser.parse_known_args()

    parser_datacard = OptionParser()
    DatacardParser.addDatacardParserOptions(parser_datacard)
    args_datacard = parser_datacard.parse_args(unknown)
    return args, args_datacard


def find_in_shapeMap(string, shapeMap):
    '''
    Tries to find a matching pattern in a shapeMap and returns the object associated to it
    (usually another shapeMap (which is really just a dictionary) or a list
    '''
    for pattern, data in shapeMap.items():
        if fnmatch(string, pattern):
            return data
    raise LookupError('No pattern matching "%s" in shapeMap. Keys are: %s' %(string, shapeMap.keys()))


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


def print_yield(data, unblind=False, float_format='%.4g', **kwargs):
    '''
    Prints to stdout a LaTeX table of event yields using Pandas DataFrame's to_latex()
    '''
    # Create list of bins that contains unique elements but mantains the order in which they were in data
    # this could not be archieved with a list comprehension (no uniqueness)
    # nor with sorting a set (the starting order would not be respeted)
    # bin_names = []
    # for _, proc_data in data.items():
    #     for bin_name in proc_data:
    #         if bin_name not in bin_names:
    #             bin_names.append(bin_name)

    # Compute total yield
    # total = {bin_name: [0, 0] for bin_name in bin_names}
    # for proc_name, proc_data in data.items():
    #     for bin_name, proc_bin_data in proc_data.items():
    #         total[bin_name] += proc_bin_data

    # Remove data_obs
    series_data = data.pop('data_obs', None)

    df = pd.DataFrame(data)
    # df.fillna(value=EventYield(0,0), inplace=True)
    # On lxplus, pandas version is 1.2.2, and fillna is bugged (does not accept object values)
    df = df.applymap(lambda x: x if not(isinstance(x, float) and isnan(x)) else EventYield(0,0))

    # Transpose: rows = samples, columns = years
    df = df.transpose()

    # Total yield of MC for each year
    df.loc['Total'] = df.sum()

    if unblind:
        if(series_data is not None):
            logging.debug('data series: %s', series_data)
            df.loc['Data'] = series_data
        else:
            logging.warning('data series is None')

    # Compute total yield for Run2
    df['Run2'] = df.sum(axis=1)

    # Convert to string using the format supplied by command line args
    formatters = [lambda x:x.to_string(fmt=float_format)]*len(df.columns)

    out_string = df.to_latex(formatters=formatters)\
                   .replace('pm', r'$\pm$')\
                   .replace('ZZGTo4LG  ', r'$\PZ\PZ\PGg\to4\Pl\PGg$      ')\
                   .replace('WZGTo3LNuG', r'$\PW\PZ\PGg\to3\Pl\PGn\PGg$  ')\
                   .replace('WZTo3LNu  ', r'$\PW\PZ\to3\Pl\PGn$          ')\
                   .replace('ZZTo4l    ', r'$\Pq\Pq\to\PZ\PZ\to4\Pl$     ')\
                   .replace('ggTo4mu   ', r'$\Pg\Pg\to\PZ\PZ\to4\PGm$    ')\
                   .replace('ggTo2e2mu ', r'$\Pg\Pg\to\PZ\PZ\to2\Pe2\PGm$')\
                   .replace('ggTo4e    ', r'$\Pg\Pg\to\PZ\PZ\to4\Pe$     ')\
                   .replace('ZZZ       ', r'$\PZ\PZ\PZ$                  ')\
                   .replace('WZZ       ', r'$\PW\PZ\PZ$                  ')\
                   .replace('WWZ       ', r'$\PW\PW\PZ$                  ')\
                   .replace('TTZJets   ', r'$\PQt\PAQt\PZ$+jets          ')\
                   .replace('ZGToLLG   ', r'$\PZ\PGg\to\Pl\Pl$           ')\
                   .replace('fake\_leptons', 'Fake leptons               ')\
                   .replace('fake\_photons', 'Fake photons               ')\

    print(out_string)


def get_yield(card, unblind=False, **kwargs):
    '''
    Retrieve the histograms in the files pointed by the shapeMap in the card
    and use their integral and error to construct a dictionary that maps
    {process: {bin: [yield +- error]}}
    '''
    out = {}
    tf_handles = {}
    total = {}
    for bin_name, bin_data in card.exp.items():
        shapeMap_bin = find_in_shapeMap(bin_name, card.shapeMap)
        # logging.debug('shapeMap for bin %s = %s', bin_name, shapeMap_bin)
        mapping = {'CHANNEL': bin_name}

        proc_names = list(bin_data.keys())
        if(unblind):
            proc_names.append('data_obs')

        for proc_name in proc_names:
            shapeMap_proc = find_in_shapeMap(proc_name, shapeMap_bin)
            # logging.debug('shapeMap for process %s = %s', proc_name, shapeMap_proc)
            mapping.update({'PROCESS': proc_name})
            filepath, rootpath, _ = shapeMap_proc
            filepath = Template(filepath).substitute(mapping)
            rootpath = Template(rootpath).substitute(mapping)
            # logging.debug('post sub  filepath=%s  rootpath=%s', filepath, rootpath)

            if not filepath in tf_handles.keys():
                tf = ROOT.TFile(filepath)
                if(not tf or not tf.IsOpen()):
                    raise FileNotFoundError(filepath)
                else:
                    tf_handles[filepath] = tf
                    logging.debug('opened %s', filepath)
            h = tf_handles[filepath].Get(rootpath)
            error = c_double(0.)
            if(h):
                integral = h.IntegralAndError(0,-1, error)
            else:
                logging.warning('Could not get "%s" from %s', rootpath, filepath)
                integral = 0

            out.setdefault(proc_name, {}).setdefault(bin_name, EventYield()) 
            out[proc_name][bin_name] += EventYield(integral, error.value)

            # logging.debug('\t%s %24s = %5.3f +. %5.3f' %(bin_name, proc_name, integral, error.value))

    for filepath, handle in tf_handles.items():
        logging.debug('Closing %s', filepath)
        handle.Close()

    return out


def rename_samples(yield_orig):
    rules_samples = {
        "ggTo([^_]+)_Contin_MCFM701(.*)": r"ggTo\1\2"
    }

    yield_new = {}
    for k_orig in tuple(yield_orig.keys()):
        # logging.debug('k_orig: %s', k_orig)
        k_new = k_orig
        for regex, sub in rules_samples.items():
            # logging.debug('    regex: %s - sub: %s', regex, sub)
            if re.match(regex, k_orig):
                k_new = re.sub(regex, sub, k_orig, count=1)
                # logging.debug('    k_new: %s', k_new)
                break
        yield_new[k_new] = yield_orig.pop(k_orig)

    return yield_new


def main():
    args, args_datacard = parse_args()
    loglevel = args.loglevel.upper() if not args.loglevel.isdigit() else int(args.loglevel)
    logging.basicConfig(format='%% %(levelname)s:%(module)s:%(funcName)s: %(message)s', level=loglevel)

    logging.debug('args          = %s', args)
    # logging.debug('args_datacard = %s', args_datacard)

    cards = []
    for datacard in args.datacards:
        logging.info('reading %s', datacard)
        with open(datacard) as f:
            cards.append( DatacardParser.parseCard(f, options=args_datacard[0]) )

    logging.info('number of cards = %d', len(cards))

    # cards[0].print_structure()
    yields = []
    for card in cards:
        yields.append(get_yield(card, **vars(args)))

    tot_yield = sum_yields(*yields)

    # logging.debug('tot_yield: %s', tot_yield)

    tot_yield = rename_samples(tot_yield)

    print_yield(tot_yield, **vars(args))

    return 0

if __name__ == '__main__':
    main()
