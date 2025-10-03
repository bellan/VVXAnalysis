#!/usr/bin/env python
from argparse import ArgumentParser
import logging
import json
import re
from math import sqrt

GROUPS = [
    r'^lumi_',
    r'^CMS_pileup_',
    r'^QCDscale_',
    r'^prop_biny',
    r'^pdf_',
    r'fake_photons_norm',
    r'fake_leptons_norm',
    r'^CMS_SMP24014_[^_]+_norm'
]


def main(args):
    logging.debug('args = %s', args)

    with open(args.inputfile) as f:
        data = json.load(f)

    logging.info('Values are expressed in percentage')
    out_template = '%-30s: '+args.format
    remaining = set(p['name'] for p in data['params']) # systs that have not been matched yet
    to_print = [] # Pairs [impact, name] to be printed at the end

    # Match syst names using regex, and compute their combined impact
    for group_re in GROUPS:
        regex = re.compile(group_re)
        gr_params = [p for p in data['params'] if regex.search(p['name'])]
        logging.debug('group "%s" = %d params', group_re, len(gr_params))

        if(len(gr_params) == 0):
            continue
        remaining -= set(p['name'] for p in gr_params)

        impact2 = sum(p['impact_r']**2 for p in gr_params)
        impact = sqrt(impact2)

        to_print.append([impact, '"'+group_re+'"'])

    logging.debug('remaining params = %d', len(remaining))

    # The systs that are not grouped
    for name in remaining:
        param = [p for p in data['params'] if p['name'] == name][0]
        to_print.append([param['impact_r'], name])

    # Sort and print
    to_print.sort(reverse=True)
    for impact, name in to_print:
        print(out_template %(name, 100*impact))

    return 0


def parse_args():
    parser = ArgumentParser(description='Table the post-fit impacts from Combine')
    parser.add_argument('inputfile', metavar='FILE', help='JSON produced by `combineTool.py -M Impacts ...')
    parser.add_argument('-f', '--format', default='%.1f', metavar='FMT', help='printf-style format string (default: %(default)s).')
    parser.add_argument('--log', dest='loglevel', metavar='LEVEL', default='WARNING', help='Level for the python logging module. Can be either a mnemonic string like DEBUG, INFO or WARNING or an integer (lower means more verbose).')

    return parser.parse_args()


if __name__ == '__main__':
    args = parse_args()
    loglevel = args.loglevel.upper() if not args.loglevel.isdigit() else int(args.loglevel)
    logging.basicConfig(format='%(levelname)s:%(module)s:%(funcName)s: %(message)s', level=loglevel)

    exit(main(args))
