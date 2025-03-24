#!/usr/bin/env python
from argparse import ArgumentParser
import logging
import json

def main(args):
    logging.debug('args: %s', args)

    with open(args.impactsfile) as f:
        impacts = json.load(f)

    # logging.debug('top: %s', impacts.keys())
    logging.debug('POIs: %s', impacts['POIs'])
    logging.debug('method: %s', impacts['method'])
    logging.debug('params: %s', [p['name'] for p in impacts['params']])

    fixes_needed = 0
    for par in impacts['params']:
        fixes_needed += fix_postfit_err(par)

    if(fixes_needed > 0):
        with open(args.impactsfile, 'w') as f:
            json.dump(impacts, f, indent=2)
        logging.info('Overwrote %s', args.impactsfile)

    return 0


def parse_args():
    parser = ArgumentParser()
    parser.add_argument('impactsfile', metavar='FILE')
    parser.add_argument('--log', dest='loglevel', metavar='LEVEL', default='WARNING', help='Level for the python logging module. Can be either a mnemonic string like DEBUG, INFO or WARNING or an integer (lower means more verbose).')

    return parser.parse_args()


def fix_postfit_err(par):
    nfixes = 0
    # par is a refernce to a sub-dict
    pre = par["prefit"]
    fit = par["fit"]
    pre_err_hi = (pre[2] - pre[1])
    pre_err_lo = (pre[1] - pre[0])
    fit_err_hi = (fit[2] - fit[1])
    fit_err_lo = (fit[1] - fit[0])
    if(fit_err_lo > pre_err_lo):
        fix = min(0.95*pre_err_lo, fit_err_hi)
        logging.debug('%s - fix fit_err_lo: %.2f -> %.2f', par['name'], fit_err_lo, fix)
        nfixes += 1
        fit[0] = fit[1] - fix
    if(fit_err_hi > pre_err_hi):
        fix = min(0.95*pre_err_hi, fit_err_lo)
        logging.debug('%s - fix fit_err_hi: %.2f -> %.2f', par['name'], fit_err_hi, fix)
        nfixes += 1
        fit[2] = fit[1] + fix

    return nfixes


if __name__ == '__main__':
    args = parse_args()
    loglevel = args.loglevel.upper() if not args.loglevel.isdigit() else int(args.loglevel)
    logging.basicConfig(format='%(levelname)s:%(module)s:%(funcName)s: %(message)s', level=loglevel)
    
    exit(main(args))
