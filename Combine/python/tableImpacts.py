#!/usr/bin/env python3
import json
import logging

def tableImpacts(fname, **kwargs):
    with open(fname, 'r') as f:
        impacts = json.load(f)

    params = sorted(impacts['params'], key=lambda p: p['impact_r'], reverse=True)
    max_name_length = max(len(param['name']) for param in params)
    logging.info('max_name_length: %d', max_name_length)

    digits = kwargs.get('digits', 1)
    fmt_string = '{:%ds} & ' %(max_name_length+1)
    if(kwargs.get('percentage')):
        fmt_string += '/'.join(2*[r'{:+.%df}\%%' %(digits  )])
    else:
        fmt_string += '/'.join(2*['{:+.%df}'   %(digits+2)])
    fmt_string += r' \\'
    logging.debug('fmt_string = %s', fmt_string)

    for param in params:
        var_up = param['r'][0] - param['r'][1]
        var_dn = param['r'][2] - param['r'][1]
        
        if(kwargs.get('percentage')):
            var_up *= 100
            var_dn *= 100

        print(fmt_string.format(param['name'].replace('_', '\_'), var_up, var_dn))


def main():
    from argparse import ArgumentParser

    parser = ArgumentParser('Formats the information on impacts as a LaTeX table. The input is a JSON produced by (a chain of) "combineTool.py -M Impacts -o impacts.json [...]"')
    parser.add_argument('impacts', help='Impacts JSON produced by Combine')
    parser.add_argument('-p', '--percentage', action='store_true', help='Show impacts as percentage')
    parser.add_argument('-n', '--digits', type=int, default=1, help='Number of decimal digits for percentages')
    parser.add_argument('--log', dest='loglevel', metavar='LEVEL', default='WARNING', help='Level for the python logging module. Can be either a mnemonic string like DEBUG, INFO or WARNING or an integer (lower means more verbose).')

    args = parser.parse_args()
    loglevel = args.loglevel.upper() if not args.loglevel.isdigit() else int(args.loglevel)
    logging.basicConfig(format='%(levelname)s:%(module)s:%(funcName)s: %(message)s', level=loglevel)

    logging.info('args: %s', args)

    tableImpacts(args.impacts, **args.__dict__)

if __name__ == '__main__':
    exit(main())
