#!/usr/bin/env python
from argparse import ArgumentParser
import os
import logging

import ROOT
from utils23 import lumi_dict
from plotUtils23 import TFileContext, InputDir, InputFile
from plotUtils import makedirs_ok

def main(args):
    # "CREATE" fails if the target file exists
    FOUT_MODE = 'CREATE' if not args.force else 'RECREATE'

    for year in args.years:
        yeardir = os.path.join(args.inputdir, year)
        if(not os.path.exists(yeardir)):
            logging.warning('skipping non-existent year "%s"', yeardir)
            continue

        # If no regions were specified, use all of the existing ones
        if(len(args.regions) == 0):
            regions = [d.split('_')[1] for d in os.listdir(yeardir) if d.startswith(args.analyzer)]
            logging.info('year: %s, regions (%d): %s', year, len(regions), regions)
        else:
            regions = args.regions
        logging.debug('analyzer: %s, year: %-11s, regions: %s', args.analyzer, year, regions)

        # Loop
        for region in regions:
            input_dir  = InputDir(args.inputdir , year=year, region=region, analyzer=args.analyzer)
            output_dir = InputDir(args.outputdir, year=year, region=region, analyzer=args.analyzer)
            if(not os.path.exists(input_dir.path())):
                logging.warning('skipping non-existent region "%s"', input_dir)
                continue

            filelist = os.listdir(input_dir.path())
            logging.debug('analyzer: %s, year: %-11s, region: %-6s, files: %d', input_dir.analyzer, input_dir.year, input_dir.region, len(filelist))

            makedirs_ok(output_dir.path())
            for fname in filelist:
                f_in  = InputFile(input_dir , fname)
                f_out = InputFile(output_dir, fname)
                with TFileContext(f_in.path(), 'READ') as tf_in, TFileContext(f_out.path(), FOUT_MODE) as tf_out:
                    ok = split_from_file(tf_in, tf_out, keyword=args.keyword)
                if(ok != 0):
                    return ok
    return 0


def parse_args():
    parser = ArgumentParser(description='Create a new set of result files (e.g. histograms from an analyzer) from a set of existing results. The histograms which begin with "SYS-<keyword>" are copied and renamed by removing the substring "-<keyword>"')
    parser.add_argument('-y', '--years'    , default=['2016preVFP', '2016postVFP', '2017', '2018'], nargs='+', choices=[k for k in lumi_dict.keys() if k.startswith('2')])
    parser.add_argument('-i', '--inputdir' , default='results', help='Input, top level directory with the results of an analyzer')
    parser.add_argument('-o', '--outputdir', default=None     , help='Output location. Default is just to append "_<keyword>" to the inputdir')
    parser.add_argument('-A', '--analyzer' , default='VVGammaAnalyzer', help='Name of the analyzer, used to compose the path of the input files (default: %(default)s)')
    parser.add_argument('-r', '--regions'  , default=[], nargs='+', metavar='REGION', help='Default is to use all of the existing directories')
    parser.add_argument('-f', '--force'    , action='store_true', help='Overwrite any existing file (default: %(default)s)')
    parser.add_argument('-k', '--keyword'  , default='triboson', help='String to use to identify plots to split (default: %(default)s')
    parser.add_argument('--log', dest='loglevel', metavar='LEVEL', default='WARNING', help='Level for the python logging module. Can be either a mnemonic string like DEBUG, INFO or WARNING or an integer (lower means more verbose).')
    args = parser.parse_args()

    if(args.outputdir is None):
        args.outputdir = args.inputdir.rstrip('/')+'_'+args.keyword

    return args


def split_from_file(tf_in, tf_out, keyword='triboson'):
    stat1 = os.fstat(tf_in .GetFd())
    stat2 = os.fstat(tf_out.GetFd())
    if(stat1.st_ino == stat2.st_ino and stat1.st_dev == stat2.st_dev):
        # Same inode on the same device --> same file!
        logging.error('The two files are the same!')
        return 1

    sysname_to_split = 'SYS-%s'%(keyword)
    keynames = (k.GetName() for k in tf_in.GetListOfKeys())
    sys_triboson_names = (k for k in keynames if k.startswith(sysname_to_split))
    tf_out.cd()
    for k in sys_triboson_names:
        h = tf_in.Get(k)
        newkey = k.replace(sysname_to_split, 'SYS')
        h.Write(newkey)
    return 0


if __name__ == '__main__':
    args = parse_args()
    loglevel = args.loglevel.upper() if not args.loglevel.isdigit() else int(args.loglevel)
    logging.basicConfig(format='%(levelname)s:%(module)s:%(funcName)s: %(message)s', level=loglevel)

    exit(main(args))
