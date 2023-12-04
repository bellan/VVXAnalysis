#!/usr/bin/env python3

################################################################################
# Checks cross sections in the samples csv uasing getXsec (which queries MCM)  # 
#                                                                              #
# Author: A. Mecca (amecca@cern.ch)                                            #
################################################################################

import pandas as pd
import json
import os
from math import isnan, floor, log10
from argparse import ArgumentParser
import logging
from getXsec import getXsec_McM, getXsec_xsecdb, getXsec_file

def myformatter(number):
    if(number+1 < 1e-5 or number < 1e-5): return "{:.0f}".format(number)
    else:                                 return "{:.3e}".format(number)

def isInvalid(xsec):
        if(xsec is None):
            return True
        elif(hasattr(xsec, '__iter__')):
            return len(xsec) == 0
        elif(xsec <= 0):
            return True
        return False

def check_csv(mycsv, year=None, types=[]):
    # Cache results to file to avoid time-expensive queries to McM
    xsec_storage_file = "xsections{}.json".format(year)
    if os.path.exists(xsec_storage_file):
        with open(xsec_storage_file) as f:
            XsecDict = json.load(f)
    else:
        XsecDict = {}
    
    jsonUpdate = False  # Do we have to update cached results?
    csvUpdate  = False  # Do we have to fix the csv?
    
    # Loop on the rows of the CSV
    for ind, series in mycsv.iterrows():
        dataset = series['dataset']
        sample = series['identifier'].strip('#').strip()
        
        # Sanitize input
        if( not type(dataset) == str or any(s in dataset for s in ('?', '-missing-')) ):
            continue
        dataset = dataset.strip()
        
        ### Get or update cross-sections ###
        XsecDict.setdefault(dataset, {})
        XsecDict[dataset]['identifier'] = sample
        
        # McM
        mcm = XsecDict[dataset].get("mcm", None)
        if(isInvalid(mcm) or 'mcm' in args.types):
            print(">>>getXsec_McM(%s)" % (dataset), end = ' --> ')
            mcm = getXsec_McM(dataset)
            XsecDict[dataset]['mcm'] = mcm
            print(mcm)
            jsonUpdate = True
        
        # xsecDB
        xsecdb = XsecDict[dataset].get("xsecdb", None)
        if(isInvalid(xsecdb) or 'xsecdb' in args.types):
            print(">>>getXsec_xsecDB(%s)" % (dataset), end=' --> ')
            xsecdb = getXsec_xsecdb(dataset)
            XsecDict[dataset]['xsecdb'] = xsecdb
            print(xsecdb)
            jsonUpdate = True
        
        # the one written to file (internal + "external" --> written in the CSV when the ntuple was produced)
        internal = XsecDict[dataset].get("internal", None)
        external = XsecDict[dataset].get("external", None)
        if(internal is None or external is None or 'internal' in args.types):
            print(">>>getXsec_file(%s)" % (sample), end=' --> ')
            Xsec_file = getXsec_file(sample, year=year)
            internal = Xsec_file['internal']
            external = Xsec_file['external']
            XsecDict[dataset].update( {"internal": internal, "external": external} )
            print('int:', internal, ' ext:', external)
            jsonUpdate = True
        
        # external --> written NOW in the CSV
        csv = series['crossSection=-1']
        csv_dict = XsecDict[dataset].get('csv', None)
        if(csv_dict is None or csv_dict != csv):
            XsecDict[dataset]['csv'] = csv
            jsonUpdate = True
        
        def is_diff(a, b): return a > 0 and b > 0 and abs(a - b) / max(a, b, 1e-10) > 0.1
        
        # print( '{:23.23s}\t{:8s}\t{:8s}\t{:8s}\t{:8s}\t{:8s}'.format(sample, myformatter(csv), myformatter(mcm), myformatter(xsecdb), myformatter(internal), myformatter(external)), end='\t' )
        # if  (mcm      > 0):
        #     if(is_diff(csv, mcm     )):
        #         print('update with McM')
        #         series['crossSection=-1'] = mcm
        #         csvUpdate = True
        #     else:
        #         print('OK (McM)')
        # elif(internal > 0):
        #     if(is_diff(csv, internal)):
        #         print('update with internal')
        #         series['crossSection=-1'] = internal
        #         csvUpdate = True
        #     else:
        #         print('OK (int)')
    
    print()  # end of loop on one year
    
    # if(csvUpdate):
    #     fname = "{:s}/src/VVXAnalysis/Producers/python/samples_{:d}UL_MC_new.csv".format(environ['CMSSW_BASE'], year)
    #     mycsv.to_csv(fname, index=False)
    #     print('INFO: Written to "{:s}"'.format(fname))
            
    
    if(jsonUpdate):
        with open(xsec_storage_file, "w") as f:
            json.dump(XsecDict, f, indent=2)

def main():
    parser = parser = ArgumentParser()
    parser.add_argument('-y', '--year' , choices=['2016preVFP', '2016postVFP', '2017', '2018'], default='2018')
    parser.add_argument('-t', '--types', choices=['mcm', 'xsecdb', 'internal'], nargs='+', default=[])
    parser.add_argument('--log', dest='loglevel', metavar='LEVEL', default='WARNING', help='Level for the python logging module. Can be either a mnemonic string like DEBUG, INFO or WARNING or an integer (lower means more verbose).')

    args = parser.parse_args()
    loglevel = args.loglevel.upper() if not args.loglevel.isdigit() else int(args.loglevel)
    logging.basicConfig(format='%(levelname)s:%(module)s:%(funcName)s: %(message)s', level=loglevel)

    # Find VVXAnalysis
    if('CMSSW_BASE' in os.environ):
        VVX_dir = os.path.join(os.environ['CMSSW_BASE'], 'src', 'VVXAnalysis')
    else:
        split_path = os.path.realpath('.').split(os.sep)
        VVX_position = split_path.index('VVXAnalysis')
        VVX_dir = os.path.join(split_path[:VVX_position+1])
    csv_dir = os.path.join(VVX_dir, 'Producers', 'python')
    csv_file = os.path.join(csv_dir, 'samples_{}UL_MC.csv'.format(year))
    logging.info("Year: %s - csv_file:", year, csv_file)

    with open(csv_file) as f:
        mycsv = pd.read_csv(f, sep=',')

    check_csv(mycsv, vars(args))

if __name__ == '__main__':
    exit(main())
