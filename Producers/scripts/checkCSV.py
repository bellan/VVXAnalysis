#!/usr/bin/env python3

################################################################################
# Checks cross sections in the samples csv uasing getXsec (which queries MCM)  # 
#                                                                              #
# Author: A. Mecca (amecca@cern.ch)                                            #
################################################################################

import pandas as pd
import json
from os import path, environ
from math import isnan
from getXsec import getXsec_McM, getXsec_file

def myformatter(number):
    if(number+1 < 1e-5 or number < 1e-5): return "{:.0f}".format(number)
    else:                                 return "{:.3e}".format(number)

for year in [2016]: #(2016, 2017, 2018):
    samples_file = "{:s}/src/VVXAnalysis/Producers/python/samples_{:d}UL_MC.csv".format(environ['CMSSW_BASE'], year)
    print("Year", year, "\tsamples_file:", samples_file)
    with open(samples_file) as f:
        mycsv = pd.read_csv(f, sep=',')
    #mycsv = mycsv[ [ind.replace('#','').replace(' ','') != '' for ind in mycsv['identifier']] ]
    mycsv = mycsv[ [type(dataset) == str for dataset in mycsv['dataset']] ]
    mycsv = mycsv[ [not any(s in dataset for s in ('?', '-missing-')) for dataset in mycsv['dataset']] ]
    
    xsec_storage_file = "xsections{:d}.json".format(year)
    if path.exists(xsec_storage_file):
        with open(xsec_storage_file) as f:
            XsecDict = json.load(f)
    else:
        XsecDict = {}
    
    toUpdate = False
    for ind, series in mycsv.iterrows():
        dataset = series['dataset']
        sample = series['identifier'].strip('#').strip(' ')
        XsecDict.setdefault(dataset, {})
        
        mcm = XsecDict[dataset].get("mcm", None)
        if(mcm is None):
            print(">>>getXsec_McM(%s)" % (dataset))
            mcm = getXsec_McM(dataset)
            XsecDict[dataset]['mcm'] = mcm
            toUpdate = True
            
        internal = XsecDict[dataset].get("internal", None)
        external = XsecDict[dataset].get("external", None)
        if(internal is None or external is None):
            print(">>>getXsec_file(%s)" % (sample))
            Xsec_file = getXsec_file(sample, year=year)
            internal = Xsec_file['internal']
            external = Xsec_file['external']
            XsecDict[dataset].update( {"internal": internal, "external": external} )
            toUpdate = True
        
        csv = series['crossSection=-1']
        csv_dict = XsecDict[dataset].get('csv', None)
        if(csv_dict is None or csv_dict != csv):
            XsecDict[dataset]['csv'] = csv
            toUpdate = True
        
        #if(isnan(csv) or abs(csv-mcm)     /max(csv,mcm,1)      > 0.1): print(sample.ljust(16), " \tmcm: ", mcm     , "\tcsv: ",csv, sep='')
        #if(isnan(csv) or abs(csv-internal)/max(csv,internal,1) > 0.1): print(sample.ljust(16), " \tint: ", internal, "\tcsv: ",csv, sep='')
        #if(isnan(csv) or abs(csv-external)/max(csv,external,1) > 0.1): print(sample.ljust(16), " \text: ", internal, "\tcsv: ",csv, sep='')
        
        if( isnan(csv) or csv == -1 or abs(csv-mcm) / max(csv,mcm,1) > 0.1 or abs(csv-internal)/max(csv,internal,1) > 0.1 or abs(csv-external)/max(csv,external,1) > 0.1 ):
            print("{:23.23s}\tmcm: {:8s}\tcsv: {:8s}".format(sample, myformatter(mcm     ), myformatter(csv)))
            print("{:23.23s}\tint: {:8s}\tcsv: {:8s}".format(sample, myformatter(internal), myformatter(csv)))
            print("{:23.23s}\text: {:8s}\tcsv: {:8s}".format(sample, myformatter(external), myformatter(csv)))
        else:
            print("{:23.23s}\tcsv: {:8s}".format(sample, myformatter(csv)))
        
        print('-'*50)
        
    print()
    
    if(toUpdate):
        with open(xsec_storage_file, "w") as f:
            json.dump(XsecDict, f, indent=2)
