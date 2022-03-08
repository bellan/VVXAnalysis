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

for year in (2016, 2017, 2018):
    samples_file = "{:s}/src/VVXAnalysis/Producers/python/samples_{:d}UL_MC.csv".format(environ['CMSSW_BASE'], year)
    print("\nYear", year, "\tsamples_file:", samples_file)
    with open(samples_file) as f:
        mycsv = pd.read_csv(f, sep=',')
    #mycsv = mycsv[ [ind.replace('#','').replace(' ','') != '' for ind in mycsv['identifier']] ]
    mycsv = mycsv[ [type(dataset) == str for dataset in mycsv['dataset']] ]
    mycsv = mycsv[ [not any(s in dataset for s in ('?', '-missing-')) for dataset in mycsv['dataset']] ]
    
    xsec_storage_file = "xsections{:d}_fromMCM.json".format(year)
    if path.exists(xsec_storage_file):
        with open(xsec_storage_file) as f:
            mcmXsec = json.load(f)
    else:
        mcmXsec = {}
        from getXsec import getXsec
        for ind, series in mycsv.iterrows():
            dataset = str(series['dataset'])
            # if(any( a in dataset for a in ['?', '-missing-'] )): continue
            mcmXsec[dataset] = getXsec(dataset)
        with open(xsec_storage_file, "w") as f:
            json.dump(mcmXsec, f, indent=2)
    
    toUpdate = False
    for ind, series in mycsv.iterrows():
        dataset = series['dataset']
        # if(any( [ a in dataset for a in ['?', '-missing-'] ] )): continue
        mcm = mcmXsec.get(dataset, None)
        if(mcm is None):
            mcm = getXsec(dataset)
            mcmXsec[dataset] = mcm
            toUpdate = True
        csv = series['crossSection=-1']
        if(isnan(csv) or abs(mcm-csv)/mcm > 0.1): print("mcm: ", mcm, " \tcsv: ", csv, " \t", dataset, sep='')
    
    if(toUpdate):
        with open(xsec_storage_file, "w") as f:
            json.dump(mcmXsec, f, indent=2)
