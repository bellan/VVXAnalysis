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

for year in [2016]:  # (2016, 2017, 2018):
    samples_file = "{:s}/src/VVXAnalysis/Producers/python/samples_{:d}UL_MC.csv".format(environ['CMSSW_BASE'], year)
    print("Year", year, "\tsamples_file:", samples_file)
    with open(samples_file) as f:
        mycsv = pd.read_csv(f, sep=',')
    
    # Sanitize input
    # mycsv = mycsv[ [ind.replace('#','').replace(' ','') != '' for ind in mycsv['identifier']] ]
    # mycsv = mycsv[ [type(dataset) == str for dataset in mycsv['dataset']] ]
    # mycsv = mycsv[ [not any(s in dataset for s in ('?', '-missing-')) for dataset in mycsv['dataset']] ]
    
    # Cache results to file to avoid time-expensive queries to McM
    xsec_storage_file = "xsections{:d}.json".format(year)
    if path.exists(xsec_storage_file):
        with open(xsec_storage_file) as f:
            XsecDict = json.load(f)
    else:
        XsecDict = {}
    
    jsonUpdate = False  # Do we have to update cached results?
    csvUpdate  = False  # Do we have to fix the csv?
    # print("{:23.23s}\t{:8s}\t{:8s}\t{:8s}\t{:8s}\t{:s}".format("SAMPLE", "CSV", "MCM", "INTERNAL", "EXTERNAL", "STATUS"))
    print("{:23.23s}\t{:8s}\t{:8s}\t{:8s}\t{:8s}".format("SAMPLE", "CSV", "MCM", "INTERNAL", "EXTERNAL"))
    
    # Loop on the rows of the CSV
    for ind, series in mycsv.iterrows():
        dataset = series['dataset']
        sample = series['identifier'].strip('#').strip(' ')
        
        # Sanitize input
        if( not type(dataset) == str or any(s in dataset for s in ('?', '-missing-')) ):
            continue

        XsecDict.setdefault(dataset, {})
        
        mcm = XsecDict[dataset].get("mcm", None)
        if(mcm is None):
            print(">>>getXsec_McM(%s)" % (dataset))
            mcm = getXsec_McM(dataset)
            XsecDict[dataset]['mcm'] = mcm
            jsonUpdate = True
            
        internal = XsecDict[dataset].get("internal", None)
        external = XsecDict[dataset].get("external", None)
        if(internal is None or external is None):
            print(">>>getXsec_file(%s)" % (sample))
            Xsec_file = getXsec_file(sample, year=year)
            internal = Xsec_file['internal']
            external = Xsec_file['external']
            XsecDict[dataset].update( {"internal": internal, "external": external} )
            jsonUpdate = True
        
        csv = series['crossSection=-1']
        csv_dict = XsecDict[dataset].get('csv', None)
        if(csv_dict is None or csv_dict != csv):
            XsecDict[dataset]['csv'] = csv
            jsonUpdate = True
        
        # Presentation of data
        # if( isnan(csv) or csv == -1 or abs(csv-mcm) / max(csv,mcm,1e-10) > 0.1 or abs(csv-internal)/max(csv,internal,1e-10) > 0.1 or abs(csv-external)/max(csv,external,1e-10) > 0.1 ):
        #     sample_status = 'PROBLEM'
        # else:
        #     sample_status = 'OK'
        # print('{:23.23s}\t{:8s}\t{:8s}\t{:8s}\t{:8s}\t{:s}'.format(sample, myformatter(csv), myformatter(mcm), myformatter(internal), myformatter(external), sample_status))
        def is_diff(a, b): return a > 0 and b > 0 and abs(a - b) / max(a, b, 1e-10) > 0.1
        
        print( '{:23.23s}\t{:8s}\t{:8s}\t{:8s}\t{:8s}'.format(sample, myformatter(csv), myformatter(mcm), myformatter(internal), myformatter(external)), end='\t' )
        if  (mcm      > 0):
            if(is_diff(csv, mcm     )):
                print('update with McM')
                series['crossSection=-1'] = mcm
                csvUpdate = True
            else:
                print('OK (McM)')
        elif(internal > 0):
            if(is_diff(csv, internal)):
                print('update with internal')
                series['crossSection=-1'] = internal
                csvUpdate = True
            else:
                print('OK (int)')
    
    print()  # end of loop on one year
    
    if(csvUpdate):
        fname = "{:s}/src/VVXAnalysis/Producers/python/samples_{:d}UL_MC_new.csv".format(environ['CMSSW_BASE'], year)
        mycsv.to_csv(fname, index=False)
        print('> Written to "{:s}"'.format(fname)) 
            
    
    if(jsonUpdate):
        with open(xsec_storage_file, "w") as f:
            json.dump(XsecDict, f, indent=2)
