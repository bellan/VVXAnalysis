#!/usr/bin/env python3

###########################################################################
# Tries to find the XSection of a sample from the dataset name using McM  #
#                                                                         #
# Author: A. Mecca (amecca@cern.ch)                                       #
###########################################################################


import sys
if('/afs/cern.ch/cms/PPD/PdmV/tools/McM/' not in sys.path): sys.path.append('/afs/cern.ch/cms/PPD/PdmV/tools/McM/')
from rest import McM
import re
import json

def getXsec(dataset, **kwargs):
    verbose = kwargs.get("verbose", 0)
    dataset_split = dataset.strip('/').split('/')
    dataset_name = dataset_split[0]
    
    mcm = McM(dev=False)
    campaign_requests = mcm.get('requests', query='dataset_name={:s}'.format(dataset_name))
    if( len(campaign_requests) == 0):
        print('Error: No campaigns found for dataset "{:s}"'.format(dataset_name))
        return -1
    
    genCampaigns = [c for c in campaign_requests if "LHEGEN" in c['prepid']]
    if(len(genCampaigns) == 0):
        print('Error: No LHEGEN campaigns found for dataset "{:s}"'.format(dataset_name))
        return -1

    if(verbose): print("LHEGEN campaigns found:", len(genCampaigns))
    if(verbose > 1):
        for prepid in [c['prepid'] for c in genCampaigns]:
            print("\t"+prepid)
    
    theCampaign = genCampaigns[0]
    # Select the correct genCampaigns if possible
    if(len(genCampaigns) >= 2 and len(dataset_split) > 1):
        match = re.search("Run(IX|IV|V?I{0,3}|\d)[^\d]+20(UL)*\d{2}", dataset_split[1])  # matches up to Run9 / RunIX
        if match is not None:
            genCampaigns = [ c for c in genCampaigns if match.group() in c['prepid'] ]
            if(verbose > 1): print("filtered campaigns {:d}:\t".format(len(genCampaigns)), [c['prepid'] for c in genCampaigns], end='\n\n')
    
    if(len(genCampaigns) >= 1):
        theCampaign = sorted(genCampaigns, key=lambda c: c['prepid'], reverse=True )[0]
    if(verbose): print("Chosen campaign: ", theCampaign['prepid'])
    
    ##### campaign found, choose which generator parameters you want to trust #####
    
    gen_params = theCampaign.get('generator_parameters', [])
    if( len(gen_params) == 0): 
        print('Error: The campaign "{}" has no generator parameters'. format(theCampaign['prepid']))
        return -1
    
    gen_params.sort( 
        key=lambda x: (  # tuples are sorted by their first element, then the second, etc.
            1 - x.get('negative_weights_fraction', 1),  # prefer the measurement with the least negative weights
            x.get('submission_details', {}).get('submission_date', '')  # in doubt, use the most recent one
        ), reverse=True )
    
    return gen_params[0].get('cross_section', -1)


if __name__ == '__main__':
    verbose = 2  # verbose must be set manually, I couldn't be bothered to use OptParse
    assert len(sys.argv) >= 2, "Must specify a dataset"
    if(verbose and len(sys.argv) > 3): print("WARN: Extra argument(s) ignored")
    
    xsec = getXsec(sys.argv[1], verbose=verbose)
    print("xsec:", xsec)
