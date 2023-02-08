#!/usr/bin/env python3

################################################################
# Tries to find the XSection of a sample from the dataset name #
#                                                              #
# Author: A. Mecca (amecca@cern.ch)                            #
################################################################


import sys
sys.path.append('/afs/cern.ch/cms/PPD/PdmV/tools/McM/')
from rest import McM
import re
import json
from argparse import ArgumentParser

mcm = McM(dev=False)

def _getXsec_genparams(gen_params):
    # print('generator_parameters:', json.dumps( [myfilter(p) for p in gen_params], indent=2 ))
    gen_params.sort( 
        key=lambda x: (  # tuples are sorted by their first element, then the second, etc.
            1 - x.get('negative_weights_fraction', 1),  # first the measurement with the least negative weights
            x.get('submission_details', {}).get('submission_date', '')  # then the most recent one
        ), reverse=True )
    # print('generator_parameters:', json.dumps( [myfilter(p) for p in gen_params], indent=2 ))
    return [ p.get('cross_section', None) for p in gen_params ]

def getXsec(dataset, **kwargs):
    verbose = kwargs.get("verbose", 0)
    dataset_split = dataset.strip('/').split('/')
    dataset_name = dataset_split[0]
    
    campaign_requests = mcm.get('requests', query='dataset_name={:s}'.format(dataset_name))
    if( len(campaign_requests) == 0):
        print('Error: No campaigns found for dataset "{:s}"'.format(dataset_name))
        return -2    # assert len(campaign_requests) > 0, 'Error: No campaigns found for dataset "{:s}"'.format(dataset_name)
    
    genCampaigns = [c for c in campaign_requests if "LHEGEN" in c['prepid']]
    if(len(genCampaigns) == 0):
        print('Error: No LHEGEN campaigns found for dataset "{:s}"'.format(dataset_name))
        return -3    # assert len(genCampaigns) > 0, 'Error: No LHEGEN campaigns found for dataset "{:s}"'.format(dataset_name)

    if(verbose): print("LHEGEN campaigns found:", len(genCampaigns))
    if(verbose > 1):
        for c in genCampaigns:
            print('\t'+c['prepid'], '-->', _getXsec_genparams(c['generator_parameters']))
    
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
    
    ##### campaign found, sort the generator parameters
    gen_params = theCampaign.get('generator_parameters', [])
    if( len(gen_params) == 0): 
        print('Error: The campaign has no generator parameters')
        return -4
    
    return _getXsec_genparams(gen_params)


if __name__ == '__main__':
    parser = ArgumentParser(description='Script that tries to find the cross section of a sample from the full DAS dataset name')
    parser.add_argument('dataset')
    parser.add_argument('--verbose', '-v', action='count', default=0)
    args = parser.parse_args()

    xsec = getXsec(args.dataset, verbose=args.verbose)
    print("cross-section(s) found:", xsec)
