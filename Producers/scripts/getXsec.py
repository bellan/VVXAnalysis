#!/usr/bin/env python3

################################################################
# Tries to find the XSection of a sample from the dataset name #
#                                                              #
# Author: A. Mecca (amecca@cern.ch)                            #
################################################################


import sys
sys.path.append('/afs/cern.ch/cms/PPD/PdmV/tools/McM/')
from rest import McM
sys.path.append('.')
from request_wrapper import RequestWrapper
import os
import subprocess
import re
import json
from argparse import ArgumentParser

mcm = None
rw  = None

def _getXsec_genparams(gen_params):
    # print('generator_parameters:', json.dumps( [myfilter(p) for p in gen_params], indent=2 ))
    gen_params.sort( 
        key=lambda x: (  # tuples are sorted by their first element, then the second, etc.
            1 - x.get('negative_weights_fraction', 1),  # first the measurement with the least negative weights
            x.get('submission_details', {}).get('submission_date', '')  # then the most recent one
        ), reverse=True )
    # print('generator_parameters:', json.dumps( [myfilter(p) for p in gen_params], indent=2 ))
    return [ p.get('cross_section', None) for p in gen_params ]


def _select_campaign(campaigns, dataset_name, dataset_campaign=None, dataset_datatier=None, **kwargs):
    verbose = kwargs['verbose']
    genCampaigns = [c for c in campaigns if "LHEGEN" in c['prepid']]
    if(len(genCampaigns) == 0):
        print('Error: No LHEGEN campaigns found for dataset "{:s}"'.format(dataset_name))
        return None   # assert len(genCampaigns) > 0, 'Error: No LHEGEN campaigns found for dataset "{:s}"'.format(dataset_name)

    if(verbose): print("LHEGEN campaigns found:", len(genCampaigns))
    if(verbose > 1):
        for c in genCampaigns:
            print('\t'+c['prepid'], '-->', _getXsec_genparams(c['generator_parameters']))
    
    theCampaign = genCampaigns[0]
    # Select the correct genCampaigns if possible
    if(len(genCampaigns) >= 2 and dataset_campaign):
        match = re.search("Run(IX|IV|V?I{0,3}|\d)[^\d]+20(UL)*\d{2}", dataset_campaign)  # matches up to Run9 / RunIX
        if match is not None:
            genCampaigns = [ c for c in genCampaigns if match.group() in c['prepid'] ]
            if(verbose > 1): print("filtered campaigns {:d}:\t".format(len(genCampaigns)), [c['prepid'] for c in genCampaigns], end='\n\n')
    
    if(len(genCampaigns) >= 1):
        theCampaign = sorted(genCampaigns, key=lambda c: c['prepid'], reverse=True )[0]
    return theCampaign


def getXsec_McM(dataset, **kwargs):
    global mcm
    if(mcm is None):
        mcm = McM(dev=False)
    
    verbose = kwargs.get("verbose", 0)
    dataset_split = dataset.strip('/').split('/')
    dataset_name = dataset_split[0]
    
    campaign_requests = mcm.get('requests', query='dataset_name={:s}'.format(dataset_name))
    if( len(campaign_requests) == 0):
        print('Error: No campaigns found for dataset "{:s}"'.format(dataset_name))
        return -2    # assert len(campaign_requests) > 0, 'Error: No campaigns found for dataset "{:s}"'.format(dataset_name)
    
    theCampaign = _select_campaign(campaign_requests, *dataset_split, verbose=verbose)
    if(theCampaign is None):
        return -3    
    if(verbose): print("Chosen campaign: ", theCampaign['prepid'])
    
    ##### campaign found, sort the generator parameters
    gen_params = theCampaign.get('generator_parameters', [])
    if( len(gen_params) == 0): 
        print('Error: The campaign has no generator parameters')
        return -4
    
    return _getXsec_genparams(gen_params)


def getXsec_file(sample, **kwargs):
    verbose = kwargs.get('verbose', 0)
    VVXpath = os.environ['CMSSW_BASE']+'/src/VVXAnalysis'
    if(not '/' in sample):
        sample = VVXpath+'/TreeAnalysis/samples/MC/'+str(kwargs.get("year", 2016))+'/'+sample
    if(not sample.endswith('.root')):
        sample += '.root'
    
    if(verbose): print('INFO: reading from:', sample)
    labels = ['internal', 'internal_std', 'external', 'external_std']
    
    if(not os.path.exists(sample)):
        return dict(zip(labels, [-1, 0, -1, 0]))
    
    command = 'cd %s ; root -l -q \'getInternalXsec.cxx+("%s")\'' % (VVXpath+'/TreeAnalysis/scripts', sample)
    result = subprocess.run(command, shell=True, universal_newlines=True, stdout=subprocess.PIPE, stderr=subprocess.DEVNULL)
    match = re.search("([-\d\.e]+;){3}[-\d\.e]+", result.stdout)
    if(match is None):
        return dict(zip(labels, [-1, 0, -1, 0]))
    
    return dict(zip(labels, [ float(s) for s in match.group().split(';')] ))


def getXsec_xsecdb(dataset, **kwargs):
    global rw
    if(rw is None):
        rw = RequestWrapper()
    
    dataset_name = dataset.strip('/').split('/')[0]
    results = rw.simple_search({'DAS':dataset_name})
    if(kwargs.get('verbose', 0)):
        print(json.dumps(results, indent=2))
    return [float(result.get('cross_section', -1)) for result in results]


if __name__ == '__main__':
    parser = ArgumentParser(description='Script that tries to find the cross section of a sample from the full DAS dataset name')
    parser.add_argument('dataset')
    parser.add_argument('--verbose', '-v', action='count', default=0, help='Increment verbosity')
    parser.add_argument('-q', '--quiet', action='store_const', dest='verbose', const=0, help='Set verbosity to 0')
    parser.add_argument('-t', '--type', choices=['mcm', 'xsecdb', 'internal'], default='mcm', help='Type of cross section to get')
    parser.add_argument('-y', '--year', choices=['2016preVFP', '2016', '2017', '2018'])
    args = parser.parse_args()
    
    if(args.type == 'mcm'):
        xsec = getXsec_McM(args.dataset, verbose=args.verbose)
        print("McM cross-section(s) found:", xsec)
    elif(args.type == 'xsecdb'):
        xsec = getXsec_xsecdb(args.dataset, verbose=args.verbose)
        print("xsecDB cross-section:", xsec)
    elif(args.type == 'internal'):
        d = getXsec_file(args.dataset, verbose=args.verbose)
        print("Internal cross-section from file:", d)
        
