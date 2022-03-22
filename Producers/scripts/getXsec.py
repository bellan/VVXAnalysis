#!/usr/bin/env python3

###########################################################################
# Tries to find the XSection of a sample from the dataset name using McM  #
#                                                                         #
# Author: A. Mecca (amecca@cern.ch)                                       #
###########################################################################


import sys
import os
import subprocess
import re
import json
if('/afs/cern.ch/cms/PPD/PdmV/tools/McM/' not in sys.path): sys.path.append('/afs/cern.ch/cms/PPD/PdmV/tools/McM/')
from rest import McM

def findCampaign(dataset, **kwargs):
    verbose = kwargs.get("verbose", 0)
    dataset_split = dataset.strip('/').split('/')
    dataset_name = dataset_split[0]
    
    mcm = McM(dev=False)
    campaign_requests = mcm.get('requests', query='dataset_name={:s}'.format(dataset_name))
    if( len(campaign_requests) == 0):
        print('Error: No campaigns found for dataset "{:s}"'.format(dataset_name))
        return None
        
    genCampaigns = [c for c in campaign_requests if "LHEGEN" in c['prepid']]
    if(len(genCampaigns) == 0):
        print('Error: No LHEGEN campaigns found for dataset "{:s}"'.format(dataset_name))
        return None

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
        
    return theCampaign

def recomputeXsec(campaign):
    #TEMP
    # print(campaign.keys())
    proddir = "fragments"
    if(not os.path.exists(proddir)):
        os.mkdir(proddir)
    os.chdir(proddir)
    gridpack = re.search("/cvmfs/[^,' ]+\.(tgz|tar\.gz)", campaign['fragment'])
    if(gridpack is None):
        print("Error parsing the fragment")
        return -1
    
    # subprocess.run(["" gridpack.group()]

    # fragname = campaign['prepid']+"_cfg.py"
    # with open(fragname, 'w') as f:
    #     f.write(campaign['fragment'].replace('5000', '100'))
    # subprocess.run('cmsRun -t -n 8 ./{} &> {}.log'.format(fragname, campaign['prepid']), shell=True)
    #/TEMP



def getXsec_McM(dataset, **kwargs):
    theCampaign = findCampaign(dataset, **kwargs)
    if(theCampaign is None):
        return -1
        
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


def getXsec_file(sample, **kwargs):
    VVXpath = os.environ['CMSSW_BASE']+'/src/VVXAnalysis'
    abspath = kwargs.get("abspath", False)
    if(not abspath):
        sample = VVXpath+'/TreeAnalysis/samples/MC/'+str(kwargs.get("year", 2016))+'/'+sample
    if(not sample.endswith('.root')):
        sample += '.root'
    
    labels = ['internal', 'internal_std', 'external', 'external_std']
    
    if(not os.path.exists(sample)):
        return dict(zip(labels, [-1, 0, -1, 0]))
    
    command = 'cd %s ; root -l -q \'getInternalXsec.cxx+("%s")\'' % (VVXpath+'/TreeAnalysis/scripts', sample)
    result = subprocess.run(command, shell=True, universal_newlines=True, stdout=subprocess.PIPE, stderr=subprocess.DEVNULL)
    match = re.search("([-\d\.e]+;){3}[-\d\.e]+", result.stdout)
    if(match is None):
        return dict(zip(labels, [-1, 0, -1, 0]))
    
    return dict(zip(labels, [ float(s) for s in match.group().split(';')] ))
    


if __name__ == '__main__':
    verbose = 0  # verbose must be set manually, I couldn't be bothered to use OptParse
    assert len(sys.argv) >= 2, "Must specify a dataset"
    if(verbose and len(sys.argv) > 3): print("WARN: Extra argument(s) ignored")
    
    campaign = findCampaign(sys.argv[1], verbose=verbose)
    if(campaign is None):
        exit(1)
    recomputeXsec(campaign)
        
    # xsec = getXsec(sys.argv[1], verbose=verbose)
    # print("xsec:", xsec)
