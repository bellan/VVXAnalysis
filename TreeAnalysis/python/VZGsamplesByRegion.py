import sys
import ROOT

##### Define type of samples ##### FIXME: make a class?

TTTo2L2Nu= [{'files':['TTTo2L2Nu'     ] , 'color':ROOT.kViolet-7, 'name':'t#bar{t}+any'}]

DY       = [{'files':['DYJetsToLL_M50'] , 'color':ROOT.kGreen-9 , 'name':'DY'     , 'skip_prompt_ph':True}]
ZG       = [{'files':['ZGToLLG'       ] , 'color':ROOT.kGreen+2 , 'name':'Z#gamma', 'skip_nonprompt_ph':True}]

VZG      =[{'files':['VZG'    ] , 'color':ROOT.kRed   , 'name':'VZ#gamma'}]
FSR      =[{'files':['FSR'    ] , 'color':ROOT.kRed-3   , 'name':'VZ+FSR'}]
data_obs =[{'files':['data_obs'          ] , 'color':ROOT.kBlack   , 'name':'Data'}]


def is3Lregion(region):
    return region in ('SR3P', 'SR3P_1L', 'SR3P_1F', 'CR000', 'CR001', 'CR010', 'CR011', 'CR100', 'CR101', 'CR110')

def is2Lregion(region):
    return region in ('SR2P', 'SR2PFJ', 'SR2P_1L', 'SR2P_1P', 'CR2P_1F', 'CRDY')

def isLepCR(region):
    return region in ('CR3P1F', 'CR2P2F', 'CR000', 'CR001', 'CR010', 'CR011', 'CR100', 'CR101', 'CR110')

def getVZGSamplesByRegion(region, MCSet, predType):
    if is2Lregion(region):
        tot = VZG + TTTo2L2Nu + DY + ZG
    else:
        raise ValueError('Don\'t know how to categorise region "%s"' %(region))

    tot.reverse()
    return tot
