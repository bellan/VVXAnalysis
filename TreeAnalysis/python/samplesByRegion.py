import sys
import ROOT

##### Define type of samples ##### FIXME: make a class?

qqZZ_pow = [{'files':['ZZTo4l'        ] , 'color':ROOT.kBlue-4  , 'name':'qq #rightarrow ZZ', 'split_prompt_ph':True, 'kfactor': 1.325/1.256}]  # 1.1  #(1.256/1.325)
qqZZ_mad = [{'files':['ZZTo4lamcatnlo'] , 'color':ROOT.kBlue-4  , 'name':'qq #rightarrow ZZ', 'split_prompt_ph':True, 'kfactor': 1.}]

ggZZ     = [{'files': ['ggTo2e2mu_Contin_MCFM701', 'ggTo4e_Contin_MCFM701', 'ggTo4mu_Contin_MCFM701'],
            'color':ROOT.kAzure-4 , 'name':'gg #rightarrow ZZ'   , 'split_prompt_ph':True, 'kfactor': 1.7}]

vbsZZ    = [{'files':['ZZ4lJJ'        ] , 'color':ROOT.kCyan-6  , 'name':'VBS'}]
HZZ      = [{'files':['HZZ'           ] , 'color':ROOT.kCyan-7  , 'name':'higgs'}]

WZ       = [{'files':['WZTo3LNu'      ] , 'color':ROOT.kOrange  , 'name':'WZ'}]
WW       = [{'files':['WWTo2L2Nu'     ] , 'color':ROOT.kYellow-4, 'name':'WW'}]

# tt   with >= 4 leptons
tt_X_4l  = [{'files':['TTZZ', 'TTWW', 'TTZJets']                , 'color':ROOT.kViolet-7, 'name':'t#bar{t}+any'}]
# t(t) with >= 3 leptons
tt_X_3l  = [{'files':tt_X_4l[0]['files']+['TZq', 'TTWJetsToLNu'], 'color':ROOT.kViolet-7, 'name':'t#bar{t}+any'}]
# t(t) with >= 2 leptons
tt_X_2l  = [{'files':tt_X_3l[0]['files']+['tW', 'TTTo2L2Nu']    , 'color':ROOT.kViolet-7, 'name':'t#bar{t}+any'}]
# single top, single lepton
t        = [{'files':['singleT'       ] , 'color':ROOT.kMagenta , 'name':'top'}]
# files missing for now: singleT, tt+gamma (/TTGJets or /TTGamma_Dilept)

ZZTo2L2Nu= [{'files':['ZZTo2L2Nu'     ] , 'color':ROOT.kCyan    , 'name':'ZZ #rightarrow 2l 2#nu'}]
ZZTo2Q2L = [{'files':['ZZTo2Q2L'      ] , 'color':ROOT.kGray    , 'name':'ZZ #rightarrow 2l 2q' }]

W        = [{'files':['WJetsToLNu'    ] , 'color':ROOT.kGreen-1 , 'name':'W+jets' }]
DY       = [{'files':['DYJetsToLL_M50'] , 'color':ROOT.kGreen-9 , 'name':'DY'     }]
ZG       = [{'files':['ZGToLLG'       ] , 'color':ROOT.kGreen+2 , 'name':'Z#gamma'}]
WG       = [{'files':['WGToLNuG'      ] , 'color':ROOT.kGray    , 'name':'W#gamma'}]

triboson = [{'files':['WWW','WWZ','WZZ','ZZZ'], 'color':ROOT.kYellow, 'name':'VVV'}]

WZG      = [{'files':['WZGTo3LNuG'    ] , 'color':ROOT.kMagenta , 'name':'WZ#gamma'}]
ZZG      = [{'files':['ZZGTo4LG'      ] , 'color':ROOT.kRed     , 'name':'ZZ#gamma', 'split_prompt_ph':True}]
ZZGTo2L2jG=[{'files':['ZZGTo2L2jG'    ] , 'color':ROOT.kRed+3   , 'name':'ZZ#gamma #rightarrow 2l 2j'}]
WZGTo2L2jG=[{'files':['WZGTo2L2jG'    ] , 'color':ROOT.kRed-5   , 'name':'WZ#gamma #rightarrow 2l 2j'}]

data_obs =  {'files':['data'          ] , 'color':ROOT.kBlack   , 'name':'Data'}


def is3Lregion(region):
    return region in ('SR3P', 'SR3P_1L', 'SR3P_1F', 'CR000', 'CR001', 'CR010', 'CR011', 'CR100', 'CR101', 'CR110')

def is2Lregion(region):
    return region in ('SR2P', 'SR2P_1L', 'SR2P_1P', 'CR2P_1F')

def isLepCR(region):
    return region in ('CR3P1F', 'CR2P2F', 'CR000', 'CR001', 'CR010', 'CR011', 'CR100', 'CR101', 'CR110')

def getSamplesByRegion(region, MCSet, predType):
    availablePredTypes = ['fromCR', 'lepCR', 'phoCR', 'fullCR', 'fullMC', 'fakeMC']  # Notes: fromCR is a legacy equivalent of lepCR; fullCR = lepCR + phoCR
    if predType not in availablePredTypes:
        raise ValueError("Wrong prediction type ("+predType+"), available: "+str(availablePredTypes))

    if(isLepCR(region) and predType in ('lepCR', 'fullCR')):
        raise ValueError('Prediction "%s" not available in "%s"' %(predType, region))

    if MCSet == 'pow':
        qqZZ = qqZZ_pow
    elif MCSet == 'mad':
        qqZZ = qqZZ_mad
    else: sys.exit("Wrong Set, choose pow or mad")

    tot = WZG + ZZG

    if   region in ['SR4P', 'SR4P_1L', 'SR4P_1F', 'CR3P1F', 'CR2P2F']:
        if   predType == 'fullMC':
            if region in ('SR4P', 'SR4P_1L', 'SR4P_1F'):
                tot += tt_X_3l
            else:
                tot += tt_X_2l
            tot += triboson + qqZZ + ggZZ + WZ + DY + ZG
        elif predType in ('lepCR', 'fromCR'):
            tot += tt_X_4l + triboson + qqZZ + ggZZ
        elif predType == 'phoCR':
            tot += qqZZ + ggZZ # only the prompt part

    elif is3Lregion(region):
        tot += tt_X_3l + triboson + ggZZ + qqZZ
        if   predType == 'fullMC':
            tot += ZZTo2Q2L + ZZTo2L2Nu + WZ + WW + ZG + DY
        elif predType in ('lepCR', 'fromCR'):
            tot += ZG
        elif predType == 'phoCR':
            pass

    elif is2Lregion(region):
        tot += ZZGTo2L2jG + WZGTo2L2jG + tt_X_2l + qqZZ + ggZZ + ZZTo2Q2L + ZZTo2L2Nu + WZ + WW + ZG
        if   predType == 'fullMC':
            tot += DY + WG
        elif predType in ('lepCR', 'fromCR'):
            raise ValueError('No fake lepton method exists for 2L regions')
        elif predType == 'phoCR':
            pass

    elif region == 'CRLFR':
        if   predType == 'fullMC':
            tot += tt_X_2l + qqZZ + ZZTo2Q2L + WZ + ZG + DY
        else:
            raise ValueError('Method "%s" not available for CRLFR'%(predType))
    else:
        raise ValueError('Don\'t know how to categorise region "%s"' %(region))

    tot.reverse()
    return tot
