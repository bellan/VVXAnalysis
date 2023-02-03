import sys
import ROOT

##### Define type of samples ##### FIXME: make a class?

qqZZ_pow = [{'files':['ZZTo4l'        ] , 'color':ROOT.kBlue-4  , 'name':'qq/qg #rightarrow ZZ(+jets)', 'kfactor': 1.}]  # 1.1  #(1.256/1.325)
qqZZ_mad = [{'files':['ZZTo4lamcatnlo'] , 'color':ROOT.kBlue-4  , 'name':'qq/qg #rightarrow ZZ(+jets)', 'kfactor': 1.}]

ggZZ     = [{'files': ['ggTo2e2mu_Contin_MCFM701', 'ggTo4e_Contin_MCFM701', 'ggTo4mu_Contin_MCFM701'],
            'color':ROOT.kAzure-4 , 'name':'gg #rightarrow ZZ(+jets)'   , 'kfactor': 1.7}]

vbsZZ    = [{'files':['ZZ4lJJ'        ] , 'color':ROOT.kCyan-6  , 'name':'VBS'}]
HZZ      = [{'files':['HZZ'           ] , 'color':ROOT.kCyan-7  , 'name':'higgs'}]

WZ       = [{'files':['WZTo3LNu'      ] , 'color':ROOT.kYellow+2, 'name':'WZ'}]
WW       = [{'files':['WWTo2L2Nu'     ] , 'color':ROOT.kYellow-4, 'name':'WW'}]

t        = [{'files':['singleT'       ] , 'color':ROOT.kMagenta , 'name':'top'}]
tX       = [{'files':['tZq','tW'      ] , 'color':ROOT.kMagenta-9, 'name':'tX'}]
tt       = [{'files':['TTTo2L2Nu'     ] , 'color':ROOT.kMagenta+2, 'name':'t#bar{t}'}]
ttX      = [{'files':['TTWJetsToLNu', 'TTZJets', 'TTGJets'] , 'color':ROOT.kViolet-7, 'name':'ttX'}]
ttXY     = [{'files':['TTWW','TTZZ'   ] , 'color':ROOT.kBlue-1  , 'name':'ttXY'}]

ZZTo2L2Nu= [{'files':['ZZTo2L2Nu'     ] , 'color':ROOT.kGray    , 'name':'ZZTo2L2Nu'}]
ZZTo2Q2L = [{'files':['ZZTo2Q2L'      ] , 'color':ROOT.kGray    , 'name':'ZZTo2Q2L' }]

W        = [{'files':['WJetsToLNu'    ] , 'color':ROOT.kGreen-1 , 'name':'W+jets' }]
DY       = [{'files':['DYJetsToLL_M50'] , 'color':ROOT.kGreen-10, 'name':'DY'     }]
ZG       = [{'files':['ZGToLLG'       ] , 'color':ROOT.kGreen-4 , 'name':'Z\gamma'}]
WG       = [{'files':['WGToLNuG'      ] , 'color':ROOT.kGray    , 'name':'W\gamma'}]

triboson = [{'files':['WWW','WWZ','WZZ','ZZZ'], 'color':ROOT.kOrange+1, 'name':'VVV'}]

ttZ      = [{'files':['TTZJets_M10_MLM'], 'color':ROOT.kOrange-5, 'name':'t#bar{t}Z'}]

WZG      = [{'files':['WZGTo3LNuG'    ] , 'color':ROOT.kMagenta , 'name':'WZ\gamma'}]
ZZG      = [{'files':['ZZGTo4LG'      ] , 'color':ROOT.kRed     , 'name':'ZZ\gamma'}]
ZZGTo2L2jG=[{'files':['ZZGTo2L2jG'    ] , 'color':ROOT.kRed+3   , 'name':'ZZ\gamma #rightarrow 2\ell 2j'}]
WZGTo2L2jG=[{'files':['WZGTo2L2jG'    ] , 'color':ROOT.kRed-5   , 'name':'WZ\gamma #rightarrow 2\ell 2j'}]

data_obs =  {'files':['data'          ] , 'color':ROOT.kBlack   , 'name':'Data'}


def is2Lregion(region):
    return region in ['SR2P', 'SR2P_1L', 'SR2P_1P', 'CR2P_1F']


def getSamplesByRegion(region, MCSet, predType):
    if predType not in ['fromCR', 'fullMC', 'fakeMC']:
        sys.exit("Wrong prediction type, fromCR from MC still needs to be added")
    
    qqZZ = {}
    if MCSet == 'pow':
        qqZZ = qqZZ_pow
    elif MCSet == 'mad':
        qqZZ = qqZZ_mad
    else: sys.exit("Wrong Set, choose pow or mad")

    tot = WZG + ZZG + qqZZ + ggZZ + triboson + ttXY #+ vbsZZ + HZZ
    if is2Lregion(region):
        tot += ZZGTo2L2jG + WZGTo2L2jG

    if   predType == 'fullMC':
        tot += WZ + DY + ZG + WG + WW + W + tt + ttX # + ZZTo2L2Nu + ZZTo2Q2L

    elif predType in ['fromCR', 'fakeMC']:
        if   region in ['SR4P', 'SR4P_1L', 'SR4P_1F']:
            pass
        elif region in ['CR3P1F', 'CR2P2F']:
            pass
        
        elif region in ['SR3P', 'SR3P_1L', 'SR3P_1F']:
            tot += WZ + ttX
        elif region in ['CR000', 'CR001', 'CR010', 'CR011', 'CR100', 'CR101', 'CRLFR']:
            tot += WZ + ttX + tX + DY + ZG
        elif region == 'CR110':
            tot += WZ + ttX + tX + DY + ZG
        elif is2Lregion(region):
            tot += WZ + DY + ZG + WG + WW + W + tt + ttX + tX + ZZTo2L2Nu + ZZTo2Q2L
        else:
            tot += WZ + DY + ZG + WG + WW + W + tt + ttX + tX + ZZTo2L2Nu + ZZTo2Q2L
    
    return tot
