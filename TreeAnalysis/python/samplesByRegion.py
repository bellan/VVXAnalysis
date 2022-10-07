import sys
import ROOT

##### Define type of samples ##### FIXME: make a class?

qqZZ_pow = [{"sample":'ZZTo4l'         , "color":ROOT.kBlue-4  , "name":'qq/qg #rightarrow ZZ(+jets)', "kfactor": 1.}]  # 1.1  #(1.256/1.325)
qqZZ_mad = [{"sample":'ZZTo4lamcatnlo' , "color":ROOT.kBlue-4  , "name":'qq/qg #rightarrow ZZ(+jets)', "kfactor": 1.1}]

ggZZ     = [{"sample":'ggTo2e2mu_Contin_MCFM701', "color":ROOT.kAzure-4 , "name":'gg #rightarrow ZZ(+jets)'   , "kfactor": 1.7},
            {"sample":'ggTo4e_Contin_MCFM701'   , "color":ROOT.kAzure-4 , "name":'gg #rightarrow ZZ(+jets)'   , "kfactor": 1.7},
            {"sample":'ggTo4mu_Contin_MCFM701'  , "color":ROOT.kAzure-4 , "name":'gg #rightarrow ZZ(+jets)'   , "kfactor": 1.7}]

vbsZZ    = [{"sample":'ZZ4lJJ'         , "color":ROOT.kCyan-6  , "name":'VBS', "kfactor": 1.0}]
HZZ      = [{"sample":'HZZ'            , "color":ROOT.kCyan-7  , "name":'higgs', "kfactor": 1.0}]
    
WZ       = [{"sample":'WZTo3LNu'       , "color":ROOT.kYellow+2, "name":'WZ', "kfactor": 1.0}]

WW       = [{"sample":'WWTo2L2Nu'      , "color":ROOT.kYellow-7, "name":'WW', "kfactor": 1.0}]

t        = [{'sample': 'singleT'       , "color":ROOT.kGray    , "name":'top'}]

tX       = [{'sample': 'tZq'           , "color":ROOT.kGray    , "name":'tX'},
            {'sample': 'tW'            , "color":ROOT.kGray    , "name":'tX'}]

tt       = [{"sample":'TTTo2L2Nu'      , "color":ROOT.kMagenta+2, "name":'t#bar{t}', "kfactor": 1.0}]

ttX      = [{"sample":'TTWJetsToLNu'   , "color":ROOT.kMagenta+2, "name":'t#bar{t}', "kfactor": 1.0},
            #{"sample":'TTGJets'        , "color":ROOT.kMagenta+2, "name":'t#bar{t}', "kfactor": 1.0}
            {"sample":'TTZJets'        , "color":ROOT.kMagenta+2, "name":'t#bar{t}', "kfactor": 1.0}]

ttXY     = [{"sample": 'TTWW'          , "color":ROOT.kBlue-1  , "name":'ttXY', "kfactor": 1.0},
            {"sample": 'TTZZ'          , "color":ROOT.kBlue-1  , "name":'ttXY', "kfactor": 1.0}]

ZZTo2L2Nu= [{"sample": 'ZZTo2L2Nu'     , "color":ROOT.kGray    , "name":'ZZTo2L2Nu'}]

ZZTo2Q2L = [{"sample": 'ZZTo2Q2L'      , "color":ROOT.kGray    , "name":'ZZTo2Q2L' }]

W        = [{"sample":'WJetsToLNu'     , "color":ROOT.kGreen-1 , "name":'W+jets', "kfactor": 1.0}]

DY       = [{"sample":'DYJetsToLL_M50' , "color":ROOT.kGreen-10, "name":'DY', "kfactor": 1.0}]

ZG       = [{"sample":'ZGToLLG'        , "color":ROOT.kGreen-4 , "name":'Z\gamma', "kfactor": 1.0}]

WG       = [{"sample":'WGToLNuG'       , "color":ROOT.kGray    , "name":'W\gamma'}]

triboson = [{"sample":'WWW'            , "color":ROOT.kOrange+1, "name":'VVV', "kfactor": 1.0},
            {"sample":'WWZ'            , "color":ROOT.kOrange+1, "name":'VVV', "kfactor": 1.0},
            {"sample":'WZZ'            , "color":ROOT.kOrange+1, "name":'VVV', "kfactor": 1.0},
            {"sample":'ZZZ'            , "color":ROOT.kOrange+1, "name":'VVV', "kfactor": 1.0}]

#ttZ      = [{"sample":'TTZJets_M10_MLM', "color":ROOT.kOrange-5, "name":'t#bar{t}Z', "kfactor": 1.0}]

WZG      = [{"sample":'WZGTo3LNuG'     , "color":ROOT.kMagenta , "name":'WZ\gamma', "kfactor": 1.0}]
ZZG      = [{"sample":'ZZGTo4LG'       , "color":ROOT.kRed     , "name":'ZZ\gamma', "kfactor": 1.0}]

data_obs = [{"sample":'data'           , "color":ROOT.kBlack   , "name":'Data', "kfactor": 1.0}]


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

    if   predType == 'fullMC':
        tot += WZ + DY + ZG + WG + WW + W + tt + ttX + ZZTo2L2Nu + ZZTo2Q2L

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
        elif region in ['SR2P', 'SR2P_1L', 'SR2P_1P', 'CR2P_1F']:
            tot += WZ + DY + ZG + WG + WW + W + tt + ttX + tX + ZZTo2L2Nu + ZZTo2Q2L
        else:
            tot += WZ + DY + ZG + WG + WW + W + tt + ttX + tX + ZZTo2L2Nu + ZZTo2Q2L
    
    return tot
