#! /usr/bin/env python2
import sys
import os
import math
import operator
import re
from copy import deepcopy
from optparse import OptionParser
import ROOT
ROOT.PyConfig.IgnoreCommandLineOptions = True
import CrossInfo
from CrossInfo import* 
from ROOT import TH1F,TCanvas, TLegend
import plotUtils  # GetPredictionsPlot, GetDataPlot
import CMS_lumi, tdrstyle
import PersonalInfo

Lumi   = 35900
regions = ['SR4P', 'CR3P1F' , 'CR2P2F' , 'SR4P_1L', 'SR4P_1P', 'CR4P_1F', 'CR4L',    
           'SR3P', 'CR110'  , 'CR101'  , 'CR011'  , 'CR100'  , 'CR001'  , 'CR010', 'CR000', 'SR3P_1L', 'SR3P_1P', 'CR3P_1F', 'CRLFR', 'CR3L',
           'SR2P', 'SR2P_1L', 'SR2P_1P', 'CR2P_1F', 
           'SR_HZZ', 'CR2P2F_HZZ', 'CR3P1F_HZZ', 'CR_HZZ']

parser = OptionParser()

parser.add_option("-r", "--region", dest="region", type="choice", choices=regions,
                  default="SR4P",
                  help="Region type are {0:s}. Default is SR.".format(', '.join(regions)))

parser.add_option("-f", "--finalstate", dest="FinalState",
                  default="4l",
                  help="Final state are 4l, 4m, 2e2m and 4e. Default is 4l")

parser.add_option("--nodata", dest="noData",
                  action="store_true",
                  default=False,
                  help="Data is True or False to plot or not the data. Default is True")

parser.add_option("-t", "--type", dest="Type",
                  default="Mass",
                  help= "type type to choose the  plot you want. Mass, Jets, DeltaEta, mjj")

parser.add_option("-S", "--Save", dest="SavePlot",
                  action="store_true",
                  default=False,
                  help="Save plot option, default is False")

parser.add_option("-m", "--mcset", dest="mcSet", type="choice", choices=['mad', 'pow'],
                  default="mad",
                  help= "Monte Carlo Set, pow for Powheg, mad for amcatnlo")

parser.add_option("-p", "--prediction-type", dest="predType",
                  default="fromCR",
                  help= "Type of prediction. fromCR = non-prompt leptons from CRs, rare background from MC; fullMC = all from MC. Default is fromCR")

parser.add_option("-l", "--lumiProj", dest="LumiProj",
                  default="",
                  help="Lumi projection")

parser.add_option("-o", "--outputDir", dest="outputDir",
                  default="last",
                  help="Directory where save plots")

parser.add_option("-A", "--Analysis", dest="Analysis", type="choice", choices=['VVXAnalyzer', 'VVGammaAnalyzer', 'ZZAnalyzer'],
                  default="VVXAnalyzer",
                  help="Analysis. Default is VVX; other options are ZZ, VBS and VVGamma")

parser.add_option("-y", "--year", dest="year",
                  default="2016",
                  help= "valid inputs are 2016, 2017, 2018")


#REMEMBER ADD DEFINTION PLOT


(options, args) = parser.parse_args()

optDoData  = not options.noData
predType   = options.predType
region     = options.region
Type       = options.Type
SavePlot   = options.SavePlot
mcSet      = options.mcSet
LumiProj   = options.LumiProj
OutputDir  = options.outputDir if options.outputDir.startswith("/") else os.path.join(PersonalInfo.personalFolder, options.outputDir)
Analysis   = options.Analysis
year       = options.year

OutputDir  = os.path.join(OutputDir, Analysis, region, "")  # Last "" ensures a trailing '/' is appended to path

tdrstyle.setTDRStyle()

ROOT.gROOT.SetBatch(True)

# VarInfo_zz = {"Mass":["m_{4l} [GeV]","m_{4\ell}",10],"Mjj":["m_{jj} [GeV]","m_{JJ}",20],"Z1Mass":["Z1 Mass","m_{2\ell}",10,],"Z2Mass":["Z2 Mass","m_{2\ell}",10,],"Z1lep0_sip":["Z1 lep 0 Sip","Sip",4],"Z1lep0_iso":["Z1 lep 0 Iso","Iso",4],"Z0lep0_pt":["Z1 lep 0 pT","p_{T}",4],"nJets":["N_{jets} (|#eta^{jet}| < 4.7)","N_{jets} (|#eta^{jet}| < 4.7)",1],"nJets_central":["N_{jets} (|#eta^{jet}| < 4.7)","N_{jets} (|#eta^{jet}| < 4.7)",1],"z":["z1","z1",1],"PtJet1":["p_{T}^{jet1} [GeV]","p_{T}^{jet}",1],"EtaJet1":["#eta^{jet1}","#eta^{jet}",9],"PtJet2":["p_{T}^{jet2} [GeV]","p_{T}^{jet}",1],"EtaJet2":["#eta^{jet2}","#eta^{jet}",10],"Z1pt":["Z1 p_{T}","p_{T}",20],"Z2pt":["Z2 p_{T}","p_{T}",10],"Z1z":["Z1 z","z_{Z_{1}}",7],"Z2z":["Z2 z","z_{Z_{2}}",7],"ptJRatio":["","#Sigma p_{T}/# Sum  ",2],"ptRatio":["","#Sum p_{T}",2],"PtZZ":["p_{T}^{4\\ell}","Sum p_{T}",20],"deltaEtaJJ":["|#eta_{jj}|","|#eta_{jj}|",2],"Dphi":["#Delta #phi_{jj}","#Delta #phi_{jj}",10],"Deta":["|#Delta#eta_{jj}|","#Delta #eta_{jj}",5],"Mjj_Central":["m_{jj}","m_{jj}",20],"Deta_Central":["#Delta #eta_{jj}","#Delta #eta_{jj}",5],"Deta2Jet":["#Delta #eta_{jj}, 2 jet","#Delta #eta_{jj} =2 jet",5],"Deta_noCentral":["#Delta #eta_{jj}, >2 jet","#Delta #eta_{jj} > 2 jet",5],"Deta_1noCentral":["#Delta #eta_{jj}, >2 jet","#Delta #eta_{jj} > 2 jet",5],"PtJet1_noCentral":["#eta Jet","#eta^{jet}",9],"EtaJet1_noCentral":["#eta Jet","#eta^{jet}",10]}

# VarInfo_vbs = {"Mass":["m_{4\ell}","m_{4\ell}",40],"Mjj":["m_{jj}","m_{JJ}",20],"Z1Mass":["Z1 Mass","m_{2\ell}",10,],"Z2Mass":["Z2 Mass","m_{2\ell}",10,],"Z1lep0_sip":["Z1 lep 0 Sip","Sip",4],"Z1lep0_iso":["Z1 lep 0 Iso","Iso",4],"Z0lep0_pt":["Z1 lep 0 pT","p_{T}",4],"nJets":["# jets","# jets",1],"nJets_central":["# jets","# jets",1],"z":["z1","z1",1],"PtJet1":["pT Jet","p_{T}^{jet}",10],"EtaJet1":["#eta Jet","#eta^{jet}",10],"PtJet2":["pT Jet","p_{T}^{jet}",10],"EtaJet2":["#eta Jet","#eta^{jet}",10],"Z1pt":["Z1 p_{T}","p_{T}",20],"Z2pt":["Z2 p_{T}","p_{T}",10],"Z1z":["Z1 z","z_{Z_{1}}",7],"Z2z":["Z2 z","z_{Z_{2}}",7],"ptJRatio":["","#Sigma p_{T}/# Sum  ",2],"ptRatio":["","#Sum p_{T}",2],"PtZZ":["p_{T}^{4\\ell}","Sum p_{T}",60],"deltaEtaJJ":["|#eta_{jj}|","|#eta_{jj}|",2],"Dphi":["#Delta #phi_{jj}","#Delta #phi_{jj}",10],"Deta":["#Delta #eta_{jj}","#Delta #eta_{jj}",5],"Mjj_Central":["m_{jj}","m_{jj}",20],"Deta_Central":["#Delta #eta_{jj}","#Delta #eta_{jj}",5]}

VarInfo_vvx = {
    "AAA_cuts"  : {'title':'Cuts', 'unblind':True, 'logy':False},
    'channel_lep':{'title':'lepton flavour', 'unblind':True}
}

if region in ['SR4P', 'CR3P1F', 'CR2P2F']:
    rebin_mZZG = 1
    if(region in ['SR4P', 'CR3P1F']):
        rebin_mZZG = 2
    elif(region == 'CR2P2F'):
        rebin_mZZG = 2
    VarInfo_vvx.update({
        'ZZ_mass' : {'title':'m_{4\ell}'     , 'rebin':1, 'unblind':True},
        'Z0_mass' : {'title':'m_{Z0}'        , 'rebin':1, 'unblind':True},
        'Z1_mass' : {'title':'m_{Z1}'        , 'rebin':1, 'unblind':True},
        'ZZ_pt'   : {'title':'p_{T}^{Z1}'    , 'rebin':1, 'unblind':True},
        'Z0_l0_pt': {'title':'p_{T}^{Z0, l0}', 'rebin':1, 'unblind':True},
        'Z0_l1_pt': {'title':'p_{T}^{Z0, l1}', 'rebin':1, 'unblind':True},
        'Z1_l0_pt': {'title':'p_{T}^{Z1, l0}', 'rebin':1, 'unblind':True},
        'Z1_l1_pt': {'title':'p_{T}^{Z1, l1}', 'rebin':1, 'unblind':True},
        'PhFRClosure_PASS_mZZG': {'title':'m_{ZZ#gamma} [GeV/c^{2}]', 'unblind':False, 'rebin':rebin_mZZG}
    })
    for name, title in [('4e','4e'), ('2e2m', '2e2\mu'), ('4m', '4\mu')]:
        VarInfo_vvx.update({
            "ZZ_mass_"+name : {'title':"m_{%s}"     %(title), 'rebin':1, 'unblind':True},
            "ZZ_pt_"  +name : {'title':"p_{T}^{%s}" %(title), 'rebin':1, 'unblind':True},
        })
    for name, title in [('ZZ', '4\ell'), ('ZZG', '4\ell\gamma')]:
        VarInfo_vvx.update({
            # name+'_mass_noG'   : {'title':'m_{%s}\:,\ no\:\gamma'                  %(title), 'rebin':1, 'unblind':True },
            name+'_mass_kinG'  : {'title':'m_{%s}\:,\ \gamma\:kin'                 %(title), 'rebin':1, 'unblind':True }, # Just kinematic selection + pixelSeed + electron veto
            name+'_mass_failG' : {'title':'m_{%s}\:,\ \gamma\:loose\,\land\:!tight'%(title), 'rebin':3, 'unblind':True },
            name+'_mass_looseG': {'title':'m_{%s}\:,\ \gamma\:loose'               %(title), 'rebin':1, 'unblind':True }, # Loose = pass 3 cuts
            name+'_mass_tightG': {'title':'m_{%s}\:,\ \gamma\:tight'               %(title), 'rebin':3, 'unblind':False}  # Tight = cutBasedIDLoose()
        })
        
elif region in ['SR3P', 'CR001', 'CR010', 'CR011', 'CR100', 'CR101', 'CR110', 'CR000']:
    VarInfo_vvx.update({
        "WZ_cutflow": {'title':'Cuts', 'logy':False},
        'ZW_massT': {'title':'mT_{3\ell\\nu}'   , 'rebin':1, 'unblind':True},
        'ZW_pt'   : {'title':'p_{T}^{3\ell\\nu}', 'rebin':1, 'unblind':True},
        'W_l_pt'  : {'title':'p_{t,l10};GeV/c'},
        'lll_mass': {'title':'m_{lll};GeV/c^{2}'},
        'paperSel_ZW_massT' : {'title':'m_{T,3l} [GeV/c^{2}]'},
        'paperSel_Z_mass'   : {'title':'m_{Z} [GeV/c^{2}]'   },
        'paperSel_W_massT'  : {'title':'m_{T,W} [GeV/c^{2}]' },
        'paperSel_ZW_pt'    : {'title':'p_{t,ZW} [GeV/c]'    },
        'paperSel_Z_l0_pt'  : {'title':'p_{t,l00} [GeV/c]'   },
        'paperSel_Z_l1_pt'  : {'title':'p_{t,l01} [GeV/c]'   },
        'paperSel_W_l_pt'   : {'title':'p_{t,l10} [GeV/c]'   },
        'paperSel_W_MET_pt' : {'title':'p_{t,MET} [GeV/c]'   },
        'paperSel_lll_mass' : {'title':'m_{lll} [GeV/c^{2}]' },
        'PhFRClosure_PASS_mWZG': {'title':'m_{WZ#gamma} [GeV/c^{2}]', 'unblind':False}
        # 'debug3L_l1_FRSF': {'title':'FR(l_{1})' }
        # 'debug3L_l2_FRSF': {'title':'FR(l_{2})' }
        # 'debug3L_l3_FRSF': {'title':'FR(l_{3})' }
        # 'debug3L_ZW_FRSF': {'title':'FR(ZW)'    }
    })
    for name, title in [('3e','3e'), ('2e1m', '2e1\mu'), ('2m1e', '2\mu1e'), ('3m', '3\mu')]:
        VarInfo_vvx.update({
            'ZW_massT_'+name : {'title':'m_{%s\\nu}'     %(title), 'rebin':1, 'unblind':True},
            'ZW_pt_'   +name : {'title':'p_{T}^{%s\\nu}' %(title), 'rebin':1, 'unblind':True},
            'W_l_pt_'  +name : {'title':'p_{t,l10};GeV/c'},
            'lll_mass_'+name : {'title':'m_{lll};GeV/c^{2}'}
        })
    for name, title in [('ZW', '3\ell\\nu'), ('ZWG', '3\ell\\nu\gamma')]:
        VarInfo_vvx.update({
            # name+'_massT_noG'   : {'title':'mT_{%s}\:,\ no\:\gamma'                %(title), 'rebin':1, 'unblind':True },
            name+'_massT_kinG'  : {'title':'mT_{%s}\:,\ \gamma\:kin'               %(title), 'rebin':1, 'unblind':True },
            name+'_massT_failG' : {'title':'mT_{%s}\:,\ \gamma\:kin\,\land\:!loose'%(title), 'rebin':1, 'unblind':True },
            name+'_massT_looseG': {'title':'mT_{%s}\:,\ \gamma\:loose'             %(title), 'rebin':1, 'unblind':True },
            name+'_massT_tightG': {'title':'mT_{%s}\:,\ \gamma\:tight'             %(title), 'rebin':1, 'unblind':False}
        })
    for name, title in [('e', 'e'), ('m','\mu')]:
        VarInfo_vvx.update({
            'l3_%s_pt'     %(name): {'title': '3^{rd} %s p_{T}'            %(title)},
            'l3_%s_Iso'    %(name): {'title': '3^{rd} %s combRelIsoFSRCorr'%(title)},
            'l3_%s_pt_MET' %(name): {'title': '3^{rd} %s p_{T}'            %(title)},
            'l3_%s_Iso_MET'%(name): {'title': '3^{rd} %s combRelIsoFSRCorr'%(title)}
            
        })

elif region in ['SR2P', 'SR2P_1L', 'SR2P_1P', 'CR2P_1F']:
    VarInfo_vvx.update({
        'AK4_N'         : {'title':'# AK4'   , 'rebin':1, 'unblind':True},
        'AK4_pt'        : {'title':'p_{T}'   , 'rebin':1, 'unblind':True},
        'AK8_N'         : {'title':'# AK8'   , 'rebin':1, 'unblind':True},
        'AK8_pt'        : {'title':'p_{T}'   , 'rebin':1, 'unblind':True},
        'Z_mass_2e'     : {'title':'m_{2e}'  , 'rebin':1, 'unblind':True},
        'Z_mass_2m'     : {'title':'m_{2\mu}', 'rebin':1, 'unblind':True},
        'Z_mass_noG'    : {'title':'m_{2\ell}\:,\ no\:\gamma'        , 'rebin':1, 'unblind':True },
        'Z_mass_kinG'   : {'title':'m_{2\ell}\:,\ kin\:\gamma'       , 'rebin':1, 'unblind':True },
        'Z_mass_failG'  : {'title':'m_{2\ell}\:,\ kin\,\land\:!loose', 'rebin':1, 'unblind':True },
        'Z_mass_looseG' : {'title':'m_{2\ell}\:,\ \gamma\:loose'     , 'rebin':1, 'unblind':False},
        'ZG_mass_kinG'  : {'title':'m_{2\ell\gamma}\:,\ \gamma\:kin'},
        'ZG_mass_looseG': {'title':'m_{2\ell\gamma}\:,\ \gamma\:loose', 'unblind':False},
        'VToJ_mass'     : {'title': 'm_{J}'},
        'VTojj_mass'    : {'title': 'm_{jj}'},
        'VToJFake_mass' : {'title': 'm_{J} fake'},
        'VTojjFake_mass': {'title': 'm_{jj} fake'}
    })
    for Vhad in ['VToJ', 'VToJFake']:
        for classifier in ['PNet', 'deepAK8', 'deepAK8MD']:
            for discriminant in ['TvsQCD', 'WvsQCD', 'ZvsQCD']:
                VarInfo_vvx.update({
                    '%s_%s_%s'%(Vhad, classifier, discriminant): {'title': discriminant+' '+classifier, 'rebin':2, 'unblind':True}
                })
    
VarInfo_VVGamma = deepcopy(VarInfo_vvx)
VarInfo_VVGamma.update({
    'kinPhotons_cutflow'   : {'title':'cut'      , 'unblind':True, 'logy':True},
    'kinPhotons_Nm1'       : {'title':'N-1 cuts' , 'unblind':True},
    'kinPhotons_MVA'       : {'title':'MVA score', 'unblind':True},
    'kinPhRes_dR'          : {'title':'#DeltaR'  , 'unblind':False},
    # 'noKinPh_all_genPh_N'  : {'title': '# #gamma_{GEN}'     },
    # 'noKinPh_all_genPh_pt' : {'title': '#gamma_{GEN} p_{T}' },
    # 'noKinPh_all_genPh_eta': {'title': '#gamma_{GEN} #eta'  },
    # 'noKinPh_rec_genPh_pt' : {'title': '#gamma_{GEN} p_{T}' },
    # 'noKinPh_rec_genPh_eta': {'title': '#gamma_{GEN} #eta'  }
})
for e in ['EB', 'EE']:
    VarInfo_VVGamma.update({
        'kinPh_sieie_' +e: {'title':'#sigma_{i#etai#eta}', 'rebin':1, 'unblind':True, 'logy':True},
        'kinPh_HoverE_'+e: {'title':'HoverE'             , 'rebin':2, 'unblind':True, 'logy':True},
        'kinPh_chIso_' +e: {'title':'chIso'              , 'rebin':1, 'unblind':True, 'logy':True, 'logx':True},
        'kinPh_neIso_' +e: {'title':'neIso'              , 'rebin':2, 'unblind':True, 'logy':True},
        'kinPh_phIso_' +e: {'title':'phIso'              , 'rebin':2, 'unblind':True, 'logy':True}
    })
# for name in ['kin', 'loose', 'medium', 'tight']:
#     VarInfo_VVGamma.update({
#         'sigmaiEtaiEta_'+name+'Photons': ['#sigma_{i#etai#eta}', 1, True]
#     })
VarInfo_VVGamma.update({
    'ph_eScale_N'  : {'title':'Number of #gamma passing selection', 'rebin':1, 'unblind':True},
    'kinPhotons_ID': {'title':'#gamma ID'                         , 'rebin':1, 'unblind':True}
})
# for name in ['all', 'kin', 'loose']:
#     VarInfo_VVGamma.update({
#         'maxG_minL_DR_'+name: {'title':'max_{#gamma}(min_{l}(#DeltaR(#gamma_{%s}, l))' %(name), 'rebin':1, 'unblind':True},
#         'minL_DR_'     +name: {'title':'min_{l}(#DeltaR(#gamma_{%s}, l)'               %(name), 'rebin':1, 'unblind':True},
#     })


if   Analysis.startswith("ZZ"     ): VarInfo = VarInfo_zz
elif Analysis.startswith("VVX"    ): VarInfo = VarInfo_vvx
elif Analysis.startswith("VVGamma"): VarInfo = VarInfo_VVGamma
else                               : VarInfo = VarInfo_vbs

#change the CMS_lumi variables  (see CMS_lumi.py)                                                                                  
#CMS_lumi.lumi_7TeV = "4.8 fb^{-1}"                               
#CMS_lumi.lumi_8TeV = "19.7 fb^{-1}"                                                                                                                                                                                                          

lumi = round(Lumi/1000.,1)
CMS_lumi.writeExtraText = 1
CMS_lumi.extraText = "Preliminary"
CMS_lumi.lumi_sqrtS = "{0}".format(lumi)+" fb^{-1} (13 TeV)\n"

iPos = 0
if( iPos==0 ): CMS_lumi.relPosX = 0.12
iPeriod = 0

try:
    os.stat(OutputDir)
except OSError as e:
    if(not e.errno == 2): raise e  # 2 = No such file or directory
    os.makedirs(OutputDir)  # mkdir() = mkdir  ;  makedirs() = mkdir -p
        
    
InputDir = os.path.join('results', year, Analysis+'_'+region, '')

if LumiProj!="":  InputDir+=LumiProj+"fbm1_"

if Type == 'all':
    variables = VarInfo.keys()
else:
    variables = [ var for var in VarInfo.keys() if re.search(Type, var) ]  # Allow for regexp to be specified from command line
#print "INFO: variables =", variables

c1 = TCanvas( 'c1', mcSet , 900, 1200 )

for Var in variables:
    c1.Clear()
    DoData = optDoData and (VarInfo[Var].get('unblind', True) or region[:2] != 'SR')
    forcePositive=True
    
    # "Temporary" hack for closure test of photon fake rate
    # if 'PhFRClosure_PASS' in Var:
    #     hMC, leg = plotUtils.GetClosureStack(region, InputDir, Var, VarInfo[Var].get('rebin', 1), forcePositive=False)
    # else:
    (hMC, leg) = plotUtils.GetPredictionsPlot(region, InputDir, Var, predType, mcSet, VarInfo[Var].get('rebin', 1), forcePositive=forcePositive)
    (graphData, histodata) = plotUtils.GetDataPlot(InputDir, Var, region, VarInfo[Var].get('rebin', 1), forcePositive=forcePositive)

    if((not hMC.GetStack()) or (not graphData)):
        continue
    
    YMaxMC = YMax = hMC.GetMaximum()
    
    YMaxData = ROOT.TMath.MaxElement(graphData.GetN(), graphData.GetEYhigh()) + ROOT.TMath.MaxElement(graphData.GetN(), graphData.GetY())
    
    hMCErr = deepcopy(hMC.GetStack().Last())
    YMaxMC = hMCErr.GetBinContent(hMCErr.GetMaximumBin()) + hMCErr.GetBinError(hMCErr.GetMaximumBin())
    
    if YMaxData > YMaxMC and DoData: YMax = YMaxData
    else: YMax = YMaxMC
    
    YMax *= 1.37
    
    HistoData = deepcopy(histodata)
    c1.cd()
    pad1 = ROOT.TPad ('hist', '', 0., 0.22, 1.0, 1.0)#0.35
    pad1.SetTopMargin    (0.10)
    pad1.SetRightMargin  (0.06)#0.10
    pad1.SetLeftMargin   (0.16)
    pad1.SetBottomMargin (1.5) 
    pad1.Draw()
    
    c1.cd()
    
    pad2 = ROOT.TPad ('rat', 'Data/MC ratio', 0., 0.0,  1., 0.28)#0.15
    pad2.SetTopMargin (0.01)
    pad2.SetRightMargin (0.06)#0.10
    pad2.SetLeftMargin (0.16)
    pad2.SetBottomMargin(0.3);
    pad2.Draw()
    if VarInfo[Var].get('logx', False):
        pad1.SetLogx()
        pad2.SetLogx()
    
    pad1.cd()
    
    if DoData: histodata.Divide(hMC.GetStack().Last())
    else:
        temp_xaxis = hMC.GetStack().Last().GetXaxis()
        histodata = ROOT.TH1F( "histodata", "", temp_xaxis.GetNbins(), temp_xaxis.GetBinLowEdge(1), temp_xaxis.GetBinUpEdge(temp_xaxis.GetNbins()) )
    
    histodata.GetYaxis().SetTitle("data/MC")
    histodata.GetYaxis().SetTitleSize(0.12)
    histodata.GetYaxis().SetTitleOffset(0.5)
    
    hMC.SetMaximum(YMax)
    hMC.Draw("hist")
    if VarInfo[Var].get('logy', False):
        pad1.SetLogy()
        hMC.GetHistogram().GetYaxis().SetMoreLogLabels()
        if(YMax < 10000):
            hMC.GetHistogram().GetYaxis().SetNoExponent()
    
    hMC.GetHistogram().GetYaxis().SetTitle("Events")
    hMC.GetHistogram().GetYaxis().SetTitleOffset(1.4)
    hMC.GetHistogram().GetYaxis().SetMaxDigits(4)
    hMC.GetHistogram().GetXaxis().SetLabelSize(0)
    
    hMCErr.SetFillStyle(3005)
    hMCErr.SetMarkerStyle(1)
    hMCErr.SetFillColor(ROOT.kBlack)
    hMCErr.Draw("sameE2")
    leg.AddEntry(hMCErr, "MC Err", "f")
    
    if DoData:
        graphData.SetMarkerStyle(20)
        graphData.SetMarkerSize(.9)
        graphData.Draw("samep")
        leg.AddEntry(graphData, "Data", "lpe")
        #HistoData.Draw("same text")
        HistoData.Draw("same")
    
    if Var == "nJets":
        hMC.GetHistogram().GetXaxis().SetTitle("N_{jets} (|#eta^{jet}| < 4.7)")
        hMC.GetHistogram().GetXaxis().SetBinLabel(1, "0 ")
        hMC.GetHistogram().GetXaxis().SetBinLabel(2, "1 ")
        hMC.GetHistogram().GetXaxis().SetBinLabel(3, "2 ")
        hMC.GetHistogram().GetXaxis().SetBinLabel(4, "3 ")
        hMC.GetHistogram().GetXaxis().SetBinLabel(5, ">3 ")
    
        histodata.GetXaxis().SetTitle("N_{jets} (|#eta^{jet}| < 4.7)")
        histodata.GetXaxis().SetBinLabel(1, "0 ")
        histodata.GetXaxis().SetBinLabel(2, "1 ")
        histodata.GetXaxis().SetBinLabel(3, "2 ")
        histodata.GetXaxis().SetBinLabel(4, "3 ")
        histodata.GetXaxis().SetBinLabel(5, ">3 ")

    x1 = leg.GetX1()
    x2 = leg.GetX2()
    shift = 0.78 - (x1+x2)/2
    leg.SetX1(x1+shift)
    leg.SetX2(x2+shift)
    leg.Draw("same")
    
    CMS_lumi.CMS_lumi(c1, iPeriod, iPos)
    
    
    pad2.cd()
        
    Line = ROOT.TLine(hMC.GetXaxis().GetXmin(), 1, hMC.GetXaxis().GetXmax(), 1) 
    Line.SetLineWidth(2)
    Line.SetLineStyle(7)
    
    # yMax_r = histodata.GetBinContent( histodata.GetMaximumBin()) + histodata.GetBinError(histodata.GetMaximumBin() )
    # yMin_r = histodata.GetBinContent( histodata.GetMinimumBin()) - histodata.GetBinError(histodata.GetMinimumBin() )
    # deltaY = (yMax_r - yMin_r)
    yMax_r = 1.5  # max(min(yMax_r + deltaY*0.1, 2), 1.1)
    yMin_r = 0.5  # min(max(yMin_r - deltaY*0.1, 0), 0.9)

    # hArea = deepcopy(hMC.GetStack().Last())  # in ratio plot, the gray area representing MC error
    # for bin in range(1, hArea.GetNbinsX()+1):
    #     r = hArea.GetBinContent(bin)
    #     if(r == 0):
    #         hArea.SetBinContent(bin, 1)
    #         hArea.SetBinError  (bin, 0)
    #     else:
    #         hArea.SetBinContent(bin, hArea.GetBinContent(bin)/r)
    #         hArea.SetBinError  (bin, hArea.GetBinContent(bin)/r)
    # if (hArea.GetXaxis().GetXmin() > 0.001 and hArea.GetXaxis().GetXmax() < 1000):
    #     hArea.GetXaxis().SetNoExponent()
    # #hArea.GetXaxis().SetMoreLogLabels()
    # hArea.SetFillColor(ROOT.kGray)
    # hArea.Draw("E3")
    
    histodata.GetXaxis().SetTitle(VarInfo[Var]['title'])
    histodata.GetXaxis().SetLabelSize(0.08)
    histodata.GetYaxis().SetLabelSize(0.08)
    histodata.GetXaxis().SetTitleSize(0.08)
    if (histodata.GetXaxis().GetXmin() > 0.001 and histodata.GetXaxis().GetXmax() < 1000):
        histodata.GetXaxis().SetNoExponent()
    histodata.SetMaximum( yMax_r )
    histodata.SetMinimum( yMin_r )
    histodata.SetMarkerStyle(20)
    #histodata.Draw("E")
    histodata.Draw("axis")
    Line.Draw()
    histodata.Draw("E same")
    
    if(not DoData):
        xm, xM = hMC.GetXaxis().GetXmin(), hMC.GetXaxis().GetXmax()
        xstart = (xm + xM)/2 - (xM - xm)/8
        text = ROOT.TText(xstart, 1.1, "BLINDED")
        text.SetTextSize(.12)
        text.Draw("same")
    
    Title=Var+"_"+mcSet #+"_"+region
    
    ROOT.gStyle.SetOptStat(0);   
    ROOT.gStyle.SetOptTitle(0)
    c1.Update()
    
    c1.SetTitle(Title)
    if SavePlot:
        # c1.SaveAs(os.path.join(OutputDir, Title+".root"))
        c1.SaveAs(os.path.join(OutputDir, Title+".png"))
        # c1.SaveAs(os.path.join(OutputDir, Title+".eps"))
        c1.SaveAs(os.path.join(OutputDir, Title+".pdf"))
    
    
    # if SavePlot:
    #     fullDir = "{Dir}/{Analysis}/{Analysis}".format(Dir=Dir, Analysis=Analysis, region=region)
    #     try:
    #         os.stat(fullDir)
    #     except OSError as e:
    #         if(not e.errno == 2): raise e  # 2 = No such file or directory. Let other exceptions pass
    #         os.makedirs(fullDir)
    #         os.system('for d in $(find {Dir} -type d) ; do [ -f $d/index.php ] || echo "cp {Dir}/index.php $d/" ; done'.format(Dir=Dir) )
        
    
    #     if(LumiProj!=""):
    #         c1.SaveAs("~/www/PlotsVV/13TeV_"+LumiProj+"fb/"+Title+".png")
    #         c1.SaveAs("~/www/PlotsVV/13TeV_"+LumiProj+"fb/"+Title+".pdf")
    #     else:
    #         c1.SaveAs(PersonalInfo.PersonalFolder+Dir+"Reco/"+Title+".png")
    #         c1.SaveAs(PersonalInfo.PersonalFolder+Dir+"Reco/"+Title+".root")
    #         c1.SaveAs(PersonalInfo.PersonalFolder+Dir+"Reco/"+Title+".pdf")
    #         c1.SaveAs("~/../../../../../eos/user/g/gpinnaan/"+Dir+"Reco/"+Title+".pdf")
    

