#! /usr/bin/env python2
import sys #,ast
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
from plotUtils import*
import CMS_lumi, tdrstyle
import PersonalInfo
Lumi   = 35900

regions = ['SR4P', 'CR3P1F' , 'CR2P2F' , 'SR4P_1L', 'SR4P_1P', 'CR4P_1F', 'CR4L',    
           'SR3P', 'CR110'  , 'CR101'  , 'CR011'  , 'CR100'  , 'CR001'  , 'CR010', 'CR000', 'SR3P_1L', 'SR3P_1P', 'CR3P_1F', 'CRLFR', 'CR3L',
           'SR2P', 'SR2P_1L', 'SR2P_1P', 'CR2P_1F', 
           'SR_HZZ', 'CR2P2F_HZZ', 'CR3P1F_HZZ', 'CR_HZZ']


parser = OptionParser(usage="usage: %prog <final state> [options]")

parser.add_option("-r", "--region", dest="region",
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

parser.add_option("-m", "--mcset", dest="mcSet",
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

parser.add_option("-A", "--Analysis", dest="Analysis",
                  default="VVXAnalyzer",
                  help="Analysis. Default is ZZ. Other option is VBS")

parser.add_option("-y", "--year", dest="year",
                  default="2016",
                  help= "valid inputs are 2016, 2017, 2018")


#REMEMBER ADD DEFINTION PLOT


(options, args) = parser.parse_args()

optDoData  = not options.noData ## Fixme, ratio plot to be removed
predType   = options.predType
region     = options.region
Type       = options.Type
Save       = options.SavePlot
mcSet      = options.mcSet
LumiProj   = options.LumiProj
OutputDir  = options.outputDir if options.outputDir.startswith("/") else PersonalInfo.personalFolder+'/'+options.outputDir
Analysis   = options.Analysis
year       = options.year

OutputDir  = os.path.join(OutputDir, Analysis, region, "")  # Last "" ensures a trailing '/' is appended to path

tdrstyle.setTDRStyle()

ROOT.gROOT.SetBatch(True)

# InfoType_zz = {"Mass":["m_{4l} [GeV]","m_{4\ell}",10],"Mjj":["m_{jj} [GeV]","m_{JJ}",20],"Z1Mass":["Z1 Mass","m_{2\ell}",10,],"Z2Mass":["Z2 Mass","m_{2\ell}",10,],"Z1lep0_sip":["Z1 lep 0 Sip","Sip",4],"Z1lep0_iso":["Z1 lep 0 Iso","Iso",4],"Z0lep0_pt":["Z1 lep 0 pT","p_{T}",4],"nJets":["N_{jets} (|#eta^{jet}| < 4.7)","N_{jets} (|#eta^{jet}| < 4.7)",1],"nJets_central":["N_{jets} (|#eta^{jet}| < 4.7)","N_{jets} (|#eta^{jet}| < 4.7)",1],"z":["z1","z1",1],"PtJet1":["p_{T}^{jet1} [GeV]","p_{T}^{jet}",1],"EtaJet1":["#eta^{jet1}","#eta^{jet}",9],"PtJet2":["p_{T}^{jet2} [GeV]","p_{T}^{jet}",1],"EtaJet2":["#eta^{jet2}","#eta^{jet}",10],"Z1pt":["Z1 p_{T}","p_{T}",20],"Z2pt":["Z2 p_{T}","p_{T}",10],"Z1z":["Z1 z","z_{Z_{1}}",7],"Z2z":["Z2 z","z_{Z_{2}}",7],"ptJRatio":["","#Sigma p_{T}/# Sum  ",2],"ptRatio":["","#Sum p_{T}",2],"PtZZ":["p_{T}^{4\\ell}","Sum p_{T}",20],"deltaEtaJJ":["|#eta_{jj}|","|#eta_{jj}|",2],"Dphi":["#Delta #phi_{jj}","#Delta #phi_{jj}",10],"Deta":["|#Delta#eta_{jj}|","#Delta #eta_{jj}",5],"Mjj_Central":["m_{jj}","m_{jj}",20],"Deta_Central":["#Delta #eta_{jj}","#Delta #eta_{jj}",5],"Deta2Jet":["#Delta #eta_{jj}, 2 jet","#Delta #eta_{jj} =2 jet",5],"Deta_noCentral":["#Delta #eta_{jj}, >2 jet","#Delta #eta_{jj} > 2 jet",5],"Deta_1noCentral":["#Delta #eta_{jj}, >2 jet","#Delta #eta_{jj} > 2 jet",5],"PtJet1_noCentral":["#eta Jet","#eta^{jet}",9],"EtaJet1_noCentral":["#eta Jet","#eta^{jet}",10]}

# InfoType_vbs = {"Mass":["m_{4\ell}","m_{4\ell}",40],"Mjj":["m_{jj}","m_{JJ}",20],"Z1Mass":["Z1 Mass","m_{2\ell}",10,],"Z2Mass":["Z2 Mass","m_{2\ell}",10,],"Z1lep0_sip":["Z1 lep 0 Sip","Sip",4],"Z1lep0_iso":["Z1 lep 0 Iso","Iso",4],"Z0lep0_pt":["Z1 lep 0 pT","p_{T}",4],"nJets":["# jets","# jets",1],"nJets_central":["# jets","# jets",1],"z":["z1","z1",1],"PtJet1":["pT Jet","p_{T}^{jet}",10],"EtaJet1":["#eta Jet","#eta^{jet}",10],"PtJet2":["pT Jet","p_{T}^{jet}",10],"EtaJet2":["#eta Jet","#eta^{jet}",10],"Z1pt":["Z1 p_{T}","p_{T}",20],"Z2pt":["Z2 p_{T}","p_{T}",10],"Z1z":["Z1 z","z_{Z_{1}}",7],"Z2z":["Z2 z","z_{Z_{2}}",7],"ptJRatio":["","#Sigma p_{T}/# Sum  ",2],"ptRatio":["","#Sum p_{T}",2],"PtZZ":["p_{T}^{4\\ell}","Sum p_{T}",60],"deltaEtaJJ":["|#eta_{jj}|","|#eta_{jj}|",2],"Dphi":["#Delta #phi_{jj}","#Delta #phi_{jj}",10],"Deta":["#Delta #eta_{jj}","#Delta #eta_{jj}",5],"Mjj_Central":["m_{jj}","m_{jj}",20],"Deta_Central":["#Delta #eta_{jj}","#Delta #eta_{jj}",5]}

InfoType_vvx = {
    "AAA_cuts"        : ["Cuts", 1, True]
    # "looseG_pt": ["p_{T}^\gamma", 1, True],
    # "failG_pt" : ["p_{T}^\gamma", 1, True],
    # "tightG_pt": ["p_{T}^\gamma", 1, True]
    # ,
}

if region in ['SR4P', 'CR3P1F', 'CR2P2F']:
    InfoType_vvx.update({
        "ZZ_mass" : ["m_{4\ell}"     , 1, True],
        "Z0_mass" : ["m_{ZZ}"        , 1, True],
        "Z1_mass" : ["m_{Z1}"        , 1, True],
        "ZZ_pt"   : ["p_{T}^{Z1}"    , 1, True],
        "Z0_l0_pt": ["p_{T}^{Z0, l0}", 1, True],
        "Z0_l1_pt": ["p_{T}^{Z0, l1}", 1, True],
        "Z1_l0_pt": ["p_{T}^{Z1, l0}", 1, True],
        "Z1_l1_pt": ["p_{T}^{Z1, l1}", 1, True]
    })
    for name, title in [('4e','4e'), ('2e2m', '2e2\mu'), ('4m', '4\mu')]:
        InfoType_vvx.update({
            "ZZ_mass_"+name : ["m_{%s}"     %(title), 1, True],
            "ZZ_pt_"  +name : ["p_{T}^{%s}" %(title), 1, True],
        })
    for name, title in [('ZZ', '4\ell'), ('ZZG', '4\ell\gamma')]:
        InfoType_vvx.update({
            name+"_mass_noG"   : ["m_{%s}\:,\ no\:\gamma"                %(title), 1, True],
            name+"_mass_kinG"  : ["m_{%s}\:,\ \gamma\:kin"               %(title), 1, True],
            name+"_mass_failG" : ["m_{%s}\:,\ \gamma\:kin\,\land\:!loose"%(title), 1, True],
            name+"_mass_looseG": ["m_{%s}\:,\ \gamma\:loose"             %(title), 1, False]
        })
        
elif region in ['SR3P', 'CR001', 'CR010', 'CR011', 'CR100', 'CR101', 'CR110']:
    InfoType_vvx.update({
        "ZW_massT": ["mT_{3\ell\\nu}"   , 1, True],
        "ZW_pt"   : ["p_{T}^{3\ell\\nu}", 1, True],
    })
    for name, title in [('3e','3e'), ('2e1m', '2e1\mu'), ('2m1e', '2\mu1e'), ('3m', '3\mu')]:
        InfoType_vvx.update({
            "ZW_massT_"+name : ["m_{%s\\nu}"     %(title), 1, True],
            "ZW_pt_"   +name : ["p_{T}^{%s\\nu}" %(title), 1, True],
        })
    for name, title in [('ZW', '3\ell\\nu'), ('ZWG', '3\ell\\nu\gamma')]:
        InfoType_vvx.update({
            name+"_massT_noG"   : ["mT_{%s}\:,\ no\:\gamma"                %(title), 1, True],
            name+"_massT_kinG"  : ["mT_{%s}\:,\ \gamma\:kin"               %(title), 1, True],
            name+"_massT_failG" : ["mT_{%s}\:,\ \gamma\:kin\,\land\:!loose"%(title), 1, True],
            name+"_massT_looseG": ["mT_{%s}\:,\ \gamma\:loose"             %(title), 1, False]
        })

# elif region in ['SR2L']:
#     InfoType_vvx.update({
#         "Z2l_mass"       : ["m_{2\ell}"              , 1, True],
#         "Z2l_mass_noG"   : ["m_{2\ell} no \gamma"    , 1, True],
#         "Z2l_mass_kinG"  : ["m_{2\ell} kin \gamma"   , 1, True],
#         "Z2l_mass_failG" : ["m_{2\ell} kin && !loose", 1, True],
#         "Z2l_mass_looseG": ["m_{2\ell} loose"        , 1, False]
#     })

InfoType_VVGamma = deepcopy(InfoType_vvx)
for name in ["kin", "loose", "medium", "tight"]:
    InfoType_VVGamma.update({
        "sigmaiEtaiEta_"+name+"Photons": ["#sigma_{i#etai#eta}", 1, True]
})
InfoType_VVGamma.update({
    "ph_eScale_count" : ["Number of #gamma passing selection", 1, True],
    "kinPhotons_ID": ["#gamma ID", 1, True]
})
for name in ['all', 'kin', 'loose']:
    InfoType_VVGamma.update({
        "maxG_minL_DR_"+name: ["max_{#gamma}(min_{l}(#DeltaR(#gamma_{%s}, l))" %(name), 1, True],
        "minL_DR_"     +name: ["min_{l}(#DeltaR(#gamma_{%s}, l)"               %(name), 1, True],
})


if   Analysis == "ZZ"             : InfoType = InfoType_zz
elif Analysis == "VVXAnalyzer"    : InfoType = InfoType_vvx
elif Analysis == "VVGammaAnalyzer": InfoType = InfoType_VVGamma
else                              : InfoType = InfoType_vbs

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
        
    
InputDir = "results/"+year+"/"+Analysis+"_"+region+"/"

if LumiProj!="":  InputDir+=LumiProj+"fbm1_"

if Type == 'all':
    variables = InfoType.keys()
else:
    variables = [ var for var in InfoType.keys() if re.search(Type, var) ]  # Allow for regexp to be specified from command line

c1 = TCanvas( 'c1', mcSet , 900, 1200 )

for Var in variables:
    c1.Clear()
    DoData = optDoData and InfoType[Var][-1]
    
    (hMC, leg) = GetPredictionsPlot(region, InputDir, Var, predType, mcSet, InfoType[Var][-2])
    (hData, histodata) = GetDataPlot(InputDir, Var, region, InfoType[Var][-2])
    
    # if(any(s in Var for s in [""]))
    if((not hMC.GetStack()) or (not hData)):
        continue
    
    YMaxMC = YMax = hMC.GetMaximum()
    
    YMaxData = ROOT.TMath.MaxElement(hData.GetN(), hData.GetEYhigh()) + ROOT.TMath.MaxElement(hData.GetN(), hData.GetY())
    
    hMCErr = copy.deepcopy(hMC.GetStack().Last())
    YMaxMC = hMCErr.GetBinContent(hMCErr.GetMaximumBin()) + hMCErr.GetBinError(hMCErr.GetMaximumBin())
    
    if YMaxData > YMaxMC and DoData: YMax = YMaxData
    else: YMax = YMaxMC
    
    YMax *= 1.37
    
    HistoData = copy.deepcopy(histodata)
    c1.cd()
    pad1 = ROOT.TPad ('hist', '', 0., 0.22, 1.0, 1.0)#0.35
    pad1.SetTopMargin    (0.10)
    pad1.SetRightMargin  (0.06)#0.10
    pad1.SetLeftMargin   (0.16)
    pad1.SetBottomMargin (1.5) 
    #pad1.SetLogy()
    pad1.Draw()
        
    c1.cd()
    
    pad2 = ROOT.TPad ('rat', 'Data/MC ratio', 0., 0.0,  1., 0.28)#0.15
    pad2.SetTopMargin (0.01)
    pad2.SetRightMargin (0.06)#0.10
    pad2.SetLeftMargin (0.16)
    pad2.SetBottomMargin(0.3);
    pad2.Draw()
    
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
        hData.SetMarkerStyle(20)
        hData.SetMarkerSize(.9)
        hData.Draw("samep")
        leg.AddEntry(hData, "Data", "lpe")
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
    
    
    leg.Draw("same")    
    
    CMS_lumi.CMS_lumi(c1, iPeriod, iPos)
    
    
    pad2.cd()
        
    Line = ROOT.TLine(hMC.GetXaxis().GetXmin(), 1, hMC.GetXaxis().GetXmax(), 1) 
    Line.SetLineWidth(2)
    histodata.GetXaxis().SetTitle(InfoType[Var][0])
    histodata.GetXaxis().SetLabelSize(0.08)
    histodata.GetYaxis().SetLabelSize(0.08)
    histodata.GetXaxis().SetTitleSize(0.08)
    
    # yMax_r = histodata.GetBinContent( histodata.GetMaximumBin()) + histodata.GetBinError(histodata.GetMaximumBin() )
    # yMin_r = histodata.GetBinContent( histodata.GetMinimumBin()) - histodata.GetBinError(histodata.GetMinimumBin() )
    # deltaY = (yMax_r - yMin_r)
    yMax_r = 1.5  # max(min(yMax_r + deltaY*0.1, 2), 1.1)
    yMin_r = 0.5  # min(max(yMin_r - deltaY*0.1, 0), 0.9)
    histodata.SetMaximum( yMax_r )
    histodata.SetMinimum( yMin_r )
    
    histodata.SetMarkerStyle(20)
    histodata.Draw("E1")
    Line.Draw("same")
    
    Title=Var+"_"+mcSet+"_"+region
    
    ROOT.gStyle.SetOptStat(0);   
    ROOT.gStyle.SetOptTitle(0)
    c1.Update()
    
    c1.SetTitle(Title)
    #c1.SaveAs(OutputDir + Title+".root")
    c1.SaveAs(OutputDir + Title+".png")
    #c1.SaveAs(OutputDir + Title+".eps")
    c1.SaveAs(OutputDir + Title+".pdf")
    
    
    # if Save:
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
    

