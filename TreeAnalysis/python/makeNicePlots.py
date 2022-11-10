#! /usr/bin/env python2

######################################################################################################################################################
# Data/MC comparison with nice style                                                                                                                 #
#                                                                                                                                                    #
# example command:                                                                                                                                   #
#   for region in SR3P ; do ./python/makeNicePlots.py -S -A VVGammaAnalyzer -t all -m pow -p fullMC -o last/EXT_fullMC -y 2016 -r $region ; done     #
#                                                                                                                                                    #
# Authors: A. Mecca, G. Pinnan (?)                                                                                                                   #
######################################################################################################################################################

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
from variablesInfo import getVariablesInfo
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
                  help="Available regions are {0:s}. Default is SR4P.".format(', '.join(regions)))

parser.add_option("-f", "--finalstate", dest="FinalState",
                  default="4l",
                  help="Final state are 4l, 4m, 2e2m and 4e. Default is 4l")

parser.add_option("--nodata", dest="noData",
                  action="store_true",
                  default=False,
                  help="Forces to NOT draw data on every plot")

parser.add_option("-u", "--unblind", dest='Unblind', action="store_true", default=False, help="Unblinds plots marked as blinded")

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

if(options.noData and options.Unblind):
    parser.error('"--nodata" and "--unblind" are mutually exclusive.')

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

VarInfo = getVariablesInfo(Analysis, region)

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
# print "INFO: variables =", variables
# exit(0)
variables.sort()

c1 = TCanvas( 'c1', mcSet , 900, 1200 )

for Var in variables:
    c1.Clear()
    DoData = optDoData and (VarInfo[Var].get('unblind', True) or region[:2] != 'SR')
    print ">>> DoData =", DoData, ", optDoData =", optDoData, ", VarInfo[Var].get('unblind', True) =", VarInfo[Var].get('unblind', True), ", region[:2] != 'SR'", region[:2] != 'SR'
    forcePositive=True
    
    # "Temporary" hack for closure test of photon fake rate
    # if 'PhFRClosure_PASS' in Var:
    #     hMC, leg = plotUtils.GetClosureStack(region, InputDir, Var, VarInfo[Var].get('rebin', 1), forcePositive=False)
    # else:
    (hMC, leg) = plotUtils.GetPredictionsPlot(region, InputDir, Var, predType, mcSet, VarInfo[Var].get('rebin', 1), forcePositive=forcePositive)
    (graphData, histodata) = plotUtils.GetDataPlot(InputDir, Var, region            , VarInfo[Var].get('rebin', 1), forcePositive=forcePositive)

    if( (not hMC.GetStack()) or (DoData and (not graphData)) ):
        continue
    
    
    if(DoData):
        YMaxData = ROOT.TMath.MaxElement(graphData.GetN(), graphData.GetEYhigh()) + ROOT.TMath.MaxElement(graphData.GetN(), graphData.GetY())
    else:
        YMaxData = 0.
    
    hMCErr = deepcopy(hMC.GetStack().Last())
    YMaxMC = hMCErr.GetBinContent(hMCErr.GetMaximumBin()) + hMCErr.GetBinError(hMCErr.GetMaximumBin())
    
    YMax = max(YMaxMC, YMaxData)
    
    YMax *= 1.37
    
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
    
    if(Var == 'AK8_pt'):
        for h in hMC.GetStack():
            h.GetXaxis().SetRangeUser(160., 500.)
        histodata.GetXaxis().SetRangeUser(160., 500.)
    
    hMC.Draw("hist")
    
    if('AAA_cuts' in Var):
        hMC.GetXaxis().SetTickLength(0.)
        YMax *= 1.3
    
    if VarInfo[Var].get('logy', False):
        YMax *= 10
        pad1.SetLogy()
        hMC.GetHistogram().GetYaxis().SetMoreLogLabels()
        if(YMax < 10000):
            hMC.GetHistogram().GetYaxis().SetNoExponent()
    hMC.SetMaximum(YMax)
    
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
        print ">>> INFO: draw data"
        leg.AddEntry(graphData, "Data", "lpe")
        #histodata.Draw("same text")
        histodata.Draw("same")
    
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
        
    Line = ROOT.TLine(hMC.GetXaxis().GetBinLowEdge(hMC.GetXaxis().GetFirst()), 1, hMC.GetXaxis().GetXmax(), 1) 
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
    

