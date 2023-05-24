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
from argparse import ArgumentParser
import ROOT
ROOT.PyConfig.IgnoreCommandLineOptions = True
import CrossInfo
from CrossInfo import* 
from ROOT import TH1F,TCanvas, TLegend
import plotUtils  # GetPredictionsPlot, GetDataPlot
from variablesInfo import getVariablesInfo
import CMS_lumi, tdrstyle
import PersonalInfo
from Colours import Evidence

lumi_dict = {'2016': 35900, '2017': 41500, '2018': 59700}
regions = ['SR4P', 'CR3P1F' , 'CR2P2F' , 'SR4P_1L', 'SR4P_1P', 'CR4P_1F', 'CR4L',    
           'SR3P', 'CR110'  , 'CR101'  , 'CR011'  , 'CR100'  , 'CR001'  , 'CR010', 'CR000', 'SR3P_1L', 'SR3P_1P', 'CR3P_1F', 'CRLFR', 'CR3L',
           'SR2P', 'SR2P_1L', 'SR2P_1P', 'CR2P_1F', 
           'SR_HZZ', 'CR2P2F_HZZ', 'CR3P1F_HZZ', 'CR_HZZ']

parser = ArgumentParser()

parser.add_argument("-r", "--region", dest="region", choices=regions,
                  default="SR4P",
                  help="Available regions are {0:s}. Default is SR4P.".format(', '.join(regions)))

parser.add_argument("-f", "--finalstate", dest="FinalState",
                  default="4l",
                  help="Final state are 4l, 4m, 2e2m and 4e. Default is 4l")

parser.add_argument("--nodata", dest="doData",
                  action="store_false",
                  default=True,
                  help="Forces to NOT draw data on every plot")

parser.add_argument("-u", "--unblind", dest='unblind', action="store_true", default=False, help="Unblinds plots marked as blinded")

parser.add_argument("-t", "--type", dest="Type",
                  default="all",
                  help= "type type to choose the  plot you want. Mass, Jets, DeltaEta, mjj")

parser.add_argument("-s", "--skip", dest="Skip",
                  default=None,
                  help= "Plots names that match this regex will be skipped")

parser.add_argument("-m", "--mcset", dest="mcSet", choices=['mad', 'pow'],
                  default="mad",
                  help= "Monte Carlo Set, pow for Powheg, mad for amcatnlo")

parser.add_argument("-p", "--prediction-type", dest="predType",
                  default="fromCR",
                  help= "Type of prediction. fromCR = non-prompt leptons from CRs, rare background from MC; fullMC = all from MC; fakeMC = use MC in CRs instead of data. Default is fromCR")

parser.add_argument("-l", "--lumiProj", dest="LumiProj",
                  default="",
                  help="Lumi projection")

parser.add_argument("-o", "--outputDir", dest="outputDir",
                  default="last",
                  help="Directory where save plots")

parser.add_argument("-A", "--Analysis", dest="Analysis", choices=['VVXAnalyzer', 'VVGammaAnalyzer', 'ZZAnalyzer'],
                  default="VVXAnalyzer",
                  help="Analysis. Default is VVX; other options are ZZ, VBS and VVGamma")

parser.add_argument("-y", "--year", dest="year",
                  default="2016",
                  help= "valid inputs are 2016, 2017, 2018")

parser.add_argument("-v", "--verbose", dest="verbosity",
                    action="count", default=1,
                    help="Increase verbosity")
parser.add_argument("--verbosity", type=int,
                    help="Set verbosity")
parser.add_argument("-q", "--quiet", dest="verbosity",
                    action="store_const", const=0,
                    help="Set verbose to minimum")

parser.add_argument("-i", "--inputDir",
                    default="results",
                    help="Directory containing the input rootfiles")

#REMEMBER ADD DEFINTION PLOT

ROOT.gInterpreter.Declare(
'''\
#import "TMath.h"
void drawtext(const char* graphName, const char* format="%.4g")
{
    Int_t i, n;
    Double_t x, y, xm1, xp1;
    TLatex *l;
    TGraph *g = (TGraph*)gPad->GetListOfPrimitives()->FindObject(graphName);
    if(!g){
        printf(">>> drawtext(%s, %s): graph not found!\\n", graphName, format);
        return;
    }

    Double_t dummy, xmin, ymin, xmax, ymax;
    gPad->GetRangeAxis(dummy, ymin, dummy, ymax);
    g->ComputeRange(xmin, dummy, xmax, dummy);
    Double_t factor = TMath::Exp((ymax/ymin)*0.08);
    Double_t step   = (ymax-ymin)*0.03;
    // printf(">>> xmin: %f, xmax: %f, ymin: %.4e, ymax: %.4e", xmin, xmax, ymin, ymax);
    // if(gPad->GetLogy())
    //     printf(" --> ratio: %.4f , factor: %.4f\\n", ymax/ymin, factor);
    // else
    //     printf(" --> diff : %.4f , step  : %.4f\\n", ymax-ymin, step  );

    n = g->GetN();
    g->GetPoint(1,xp1,y);
    for (i=0; i<n; i++) {
        g->GetPoint(i,x,y);
        Double_t xwidth = (i == 0 ? xp1-x : x-xm1);
        Double_t xposition = x - xwidth/3;
        Double_t yposition = gPad->GetLogy() ? y*factor : y+step;
        //printf("\\t%2d: y: %.4g, yposition: %.4g\\n", i, y, yposition);
        //printf("\\t%2d: x: %f, xm1: %f, width: %f\\n", i, x, xm1, xwidth);
        l = new TLatex(xposition, yposition, Form(format,y));
        l->SetTextSize(0.02 + 0.01/n);
        //l->SetTextFont(42);
        //l->SetTextAlign(21);
        l->Paint();
        xm1 = x;
    }
}
'''
)

options = parser.parse_args()

if(not options.doData and options.unblind):
    parser.error('"--nodata" and "--unblind" are mutually exclusive.')

optDoData  = options.doData
predType   = options.predType
region     = options.region
Type       = options.Type
mcSet      = options.mcSet
LumiProj   = options.LumiProj
OutputDir  = options.outputDir if options.outputDir.startswith("/") else os.path.join(PersonalInfo.personalFolder, options.outputDir)
Analysis   = options.Analysis
year       = options.year

InputDir   = os.path.join(options.inputDir, year, Analysis+'_'+region, '')
OutputDir  = os.path.join(OutputDir, Analysis, year, predType, region, "")  # Last "" ensures a trailing '/' is appended to path
try:
    os.stat(OutputDir)
except OSError as e:
    if(not e.errno == 2): raise e  # 2 = No such file or directory
    os.makedirs(OutputDir)  # mkdir() = mkdir  ;  makedirs() = mkdir -p


tdrstyle.setTDRStyle()
ROOT.gROOT.SetBatch(True)

if LumiProj != "":
    InputDir+=LumiProj+"fbm1_"
    lumi = LumiProj
else:
    lumi = lumi_dict[year]
lumi = round(lumi/1000.,1)
CMS_lumi.writeExtraText = 1
CMS_lumi.extraText = "Preliminary"
CMS_lumi.lumi_sqrtS = "{0} fb^{{-1}} (13 TeV)\n".format(lumi)

iPos = 0
if( iPos==0 ): CMS_lumi.relPosX = 0.12
iPeriod = 0


VarInfo = getVariablesInfo(Analysis, region)

if Type == 'all':
    variables = VarInfo.keys()
else:
    variables = [ var for var in VarInfo.keys() if re.search(Type, var) ]  # Allow for regexp to be specified from command line
    if len(variables) == 0:
        print 'WARN: no variables matching regex "{}" for {} in {}'.format(Type, Analysis, region)
        exit(0)

if options.Skip is not None:
    variables = [ var for var in variables if not re.search(options.Skip, var) ]
    if len(variables) == 0:
        print 'WARN: using regex "{}" all variables are skipped'
        exit(0)

if(options.verbosity >= 2):
    print 'INFO: variables =', variables
variables.sort()

c1 = TCanvas( 'c1', mcSet , 900, 1200 )

for Var in variables:
    info = VarInfo[Var]
    info.update({'name':Var})
    c1.Clear()
    DoData = optDoData and (info.get('unblind', True) or options.unblind or region[:2] != 'SR')

    forcePositive=True
    
    # "Temporary" hack for closure test of photon fake rate
    if 'PhFRClosure' in Var and 'PASS' in Var:
        hMC, leg = plotUtils.GetClosureStack(region, InputDir, info, forcePositive=False, verbosity=options.verbosity)
    else:
        (hMC, leg) = plotUtils.GetPredictionsPlot(region, InputDir, info, predType, mcSet, forcePositive=forcePositive, verbosity=options.verbosity)

    if(not hMC.GetStack()):
        print Evidence('ERROR'), 'skipping', Var, 'because: no MC'
        continue

    if(DoData):
        (graphData, histodata) = plotUtils.GetDataPlot(region, InputDir, info, forcePositive=forcePositive, verbosity=options.verbosity)

        if(not (graphData and histodata)):
            print Evidence('ERROR'), 'skipping', Var, 'because: no data'
            continue
    
    hMCErr = deepcopy(hMC.GetStack().Last())
    
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
    if info.get('logx', False):
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
    
    # Maximum and minimum of upper plot
    YMax = info.get('ymax', False)
    if(not YMax):
        YMaxData = ROOT.TMath.MaxElement(graphData.GetN(), graphData.GetEYhigh()) + ROOT.TMath.MaxElement(graphData.GetN(), graphData.GetY()) if DoData else 0.
        YMaxMC = hMCErr.GetBinContent(hMCErr.GetMaximumBin()) + hMCErr.GetBinError(hMCErr.GetMaximumBin())
        YMax = max(YMaxMC, YMaxData)
        YMax *= 1.37
        
        if info.get('logy', False):
            YMax *= 10
            pad1.SetLogy()
            hMC.GetHistogram().GetYaxis().SetMoreLogLabels()
            if(YMax < 100000):
                hMC.GetHistogram().GetYaxis().SetNoExponent()
        if('AAA_cuts' in Var):
            YMax *= 1.3

    hMC.SetMaximum(YMax)
    YMin = info.get('ymin', False)
    if(YMin):
        hMC.SetMinimum(YMin)
    
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
        if(info.get('text')):
            texec = ROOT.TExec("texec", 'drawtext("{}");'.format(graphData.GetName()))
            graphData.GetListOfFunctions().Add(texec)
            graphData.Draw("samep text")
        else:
            graphData.Draw("samep")
        leg.AddEntry(graphData, "Data", "lpe")
    
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
    
    if(info.get('title')):
        histodata.GetXaxis().SetTitle(info['title'])
    histodata.GetXaxis().SetLabelSize(0.08)
    histodata.GetYaxis().SetLabelSize(0.08)
    histodata.GetXaxis().SetTitleSize(0.08)
    if (histodata.GetXaxis().GetXmin() > 0.001 and histodata.GetXaxis().GetXmax() < 1000):
        histodata.GetXaxis().SetNoExponent()
    histodata.SetMarkerStyle(20)
    
    # Dummy TH1 to act as a frame for the ratio plot
    frame_ratio = deepcopy(histodata)
    frame_ratio.GetYaxis().SetRangeUser(yMin_r, yMax_r)
    frame_ratio.Draw("axis")
    
    Line.Draw()
    histodata.Draw("E0 same")
    
    
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

    # c1.SaveAs(os.path.join(OutputDir, Title+".root"))
    c1.SaveAs(os.path.join(OutputDir, Title+".png"))
    # c1.SaveAs(os.path.join(OutputDir, Title+".eps"))
    c1.SaveAs(os.path.join(OutputDir, Title+".pdf"))
    
    del histodata

