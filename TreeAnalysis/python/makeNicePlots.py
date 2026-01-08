#! /usr/bin/env python2

######################################################################################################################################################
# Data/MC comparison with nice style                                                                                                                 #
#                                                                                                                                                    #
# example command:                                                                                                                                   #
#   for region in SR3P ; do ./python/makeNicePlots.py -S -A VVGammaAnalyzer -t all -m pow -p fullMC -o last/EXT_fullMC -y 2016 -r $region ; done     #
#                                                                                                                                                    #
# Authors: A. Mecca, G. L. Pinna Angioni (?)                                                                                                         #
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
from plotUtils23 import PlotNotFoundError, InputDir, TFileContext
from plotUtils23  import GetPredictionsPlot, GetDataPlot, GetClosureStack
from plotUtils23 import graph_and_ratio, get_range_tga, clamp_expnd_r
from utils23 import lumi_dict
from variablesInfo import getVariablesInfo
import cmsstyle
import PersonalInfo
from Colours import Evidence, Warn

from array import array

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
                  help= "Regular expression; only plot names that match it will be used. Default: all")

parser.add_argument("-s", "--skip", dest="Skip",
                  default=None,
                  help= "Plots names that match this regex will be skipped")

parser.add_argument("-m", "--mcset", dest="mcSet", choices=['mad', 'pow'],
                  default="mad",
                  help= "Monte Carlo Set, pow for Powheg, mad for amcatnlo")

parser.add_argument("-p", "--prediction-type", dest="predType",
                  default="fullMC",
                  help= "Type of prediction. lepCR = non-prompt leptons from CRs, rare background from MC; fullMC = all from MC; fakeMC = use MC in CRs instead of data")

parser.add_argument("-l", "--lumiProj", dest="LumiProj",
                  default="",
                  help="Lumi projection")

parser.add_argument("-o", "--outputDir", dest="outputDir",
                  help="Directory where save plots. Default is based on inputDir")

parser.add_argument("-A", "--Analysis", dest="Analysis", choices=['VVXAnalyzer', 'VVGammaAnalyzer', 'ZZAnalyzer'],
                  default="VVXAnalyzer",
                  help="Analysis. Default is VVX; other options are ZZ, VBS and VVGamma")

parser.add_argument("-y", "--year", dest="year",
                  default="2016",
                  help= "valid inputs are 2016preVFP, 2016postVFP, 2017, 2018, Run2")

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

parser.add_argument('--skip-missing', action='store_true',
                    help='Don\'t crash if a plot is missing; instead continue with the others')

parser.add_argument('--allow-empty-data', action='store_true',
                    help='If the data histogram is empty, create and use an empty one instead')

parser.add_argument('--draw-label', action='store_true', dest='draw_label',
                    default=True,
                    help='Draw a textbox with the name of the region in the plot (default: %(default)s)')

parser.add_argument('--no-draw-label', action='store_false', dest='draw_label',
                    help='Set "%(dest)s" to false')

parser.add_argument('--region-label', default=None,
                    help='Specify manually the region label to be drawn. If it contains a "%%s", the region name will be subsituted')

parser.add_argument('--force-positive'   , action='store_true' , dest='forcePositive', help='Do `Scale(-1)` in regions with negative fake lepton transfer factor (default = %(default)s)')
parser.add_argument('--no-force-positive', action='store_false', dest='forcePositive')

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
    Double_t factor = TMath::Exp(fabs(ymax/ymin)*0.08);
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
Analysis   = options.Analysis
year       = options.year

inputDir   = InputDir(basedir=options.inputDir, year=year, region=region, analyzer=Analysis)
if(options.outputDir is None):
    split = options.inputDir.split('_')
    outputPrefix = os.path.join(PersonalInfo.personalFolder, '_'.join(split[1:]) if len(split) > 1 else 'last')
else:
    outputPrefix = os.path.join('' if options.outputDir.startswith('/') else PersonalInfo.personalFolder, options.outputDir)
OutputDir  = os.path.join(outputPrefix, Analysis, year, predType, region)

try:
    os.stat(OutputDir)
except OSError as e:
    if(not e.errno == 2): raise e  # 2 = No such file or directory
    os.makedirs(OutputDir)  # mkdir() = mkdir  ;  makedirs() = mkdir -p


cmsstyle.setCMSStyle()
ROOT.gROOT.SetBatch(True)

if LumiProj != "":
    inputDir.basedir += LumiProj+"fbm1_"
    lumi = LumiProj
else:
    lumi = lumi_dict[year]['value']
lumi = lumi/1000.
cmsstyle.SetExtraText('')
cmsstyle.SetEnergy(13, unit='TeV')
cmsstyle.SetLumi('{:.3g}'.format(lumi))

# Change the thickness of the MC stat error band
ROOT.gStyle.SetHatchesLineWidth(2)
ROOT.gStyle.SetHatchesSpacing(1)

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


missing_plots = []
for Var in variables:
    info = VarInfo[Var]
    info.update({'name':Var})

    DoData = optDoData and (info.get('unblind', True) or options.unblind or region[:2] != 'SR')
    
    # "Temporary" hack for closure test of photon fake rate
    if False: #'PhFRClosure' in Var and 'PASS' in Var:
        hMC, leg = GetClosureStack(region, inputDir.get_path(), info, forcePositive=options.forcePositive, verbosity=options.verbosity)
    else:
        if info.get('special'):
            info['name'] = info['stack']['plot']
        try:
            (hMC, leg) = GetPredictionsPlot(inputDir, info, predType, mcSet, forcePositive=options.forcePositive, verbosity=options.verbosity)
        except PlotNotFoundError as e:
            if(options.skip_missing):
                missing_plots.append(e)
                continue
            else:
                raise e

    if(not hMC.GetStack()):
        print Evidence('ERROR'), 'skipping', Var, 'because: no MC'
        continue

    if(DoData):
        if info.get('special'):
            info['name'] = info['data']['plot']
        try:
            histodata = GetDataPlot(inputDir, info, forcePositive=options.forcePositive, verbosity=options.verbosity)
        except PlotNotFoundError as e:
            if(options.skip_missing):
                missing_plots.append(e)
                continue
            elif(options.allow_empty_data):
                print(Warn('ERROR')+': missing data histogram for "%s"' %(info['name']))
                missing_plots.append(e)
                # Copy the MC histogram and set all bins to 0
                histodata = hMC.GetStack().First().Clone()
                for b in range(0, histodata.GetNbinsX()+2):
                    histodata.SetBinContent(b, 0)
                    histodata.SetBinError  (b, 0)
            else:
                raise e

        if(not histodata):
            print Evidence('ERROR'), 'skipping', Var, 'because: no data'
            continue

    hStackSum = hMC.GetStack().Last()
    if(not DoData):
        histodata = ROOT.TH1F(hMC.GetStack().Last())
        histodata.SetName("histodata")
        histodata.Reset()

    # Check for underflow or overflow
    overflow_fraction  = hStackSum.GetBinContent(hStackSum.GetNbinsX()+1) / hStackSum.Integral(0, -1)
    underflow_fraction = hStackSum.GetBinContent(0                      ) / hStackSum.Integral(0, -1)
    has_overflow  = overflow_fraction  > 0.1 # Overflow  is > 10% of total
    has_underflow = underflow_fraction > 0.1 # Underflow is > 10% of total
    if(has_overflow ):
        if(options.verbosity >= 1):
            print Warn('WARN'), 'overflow (%.1f %%)'  %(100*overflow_fraction )
    if(has_underflow):
        if(options.verbosity >= 1):
            print Warn('WARN'), 'underflow (%.1f %%)' %(100*underflow_fraction)

    # X range
    draw_overflow  = info.get('draw_overflow' , False)
    draw_underflow = info.get('draw_underflow', False)
    xaxis  = hStackSum.GetXaxis()
    bx_min = (0 if draw_underflow else 1)
    bx_max = (1 if draw_overflow  else 0) + xaxis.GetNbins()

    xmin_info = info.get('xmin')
    xmax_info = info.get('xmax')
    if(xmin_info is not None):
        bx_min = xaxis.FindFixBin(xmin_info + abs(xmin_info)*1e-6) # in case the requested xmin is a bin edge, get the right bin
        if('draw_underflow' in info and not draw_underflow):
            print(Warn('WARN') + ' xmin overrides draw_underflow')
        draw_underflow = (bx_min == 0)
    if(xmax_info is not None):
        bx_max = xaxis.FindFixBin(xmax_info - abs(xmax_info)*1e-6) # in case the requested xmax is a bin edge, get the left bin
        if('draw_overflow'  in info and not draw_overflow ):
            print(Warn('WARN') + ' xmax overrides draw_overflow')
        draw_overflow  = (bx_max == xaxis.GetNbins()+1)

    x_min  = xaxis.GetBinLowEdge(bx_min)
    x_max  = xaxis.GetBinLowEdge(bx_max+1)

    # Fill an array of bin edges. This is needed to handle custom binnings
    xedges = array('d', (xaxis.GetNbins())*[0.])
    xaxis.GetLowEdge(xedges)
    xedges.append(xaxis.GetBinLowEdge(xaxis.GetNbins()+1)) # add right edge of last bin
    if(draw_underflow): xedges.insert(0, xaxis.GetBinLowEdge(0))
    if(draw_overflow ): xedges.append(   xaxis.GetBinLowEdge(xaxis.GetNbins()+2))
    # Filter the bin edges it so that only those within the requested limits remain
    xedges = array('d', [e for e in xedges if e >= x_min and e <= x_max])# xedges[bx_min-1:bx_max+1]

    # TGraphs to draw in the upper plot (data) and in the ratio plot
    graphData, tgaData = graph_and_ratio(histodata, hStackSum, xedges=xedges, bx_min=bx_min, bx_max=bx_max,
                                         unblind=DoData, remove_zeros=True)
    if(info.get('special')):
        # Avoid the large error bars that are associated with poisson errors
        graphData = ROOT.TGraphAsymmErrors(histodata)

    # Y range - upper plot
    y_max = info.get('ymax', False)
    if(not y_max):
        y_max_data = ROOT.TMath.MaxElement(graphData.GetN(), graphData.GetEYhigh()) + ROOT.TMath.MaxElement(graphData.GetN(), graphData.GetY()) if DoData else 0.
        y_max_MC = hStackSum.GetBinContent(hStackSum.GetMaximumBin()) + hStackSum.GetBinError(hStackSum.GetMaximumBin())
        y_max = max(y_max_MC, y_max_data)
        y_max *= info.get('scale_ymax', 1.37)

        if info.get('logy', False):
            y_max *= 10
    y_min = info.get('ymin', 0 if not info.get('logy') else hMC.GetMinimum())

    # Ratio range
    if(DoData and tgaData.GetN() > 0):
        y_min_r, y_max_r = get_range_tga(tgaData, include_err=True)
    else:
        y_max_r = 1.
        y_min_r = 1.
    deltaY = (y_max_r - y_min_r)
    y_min_r, y_max_r = clamp_expnd_r(y_min_r, y_max_r,
                                     min_lo=0. , max_lo=0.5,
                                     min_hi=1.5, max_hi=15.)
    y_max_r = info.get('ratio_ymax', y_max_r)
    y_min_r = info.get('ratio_ymin', y_min_r)

    # Make the canvas
    canvas = cmsstyle.cmsDiCanvas(
        'canvas',
        x_min=x_min,
        x_max=x_max,
        y_min=y_min,
        y_max=y_max,
        r_min=y_min_r,
        r_max=y_max_r,
        nameXaxis=info.get('title', ''),
        nameYaxis='Events',
        nameRatio='Data/Pred.',
        square=True,
        iPos=0,
        extraSpace=0.02,
    )
    # as of cmsstyle 0.4.2, this resets the TStyle, so any change to e.g. the axis title offset must be done after this call
    pad1 = canvas.GetPad(1)
    pad2 = canvas.GetPad(2)
    hFrameUp = pad1.FindObject("hframe")
    hFrameDn = pad2.FindObject("hframe")

    # Log scale
    if info.get('logy', False):
        pad1.SetLogy()
        if(y_min > 0 and math.log10(y_max/y_min) < 4):
            hFrameUp.GetYaxis().SetMoreLogLabels()
        hFrameUp.GetYaxis().SetLabelOffset(0.010)
        if(y_max < 10000):
            hFrameUp.GetYaxis().SetNoExponent()
            cmsstyle.UpdatePad(pad1)

    if info.get('logx', False):
        pad1.SetLogx()
        pad2.SetLogx()

    # Draw the THStack
    pad1.cd()
    hMC.Draw("hist same")

    # The THStack axis cannot be modified before it has been drawn
    hMC.GetXaxis().SetRange(bx_min, bx_max)
    if(DoData):
        histodata.GetXaxis().SetRange(bx_min, bx_max)

    # Error band in the upper canvas
    hMCErr = deepcopy(hStackSum)

    hMCErr.SetFillStyle(3345)
    hMCErr.SetMarkerStyle(1)
    hMCErr.SetFillColor(ROOT.kBlack)
    hMCErr.Draw("sameE2")
    leg.AddEntry(hMCErr, "Stat. only", "f")
    
    if DoData:
        if(info.get('text')):
            texec = ROOT.TExec("texec", 'drawtext("{}");'.format(graphData.GetName()))
            graphData.GetListOfFunctions().Add(texec)
            graphData.Draw("samep text")
        else:
            graphData.Draw("samep")
        leg.AddEntry(graphData, info.get('data', dict()).get('legend', 'Data'), "lpe")

    x1 = leg.GetX1()
    x2 = leg.GetX2()
    shift = 0.78 - (x1+x2)/2
    leg.SetX1(x1+shift)
    leg.SetX2(x2+shift)
    leg.Draw("same")

    if(options.draw_label):
        region_text = ROOT.TText()
        region_text.SetNDC()
        if options.region_label is not None:
            try:
                text = options.region_label %(region)
            except TypeError:
                text = options.region_label
        else:
            text = region
        region_text.SetText(pad1.GetLeftMargin()+0.05, 1-pad1.GetTopMargin()-0.1, text)
        region_text.SetTextSize(.05)
        region_text.Draw('same')

    # Ratio plot
    pad2.cd()

    Line = ROOT.TLine(x_min, 1, x_max, 1)
    Line.SetLineWidth(2)
    Line.SetLineStyle(7)
    
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

    # Fixes to X axis (ratio pad)
    model_axis = histodata.GetXaxis()
    draw_axis  = hFrameDn.GetXaxis()
    draw_axis.SetLabelSize(0.1  ) # cmsstyle defaults to 0.1171875
    draw_axis.SetTitleSize(0.12 ) # cmsstyle defaults to 0.140625
    draw_axis.SetTitleOffset(1.05)# cmsstyle defaults to 0.9

    # No exponent on x axis if the range is small enough
    if (draw_axis.GetXmin() > 0.001 and draw_axis.GetXmax() < 1000):
        draw_axis.SetNoExponent()

    # Deal with alphanumeric labels by resetting and redrawing the xaxis
    if(model_axis.IsAlphanumeric()):
        draw_axis.SetAlphanumeric()
        draw_axis.Set(model_axis.GetNbins(), model_axis.GetXmin(), model_axis.GetXmax())
        for b in range(1, model_axis.GetNbins()+1):
            draw_axis.SetBinLabel(b, model_axis.GetBinLabel(b))
        cmsstyle.GetcmsCanvasHist(pad2).Draw()
        canvas.RedrawAxis()
        cmsstyle.UpdatePad(canvas)
        pad2.cd()

    Line.Draw()
    tgaData.Draw("PE0 same")

    if(not DoData):
        xm, xM = hMC.GetXaxis().GetXmin(), hMC.GetXaxis().GetXmax()
        xstart = (xm + xM)/2 - (xM - xm)/8
        text = ROOT.TText(xstart, 1.1, "BLINDED")
        text.SetTextSize(.12)
        text.Draw("same")
    
    Title=Var+"_"+mcSet #+"_"+region
    if(not DoData): Title+='_blind'
    
    canvas.SetTitle(Title)

    for ext in ('png', 'pdf'): #, 'root', 'eps'
        canvas.SaveAs(os.path.join(OutputDir, Title+'.'+ext))

    with TFileContext(os.path.join(OutputDir, Title+'.root'), 'recreate') as f:
        f.cd()
        leg.Write("legend")
        hMC.Write("MC")
        graphData.Write("data")

    del histodata, tgaData

for missing_plot in missing_plots:
    print(missing_plot)
