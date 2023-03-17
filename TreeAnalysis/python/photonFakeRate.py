#!/usr/bin/env python2

##############################################
#                                         
# Creates files for photon fake rate
# 
##############################################

from __future__ import print_function
from os import path, environ
import copy
import ROOT
# from readSampleInfo import *
# from Colours import *
from math import log10
from ctypes import c_int
from array import array

from plotUtils import TFileContext, makedirs_ok

ROOT.gStyle.SetPaintTextFormat(".2f")

def getPlots(inputdir, sample, plots):
    fname = path.join(inputdir, sample+".root")
    retrieved = []
    with TFileContext(fname) as rFile:
        for plot in plots:
            h = rFile.Get(plot)
            if(not h):
                print('Warning: Could not get "%s" from "%s"' % (plot, fname))
                retrieved.append(None)
            else:
                retrieved.append( copy.deepcopy(h) )
    return retrieved


def iterate_bins(h, **kwargs):
    loX = 0 if kwargs.get("underX") else 1
    loY = 0 if kwargs.get("underY") else 1
    loZ = 0 if kwargs.get("underZ") else 1
    upX = 0 if kwargs.get("overfX") else 1
    upY = 0 if kwargs.get("overfY") else 1
    upZ = 0 if kwargs.get("overfZ") else 1
    for i in range(loX, h.GetNbinsX() + 1 + upX):
        for j in range(loY, h.GetNbinsY() + 1 + upY):
            for k in range(loZ, h.GetNbinsZ() + 1 + upZ):
                yield h.GetBin(i, j, k)


def fix_neg_bins(h):
    for b in iterate_bins(h, underX=True, underY=True, underZ=False, overX=True, overY=True, overZ=False):
        if(h.GetBinContent(b) < 0):
            h.SetBinContent(b, 0.)


# def _show_overflow_axis(ax, u_label=None, o_label=None):
#     ax.SetRange(0, ax.GetNbins() + 1)
#     if(u_label is None): u_label = '0'
#     if(o_label is None): o_label = '\infty'
#     ax.ChangeLabel(0                          , -1., -1., -1, -1, -1, u_label)
#     ax.ChangeLabel(ax.GetNdivisions()%100 - 1 , -1., -1., -1, -1, -1, o_label)

# def show_overflow(h, **kwargs):
#     doX = kwargs.get("doX") or any([ k in ['ux_label', "ox_label"] for k in kwargs.keys() ])
#     doY = kwargs.get("doY") or any([ k in ['uy_label', "oy_label"] for k in kwargs.keys() ])
#     doZ = kwargs.get("doZ") or any([ k in ['uz_label', "oz_label"] for k in kwargs.keys() ])
#     if(doX):
#         _show_overflow_axis( h.GetXaxis(), kwargs.get("ux_label"), kwargs.get("ox_label") )
#     if(doY and h.Class().InheritsFrom("TH2")):
#         _show_overflow_axis( h.GetYaxis(), kwargs.get("uy_label"), kwargs.get("oy_label") )
#     if(doZ and h.Class().InheritsFrom("TH3")):
#         _show_overflow_axis( h.GetZaxis(), kwargs.get("uz_label"), kwargs.get("oz_label") )


def linesAndTicks(h, logx=False, logy=False):
    lines = []
    ticks = []
    xmin, xmax = h.GetXaxis().GetXmin(), h.GetXaxis().GetXmax()
    ymin, ymax = h.GetYaxis().GetXmin(), h.GetYaxis().GetXmax()
    xsize = xmin + 0.03*(xmax - xmin)  # Ticks are 3 % of axis length
    if(logx): xsize = xmin + 10**(0.03*log10(xmax/xmin))
    ysize = ymin + 0.03*(ymax - ymin)
    if(logy): ysize = ymin + 10**(0.03*log10(ymax/ymin))
        
    for i in range(1, h.GetNbinsX()):
        edge = h.GetXaxis().GetBinUpEdge(i)
        
        lines.append( ROOT.TLine(edge, ymin, edge, ymax) )
        lines[-1].SetLineStyle(2); lines[-1].SetLineColor(ROOT.kGray +2); lines[-1].SetLineWidth(1)

        ticks.append( ROOT.TLine(edge, ymin, edge, ysize) )
        ticks[-1].SetLineStyle(1); ticks[-1].SetLineColor(ROOT.kBlack); ticks[-1].SetLineWidth(1)
    
    for i in range(1, h.GetNbinsY()):
        edge = h.GetYaxis().GetBinUpEdge(i)
        
        lines.append( ROOT.TLine(xmin, edge, xmax, edge) )
        lines[-1].SetLineStyle(2); lines[-1].SetLineColor(ROOT.kGray +2); lines[-1].SetLineWidth(1)

        ticks.append( ROOT.TLine(xmin, edge, xsize, edge) )
        ticks[-1].SetLineStyle(1); ticks[-1].SetLineColor(ROOT.kBlack); ticks[-1].SetLineWidth(1)
    
    return lines, ticks

def beautify(canvas, hist, logx=False, logy=False, logz=False):
    canvas.cd()
    # hist.GetXaxis().SetNdivisions(0)
    hist.GetXaxis().SetRange(1, hist.GetXaxis().GetNbins())
    hist.Draw("colz texte")

    if(logx):
        canvas.SetLogx()
        hist.GetXaxis().SetNoExponent()
        hist.GetXaxis().SetMoreLogLabels()
    if(logy):
        canvas.SetLogy()
        hist.GetYaxis().SetNoExponent()
        hist.GetYaxis().SetMoreLogLabels()
    if(logz):
        canvas.SetLogz()
        hist.GetZaxis().SetNoExponent()
        hist.GetZaxis().SetMoreLogLabels()

    lines, ticks = linesAndTicks(hist, logx, logy)
    to_preserve = lines + ticks
    for l in to_preserve:
        l.Draw("same")
    return to_preserve


def get_path_results():
    if(environ.get("CMSSW_BASE")): path_base = path.join(environ['CMSSW_BASE'], "src/VVXAnalysis/TreeAnalysis/results")
    else:                          path_base = "results" # "~/Work/Analysys/rsync_results/Updated/EXT"
    assert path.exists(path_base), 'Please call this script from "VVXAnalysis/TreeAnalysis" or do `cmsenv`'
    return path_base
    


def fakeRateABCD(sample_data, sample_prompt, analyzer="VVGammaAnalyzer", year=2016, region="CRLFR", logx=False, logy=False, fixNegBins=False):
    print("Fake Rate ABCD: data     =", sample_data  )
    print("Fake Rate ABCD: promptMC =", sample_prompt)
    path_in = "{:s}/{}/{:s}_{:s}".format(get_path_results(), year, analyzer, region)
    outfname = "data/ABCD_FR_{:s}-{:s}.root".format(sample_data["name"], sample_prompt["name"])
    print('Output (with prompt MC subtraction) in: "{:s}"'.format(outfname))
    
    ##### First: only data ####
    
    hdata_A, hdata_B, hdata_C, hdata_D = getPlots(path_in, sample_data  ["file"], [ "PhFR_%s" % (c) for c in 'ABCD' ] )
    assert (hdata_A and hdata_B and hdata_C and hdata_D), "Could't get all 4 data histograms!"

    if(sample_data.get("fixNegBins", False)):
        for h in [hdata_A, hdata_B, hdata_C, hdata_D]:
            fix_neg_bins(h)
    
    ##### Second: data - promptMC #####
    if(sample_prompt is None):
        return

    hprompt_A, hprompt_B, hprompt_C, hprompt_D = getPlots(path_in, sample_prompt["file"], [ "PhFR_%s_prompt" % (c) for c in 'ABCD' ] )
    assert (hprompt_A and hprompt_B and hprompt_C and hprompt_D), "Could't get all 4 promptMC histograms!"

    if(sample_prompt.get("fixNegBins", False)):
        for h in (hprompt_A, hprompt_B, hprompt_C, hprompt_D):
            fix_neg_bins(h)

    hdata_A.Add(hprompt_A, -1)
    hdata_B.Add(hprompt_B, -1)
    hdata_C.Add(hprompt_C, -1)
    hdata_D.Add(hprompt_D, -1)

    ## Fake rate = C/D ##
    cFR = ROOT.TCanvas( "cFR", "Fake Rate {:s} - {:s}".format(sample_data["name"], sample_prompt["name"]), 200, 10, 1200, 900 )
    cFR.cd()
    
    hFR = hdata_C
    hFR.Divide(hdata_D)
    hFR.SetTitle( "Photon Fake Rate: C/D (from {} - {})".format(sample_data["title"], sample_prompt["title"]) )
    hFR.SetName( "PhFR" )
    
    if(fixNegBins):
        fix_neg_bins(hFR)
    
    with TFileContext(outfname, "UPDATE") as tf:
        tf.cd()
        hFR.Write(hFR.GetName(), ROOT.TObject.kOverwrite)
        
    to_preserve_FR = beautify(cFR, hFR, logx, logy)
    
    for ext in ["png"]:
        cFR.SaveAs( "Plot/PhFR/ABCD_FR_{:s}-{:s}.{:s}".format(sample_data["name"], sample_prompt["name"], ext) )
    del cFR, to_preserve_FR
    
    ## Estimate = B*C/D ##
    cES = ROOT.TCanvas( "cES", "Estimate {} - {}".format(sample_data["name"], sample_prompt["name"]) , 200, 10, 1200, 900 )
    cES.cd()

    hES = copy.deepcopy(hFR)
    hES.Multiply(hdata_B)
    hES.SetTitle( "Fake Photon estimate: B*C/D (from {:s} - {:s})".format(sample_data["title"], sample_prompt["title"]) )
    hES.SetName( "FakeEstimate" )

    if(fixNegBins):
        for b in iterate_bins(hES):
            if(hES.GetBinContent(b) < 0): hES.SetBinContent(b, 0.)

    with TFileContext(outfname, "UPDATE") as tf:
        tf.cd()
        hES.Write(hES.GetName(), ROOT.TObject.kOverwrite)
    
    to_preserve_ESdata = beautify(cES, hES, logx, logy, logz=True)
    
    for ext in ["png"]:
        cES.SaveAs( "Plot/PhFR/ABCD_Estimate_{:s}-{:s}.{:s}".format(sample_data["name"], sample_prompt["name"], ext) )
    del cES, to_preserve_ESdata

    ## k-Factor = (C/D)/(A/B)
    cKF = ROOT.TCanvas('cKF', 'kFactor {} - {}'.format(sample_data['name'], sample_prompt['name']), 1200,900)
    cKF.cd()
    
    hKF = copy.deepcopy(hES)
    hKF.Divide(hdata_A)
    hKF.SetTitle('k-Factor: (C/D)/(A/B) (from {:s} - {:s})'.format(sample_data['title'], sample_prompt['title']))
    hKF.SetName('kFactor')

    with TFileContext(outfname, 'UPDATE') as tf:
        tf.cd()
        hKF.Write(hKF.GetName(), ROOT.TObject.kOverwrite)
    
    to_preserve_KF = beautify(cKF, hKF, logx, logy)
    cKF.SaveAs( 'Plot/PhFR/ABCD_kFactor_{:s}-{:s}.{:s}'.format(sample_data['name'], sample_prompt['name'], 'png') )
    del cKF, to_preserve_KF
    
    return hFR, hES, hKF


    return hFR, hES


def fakeRateABCD_noSubtract(sample, analyzer="VVGammaAnalyzer", year=2016, region="CRLFR", logx=False, logy=False, fixNegBins=False):
    print("Fake Rate ABCD: sample   =", sample  )
    path_in = "{:s}/{}/{:s}_{:s}".format(get_path_results(), year, analyzer, region)
    outfname = "data/ABCD_FR_{:s}.root".format(sample["name"])
    print('Output (without MC subtraction) in: "{:s}"'.format(outfname))
    
    hA, hB, hC, hD = getPlots(path_in, sample["file"], [ "PhFR_%s_nonprompt" % (c) for c in 'ABCD' ] )
    assert (hA and hB and hC and hD), "Could't get all 4 data histograms!"

    if(sample.get("fixNegBins", False) and fixNegBins):
        for h in [hA, hB, hC, hD]:
            fix_neg_bins(h)

    ## Fake rate = C/D ##
    cFR = ROOT.TCanvas( "cFR_ABCD_nosub", "Fake Rate "+sample["name"], 200, 10, 1200, 900 )
    cFR.cd()
    
    hFR = copy.deepcopy(hC)
    hFR.Divide(hD)
    hFR.SetTitle( "Photon Fake Rate: C/D (from {})".format(sample["title"]) )
    hFR.SetName( "PhFR" )

    with TFileContext(outfname, "UPDATE") as tf:
        tf.cd()
        hFR.Write(hFR.GetName(), ROOT.TObject.kOverwrite)

    to_preserve_FR = beautify(cFR, hFR, logx, logy)
    # Probably python/ROOT gc deletes these graphic objects and they disappear from the canvas. Holding a reference to them seems to fix it
        
    for ext in ["png"]:
        cFR.SaveAs( "Plot/PhFR/ABCD_FR_{:s}.{:s}".format(sample["name"], ext) )
    del cFR, to_preserve_FR

    ## Estimate = B*C/D ##
    cES = ROOT.TCanvas( "cES_nosub", "Estimate "+sample["name"], 200, 10, 1200, 900 )
    cES.cd()

    hES = copy.deepcopy(hFR)  # rename for clarity
    hES.Multiply(hB)
    hES.SetTitle( "Fake Photon estimate: B*C/D (from {:s})".format(sample["title"]) )
    hES.SetName ( "FakeEstimate" )

    with TFileContext(outfname, "UPDATE") as tf:
        tf.cd()
        hES.Write(hES.GetName(), ROOT.TObject.kOverwrite)

    to_preserve_ES = beautify(cES, hES, logx, logy, logz=True)
    for ext in ["png"]:
        cES.SaveAs( "Plot/PhFR/ABCD_Estimate_{:s}.{:s}".format(sample["name"], ext) )
    del cES, to_preserve_ES

    ## k-Factor = (C/D)/(A/B)
    cKF = ROOT.TCanvas('cKF_nosub', 'kFactor '+sample['name'], 1200,900)
    cKF.cd()

    hKF = copy.deepcopy(hES)
    hKF.Divide(hA)
    hKF.SetTitle('k-Factor: (C/D)/(A/B) from {:s}'.format(sample['title']))
    hKF.SetName('kFactor')

    with TFileContext(outfname, 'UPDATE') as tf:
        tf.cd()
        hKF.Write(hKF.GetName(), ROOT.TObject.kOverwrite)
    
    to_preserve_KF = beautify(cKF, hKF, logx, logy)
    cKF.SaveAs( 'Plot/PhFR/ABCD_kFactor_{:s}.{:s}'.format(sample['name'], 'png') )
    del cKF, to_preserve_KF
    
    return hFR, hES, hKF


def fakeRateLtoT_noSubtract(sample, plotname='PhFR_LtoT_%s_nonprompt', analyzer='VVGammaAnalyzer', year=2016, region='CRLFR', logx=False, logy=False, fixNegBins=False):
    print('Fake Rate LtoT: sample   =', sample)
    path_in = '{:s}/{}/{:s}_{:s}'.format(get_path_results(), year, analyzer, region)
    outfname = 'data/LtoT_FR_{:s}.root'.format(sample['name'])
    print('Output in: "{:s}"'.format(outfname))

    h5, h4, h3 = getPlots(path_in, sample['file'], [ plotname % (s) for s in [5, 4, 3] ] )
    assert (h5 and h4 and h3), "Could't get the 3 histograms!"
    
    hPASS = h5
    hFAIL = h3
    hFAIL.Add(h4)
    hFAIL.Add(hPASS)

    if(sample.get('fixNegBins', False) and fixNegBins):
        for h in [hPASS, hFAIL]:
            fix_neg_bins(h)

    ## Fake rate = PASS/FAIL ##
    cFR = ROOT.TCanvas( 'cFR_LtoT_nosub', 'Fake Rate '+sample['name'], 1200, 900 )
    cFR.cd()
    
    hFR = copy.deepcopy(hPASS)
    hFR.Divide(hFAIL)
    hFR.SetTitle( 'Photon Fake Rate: 5/(3+4+5) (from {:s})'.format(sample['title']) )
    hFR.SetName( 'PhFR_LtoT' )

    with TFileContext(outfname, 'UPDATE') as tf:
        tf.cd()
        hFR.Write(hFR.GetName(), ROOT.TObject.kOverwrite)

    to_preserve = beautify(cFR, hFR, logx, logy)
    
    for ext in ['png']:
        cFR.SaveAs( 'Plot/PhFR/LtoT_FR_{:s}.{:s}'.format(sample['name'], ext) )
    del cFR, to_preserve

    return hFR

    
def fakeRateLtoT(sample_data, sample_prompt, analyzer='VVGammaAnalyzer', year=2016, region='CRLFR', logx=False, logy=False, fixNegBins=False):
    print('Fake Rate LtoT: data     =', sample_data  )
    print('Fake Rate LtoT: promptMC =', sample_prompt)
    path_in = '{:s}/{}/{:s}_{:s}'.format(get_path_results(), year, analyzer, region)
    outfname_full = 'data/LtoT_FR_{:s}-{:s}.root'.format(sample_data['name'], sample_prompt['name'])
    print('Output (with prompt MC subtraction) in: "{:s}"'.format(outfname_full))
    
    ##### First: take data ####
    hdata_5  , hdata_4  , hdata_3   = getPlots(path_in, sample_data['file']  , [ 'PhFR_LToT_%d'        % (s) for s in [5, 4, 3] ] )
    hprompt_5, hprompt_4, hprompt_3 = getPlots(path_in, sample_prompt['file'], [ 'PhFR_LToT_%d_prompt' % (s) for s in [5, 4, 3] ] )
    assert (hdata_5   and hdata_4   and hdata_3  ), "Could't get the 3 data histograms!"
    assert (hprompt_5 and hprompt_4 and hprompt_3), "Could't get the 3 prompt MC histograms"
    
    hdata_PASS   = hdata_5
    hdata_FAIL   = hdata_3
    hdata_FAIL  .Add(hdata_4  )
    hprompt_FAIL = hprompt_3
    hprompt_PASS = hprompt_5
    hprompt_FAIL.Add(hprompt_4)
    
    hdata_FAIL  .Add(hdata_PASS  )
    hprompt_FAIL.Add(hprompt_PASS)

    if(sample_data.get('fixNegBins', False) and fixNegBins):
        for h in [hdata_PASS, hdata_FAIL]:
            fix_neg_bins(h)
    if(sample_prompt.get('fixNegBins', False) and fixNegBins):
        for h in (hprompt_PASS, hprompt_FAIL):
            fix_neg_bins(h)

    hdata_PASS.Add(hprompt_PASS, -1)
    hdata_FAIL.Add(hprompt_FAIL, -1)

    ## Fake rate = PASS/FAIL ##
    cFR = ROOT.TCanvas( 'cFR_LtoT', 'Fake Rate {:s} - {:s}'.format(sample_data['name'], sample_prompt['name']), 1200, 900 )
    cFR.cd()
    
    hFR = hdata_PASS
    hFR.Divide(hdata_FAIL)
    hFR.SetTitle( 'Photon Fake Rate: 5/(3+4+5) (from {} - {})'.format(sample_data['title'], sample_prompt['title']) )
    hFR.SetName( 'PhFR' )
    
    if(fixNegBins):
        fix_neg_bins(hFR)
    
    with TFileContext(outfname_full, 'UPDATE') as tf:
        tf.cd()
        hFR.Write(hFR.GetName(), ROOT.TObject.kOverwrite)
        
    to_preserve_FR = beautify(cFR, hFR, logx, logy)
        
    for ext in ['png']: 
        cFR.SaveAs( 'Plot/PhFR/LtoT_FR_{:s}-{:s}.{:s}'.format(sample_data['name'], sample_prompt['name'], ext) )
    del cFR, to_preserve_FR

    return hFR


def plotRatio(h1, h2, picture_name="ratio", title="ratio"):
    assert h1, "h1 missing"
    assert h2, "h2 missing"

    ratio = copy.deepcopy(h1)
    ratio.Divide(h2)
    ratio.SetTitle(title)
    
    c = ROOT.TCanvas("cratio_{:s}".format(picture_name), "ratio", 1200, 900)
    c.cd()
    
    ratio.Draw("colz texte")
    to_preserve = beautify(c, ratio, True, False)
    # ratio.GetYaxis().SetRange(1, ratio.GetXaxis().GetNbins())
    c.SaveAs('Plot/PhFR/{:s}.png'.format(picture_name))
    del to_preserve, c

    return ratio


def plotProfiled(h2, picture_name=None, title='profile', direction='X'):
    if(not h2.Class().InheritsFrom("TH2")):
        print('ERROR: "{:s}" does not inherit from TH2'.format(h2.GetName()))
        return
    if(picture_name is None):
        picture_name = h2.GetName()

    _style = [ {"color":ROOT.kRed, "marker":ROOT.kFullCircle}, {"color":ROOT.kBlue, "marker":ROOT.kFullTriangleUp}, {"color":ROOT.kGreen, "marker":ROOT.kFullTriangleDown}, {"color":ROOT.kBlack, "marker":ROOT.kFullSquare}, {"color":ROOT.kMagenta, "marker":ROOT.kFullCrossX}, {"color":ROOT.kCyan, "marker":ROOT.kFullStar}]
    
    cprof = ROOT.TCanvas("c_"+h2.GetName(), h2.GetName(), 1200, 900)
    cprof.cd()
    if  (direction=='Y'):
        axis_2D_proj = h2.GetYaxis()
        axis_2D_draw = h2.GetXaxis()
        nbins_proj = h2.GetNbinsY()
        nbins_draw = h2.GetNbinsX()
    elif(direction=='X'):
        axis_2D_proj = h2.GetXaxis()
        axis_2D_draw = h2.GetYaxis()
        nbins_proj = h2.GetNbinsX()
        nbins_draw = h2.GetNbinsY()
    
    h1s = []
    draw_edges = array('d')
    for i in range(1, axis_2D_draw.GetNbins()+2):
        draw_edges.append(axis_2D_draw.GetBinLowEdge(i))
    draw_edges.append(axis_2D_draw.GetBinUpEdge(axis_2D_draw.GetNbins()+1))
    # print('edges:',draw_edges)

    zmax, zmin = 0, 0
    legend = ROOT.TLegend(.7,.7,.9,.9)
    
    for j in range(1, axis_2D_proj.GetNbins()+1):
        h1 = ROOT.TH1D("{:s}_bin{:d}".format(h2.GetName(), j), title, nbins_draw+1, draw_edges)#, axis_2D_draw.GetNbins(), axis_2D_draw.GetXbins().GetArray())
        h1.GetXaxis().SetTitle(axis_2D_draw.GetTitle())
        for i in range(0, h1.GetNbinsX()+2):
            bin_2d = h2.GetBin(i, j)
            content = h2.GetBinContent(bin_2d)
            error = h2.GetBinError  (bin_2d)
            zmax = max(zmax, content+error)
            zmin = min(zmin, content-error)
            # print("\t>>> bin:", i, "bin_2d:", bin_2d, "content: {:.2f} +- {:.2f}".format(content, error))
            h1.SetBinContent(i, content)
            h1.SetBinError  (i, error)
        h1s.append(h1)
        legend.AddEntry( h1, "{}<{:s}<{}".format(axis_2D_proj.GetBinLowEdge(j), axis_2D_proj.GetTitle().split(' ')[0], axis_2D_proj.GetBinUpEdge(j)) )

    for i,h in enumerate(h1s):
        h.SetLineColor(  _style[i]["color"] )
        h.SetMarkerColor(_style[i]["color"] )
        h.SetMarkerStyle(_style[i]["marker"])
    
    h1s[0].GetYaxis().SetRangeUser(zmin - (zmax-zmin)/10, zmax + (zmax-zmin)/10)
    h1s[0].Draw()
    for h1 in h1s[1:]:
        h1.Draw("P0 E1 same")
    
    legend.Draw("SAME")
    
    cprof.SaveAs( "Plot/PhFR/{:s}.png".format(picture_name) )


if __name__ == "__main__":
    sampleList = {
        "data"     : {"file": 'data'},
        "ZGToLLG"  : {"file": 'ZGToLLG', "title": "Z#gamma_{MC}", "fixNegBins": True},
        "Drell-Yan": {"file": 'DYJetsToLL_M50', "fixNegBins": True}
    }
    for key, sample in sampleList.items():
        sample.setdefault("name" , key)
        sample.setdefault("title", sample["name"])

    from argparse import ArgumentParser
    parser = ArgumentParser()
    parser.add_argument("-y", "--year"   , type=int, default=2016, choices=[2016, 2017, 2018])
    parser.add_argument("-m", "--method" , nargs="+", choices=["LtoT", "ABCD"], default="LtoT")
    parser.add_argument(      "--no-mc"  , dest="do_mc"  , action="store_false", help="Skip MC plots"  )
    parser.add_argument(      "--no-data", dest="do_data", action="store_false", help="Skip data plots")
    args = parser.parse_args()

    ROOT.gStyle.SetOptStat(0)
    ROOT.gROOT.SetBatch(True)

    makedirs_ok('data')
    makedirs_ok('Plot/PhFR')

    if("ABCD" in args.method):
        if(args.do_data):
            print('>>>ABCD data')
            hFR_ABCD_data, hES_data, hKF_data= fakeRateABCD(sampleList["data"], sampleList["ZGToLLG"], year=args.year, logx=True, fixNegBins=False)
            print()
        if(args.do_mc):
            print('>>>ABCD MC')
            hFR_ABCD_MC  , hES_MC  , hKF_MC  = fakeRateABCD_noSubtract(sampleList["Drell-Yan"]       , year=args.year, logx=True, fixNegBins=True)
            print()
        if(args.do_data and args.do_mc):
            print('>>>ABCD ratio')
            plotRatio(hFR_ABCD_data, hFR_ABCD_MC, picture_name="ABCD_ratio_data-ZG_over_DY_{}".format(args.year), title="Ratio FR(data-Z#gamma)/FR(DY)")
            print()

    if("LtoT" in args.method):
        if(args.do_data):
            print('>>>LtoT data')
            hFR_LtoT_data = fakeRateLtoT(sampleList["data"], sampleList["ZGToLLG"], year=args.year, logx=True, fixNegBins=False)
            plotProfiled(hFR_LtoT_data, picture_name='LtoT_FR_profiledX_data-ZGToLLG', title='FR(#gamma) vs #eta' , direction='X')
            plotProfiled(hFR_LtoT_data, picture_name='LtoT_FR_profiledY_data-ZGToLLG', title='FR(#gamma) vs p_{T}', direction='Y')
            print()
        if(args.do_mc):
            print('>>>LtoT MC')
            hFR_LtoT_MC   = fakeRateLtoT_noSubtract(sampleList["Drell-Yan"]       , year=args.year, logx=True, fixNegBins=True)
            print()
        if(args.do_data and args.do_mc):
            print('>>>LtoT ratio')
            ratioLtoT = plotRatio(hFR_LtoT_data, hFR_LtoT_MC, picture_name="LtoT_ratio_data-ZG_over_DY_{}".format(args.year), title="Ratio FR(data-Z#gamma)/FR(DY)")
            print()
