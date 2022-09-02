#!/usr/bin/env python

##############################################
#                                         
# Creates files for photon fake rate
# 
##############################################

from __future__ import print_function
import copy, sys, os
import ROOT
if '-b' in sys.argv:
    ROOT.gROOT.SetBatch(True)
# from readSampleInfo import *
# from Colours import *
# import ctypes
from math import log10

from plotUtils import TFileContext, makedirs_ok

ROOT.gStyle.SetPaintTextFormat(".2f")

def getPlots(inputdir, sample, plots):
    fname = os.path.join(inputdir, sample+".root")
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
    hist.GetXaxis().SetRange(1, hist.GetXaxis().GetNbins()+1)
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


def get_path_base():
    if(os.environ.get("CMSSW_BASE")): path_base = os.path.join(os.environ['CMSSW_BASE'], "src/VVXAnalysis/TreeAnalysis/results")
    else:                             path_base = "results" # "rsync_results/Updated/EXT"
    assert os.path.exists(path_base), 'Please call this script from "VVXAnalysis/TreeAnalysis" or do `cmsenv`'
    return path_base
    


def fakeRate(sample_data, sample_prompt, analyzer="VVXAnalyzer", region="CRLFR", logx=False, logy=False, fixNegBins=False):
    print("Fake Rate: data     =", sample_data  )
    print("Fake Rate: promptMC =", sample_prompt)
    makedirs_ok("data")
    makedirs_ok("Plots")
    outfname_data = "data/PhFR_from_{:s}.root".format(sample_data["name"])
    print('Output (without MC subtraction) in: "{:s}"'.format(outfname_data))
    
    path_in = "{:s}/2016/{:s}_{:s}".format(get_path_base(), analyzer, region)
    
    ##### First: only data ####
    
    hdata_A, hdata_B, hdata_C, hdata_D = getPlots(path_in, sample_data  ["file"], [ "PhFR_%s" % (c) for c in 'ABCD' ] )
    assert (hdata_A and hdata_B and hdata_C and hdata_D), "Could't get all 4 data histograms!"

    if(sample_data.get("fixNegBins", False)):
        for h in [hdata_A, hdata_B, hdata_C, hdata_D]:
            fix_neg_bins(h)

    ## Fake rate = C/D ##
    cFR_data = ROOT.TCanvas( "cFR_data", "Fake Rate "+sample_data["name"], 200, 10, 1200, 900 )
    cFR_data.cd()
    
    hFR_data = copy.deepcopy(hdata_C)
    hFR_data.Divide(hdata_D)
    hFR_data.SetTitle( "Photon Fake Rate: C/D (from {})".format(sample_data["title"]) )
    hFR_data.SetName( "PhFR" )

    with TFileContext(outfname_data, "UPDATE") as tf:
        tf.cd()
        hFR_data.Write(hFR_data.GetName(), ROOT.TObject.kOverwrite)

    to_preserve_FRdata = beautify(cFR_data, hFR_data, logx, logy)
    # Probably python/ROOT gc deletes these graphic objects and they disappear from the canvas. Holding a reference to them seems to fix it
        
    for ext in ["png"]:
        cFR_data.SaveAs( "Plots/PhFR/PhFR_from_{:s}.{:s}".format(sample_data["name"], ext) )
    del cFR_data, to_preserve_FRdata

    ## Estimate = B*C/D ##
    cES_data = ROOT.TCanvas( "cES_data", "Estimate "+sample_data["name"], 200, 10, 1200, 900 )
    cES_data.cd()

    hES_data = hFR_data  # rename for clarity
    hES_data.Multiply(hdata_B)
    hES_data.SetTitle( "Fake Photon estimate: B*C/D (from {:s})".format(sample_data["title"]) )
    hES_data.SetName ( "FakeEstimate" )

    with TFileContext(outfname_data, "UPDATE") as tf:
        tf.cd()
        hES_data.Write(hES_data.GetName(), ROOT.TObject.kOverwrite)

    to_preserve_ESdata = beautify(cES_data, hES_data, logx, logy, logz=True)
    
    for ext in ["png"]:
        cES_data.SaveAs( "Plots/PhFR/PhEstimate_from_{:s}.{:s}".format(sample_data["name"], ext) )
    del cES_data, to_preserve_ESdata, hES_data
    
    ##### Second: data - promptMC #####
    if(sample_prompt is None):
        return

    outfname_full = "data/PhFR_from_{:s}-{:s}.root".format(sample_data["name"], sample_prompt["name"])
    print('Output (with prompt MC subtraction) in: "{:s}"'.format(outfname_full))
    
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
    
    with TFileContext(outfname_full, "UPDATE") as tf:
        tf.cd()
        hFR.Write(hFR.GetName(), ROOT.TObject.kOverwrite)
        
    to_preserve_FR = beautify(cFR, hFR, logx, logy)
        
    for ext in ["png"]: 
        cFR.SaveAs( "Plots/PhFR/PhFR_from_{:s}-{:s}.{:s}".format(sample_data["name"], sample_prompt["name"], ext) )
    del cFR, to_preserve_FR
    
    ## Estimate = B*C/D ##
    cES = ROOT.TCanvas( "cES", "Estimate {} - {}".format(sample_data["name"], sample_prompt["name"]) , 200, 10, 1200, 900 )
    cES.cd()

    hES = hFR  # rename for clarity
    hES.Multiply(hdata_B)
    hES.SetTitle( "Fake Photon estimate: B*C/D (from {:s} - {:s})".format(sample_data["title"], sample_prompt["title"]) )
    hES.SetName( "FakeEstimate" )

    if(fixNegBins):
        for b in iterate_bins(hES):
            if(hES.GetBinContent(b) < 0): hES.SetBinContent(b, 0.)

    with TFileContext(outfname_full, "UPDATE") as tf:
        tf.cd()
        hES.Write(hES.GetName(), ROOT.TObject.kOverwrite)
    
    to_preserve_ESdata = beautify(cES, hES, logx, logy, logz=True)
    
    for ext in ["png"]:
        cES.SaveAs( "Plots/PhFR/PhEstimate_from_{:s}-{:s}.{:s}".format(sample_data["name"], sample_prompt["name"], ext) )
    del cES, to_preserve_ESdata, hES
    print()

    

def kFactor(sample, analyzer='VVGammaAnalyzer', region='CRLFR', logx=False, logy=False, fixNegBins=None):
    print("kFactor: MC sample =", sample)  # sample = nonprompt MC
    
    if(fixNegBins is None): fixNegBins = sample.get("fixNegBins", False)
    
    path_in = "{:s}/2016/{:s}_{:s}".format(get_path_base(), analyzer, region)
    hA, hB, hC, hD = getPlots( path_in, sample['file'], [ "PhFR_%s_nonprompt" % (c) for c in 'ABCD' ])
    assert (hA and hB and hC and hD), "Could not get all 4 histograms!"
    
    hKF = copy.deepcopy(hA)
    hKF.Divide(hB)
    hKF.Divide(hC)
    hKF.Multiply(hD)

    if(fixNegBins): fix_neg_bins(hKF)
    
    hKF.SetTitle("kFactor: {} (A/B)/(C/D)".format(sample["name"]))
    
    with TFileContext("data/PhKF_from_{:s}.root".format(sample["name"]), "RECREATE") as tf:
        tf.cd()
        hKF.Write(hKF.GetName(), ROOT.TObject.kOverwrite)
    
    c = ROOT.TCanvas("c", "kFactor", 1200, 900)

    beautify(c, hKF, logx, logy)

    hKF.Draw("colz texte")
    lines, ticks = linesAndTicks(hKF, logx, logy)
    to_preserve = lines + ticks
    for l in to_preserve:
        l.Draw("same")
    
    for ext in ["png"]:
        c.SaveAs( "Plots/PhFR/kFactor_from_{:s}.{:s}".format(sample["name"], ext) )
    
    del c, to_preserve

if __name__ == "__main__":
    sampleList = {
        "data"     : {"file": 'data'},
        "ZGToLLG"  : {"file": 'ZGToLLG', "title": "Z#gamma_MC", "fixNegBins": True},
        "Drell-Yan": {"file": 'DYJetsToLL_M50', "fixNegBins": True}
    }
    for key, sample in sampleList.items():
        sample.setdefault("name" , key)
        sample.setdefault("title", sample["name"])

    ROOT.gStyle.SetOptStat(0)

    fakeRate(sampleList["data"], sampleList["ZGToLLG"], analyzer="VVGammaAnalyzer", region="CRLFR", logx=True, fixNegBins=True)
    # fakeRate(sampleList["Drell-Yan"], None            , analyzer="VVGammaAnalyzer", region="CRLFR", logx=True, fixNegBins=True)

    kFactor(sampleList["Drell-Yan"], analyzer='VVGammaAnalyzer', region='CRLFR', logx=True, logy=False)
    # kFactor(sampleList["data"]     , analyzer='VVGammaAnalyzer', region='CRLFR', logx=True, logy=False)
