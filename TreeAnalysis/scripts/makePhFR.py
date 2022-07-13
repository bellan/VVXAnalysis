#!/usr/bin/env python

from __future__ import print_function
import copy, sys, os
import ROOT
if '-b' in sys.argv:
    ROOT.gROOT.SetBatch(True)
sys.path.append("VVXAnalysis/TreeAnalysis/python")
# from readSampleInfo import *
# from Colours import *
# import ctypes
import re
import math


class TFileContext(object):
    def __init__(self, *args):
        self.tfile = ROOT.TFile(*args)
    def __enter__(self):
        return self.tfile
    def __exit__(self, type, value, traceback):
        self.tfile.Close()


def getPlots(inputdir, sample, plots):
    fname = os.path.join(inputdir, sample+".root")
    retrieved = []
    with TFileContext(fname) as rFile:
        for plot in plots:
            h = rFile.Get(plot)
            if(not h):
                print('Warning: Could not get "%s" from "%s"' % (plot, fname))
                retrieved.append(None)
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


def linesAndTicks(h, logx=False, logy=False):
    lines = []
    ticks = []
    xmin, xmax = h.GetXaxis().GetXmin(), h.GetXaxis().GetXmax()
    ymin, ymax = h.GetYaxis().GetXmin(), h.GetYaxis().GetXmax()
    xsize = xmin + 0.03*(xmax - xmin)  # Ticks are 3 % of axis length
    if(logx): xsize = xmin + 10**(0.03*math.log10(xmax/xmin))
    ysize = ymin + 0.03*(ymax - ymin)
    if(logy): ysize = ymin + 10**(0.03*math.log10(ymax/ymin))
        
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

def fakeRate(sample_data, sample_prompt, analyzer="VVXAnalyzer", region="CRLFR", logx=False, logy=False, fixNegBins=False):
    print("Fake Rate: data     =", sample_data  )
    print("Fake Rate: promptMC =", sample_prompt)
    
    path_base = "results"  # os.environ["CMSSW_BASE"]+"src/VVXAnalysis/TreeAnalysis/results"
    path = "{:s}/2016/{:s}_{:s}".format(path_base, analyzer, region)
    
    hdata_A  , hdata_B  , hdata_C  , hdata_D   = getPlots(path, sample_data  ["file"], [ "PhFR_%s" % (c) for c in 'ABCD' ] )
    assert (hdata_A   and hdata_B   and hdata_C   and hdata_D  ), "Could't get all 4 data histograms!"

    if(sample_data.get("fixNegBins", False)):
        for h in (hdata_A, hdata_B, hdata_C, hdata_D):
            for b in iterate_bins(h):
                if(h.GetBinContent(b) < 0): h.SetBinContent(b, 0.)

    ##### First: only data ####

    ## Fake rate = C/D ##
    cFR_data = ROOT.TCanvas( "cFR_data", "Fake Rate "+sample_data["name"], 200, 10, 1200, 900 )
    cFR_data.cd()
    
    hFR_data = copy.deepcopy(hdata_C)
    # hFR_data.Sumw2()
    hFR_data.Divide(hdata_D)
    hFR_data.SetTitle( "Photon Fake Rate: C/D (from {})".format(sample_data["title"]) )

    to_preserve_FRdata = beautify(cFR_data, hFR_data, logx, logy)
    # Probably python/ROOT gc deletes these graphic objects and they disappear from the canvas. Holding a reference to them seems to fix it
        
    for ext in ["png"]:
        cFR_data.SaveAs( "Plot/PhFR_from_{:s}.{:s}".format(sample_data["name"], ext) )
    del cFR_data, to_preserve_FRdata

    ## Estimate = B*C/D ##
    cES_data = ROOT.TCanvas( "cES_data", "Estimate "+sample_data["name"], 200, 10, 1200, 900 )
    cES_data.cd()

    hES_data = hFR_data  # rename for clarity
    hES_data.Multiply(hdata_B)
    hES_data.SetTitle( "Fake Photon estimate: B*C/D (from {:s})".format(sample_data["title"]) )

    to_preserve_ESdata = beautify(cES_data, hES_data, logx, logy, logz=True)
    
    for ext in ["png"]:
        cES_data.SaveAs( "Plot/PhEstimate_from_{:s}.{:s}".format(sample_data["name"], ext) )
    del cES_data, to_preserve_ESdata, hES_data

    
    ##### Second: data - promptMC #####
    if(sample_prompt is None):
        return

    hprompt_A, hprompt_B, hprompt_C, hprompt_D = getPlots(path, sample_prompt["file"], [ "PhFR_%s" % (c) for c in 'ABCD' ] )
    assert (hprompt_A and hprompt_B and hprompt_C and hprompt_D), "Could't get all 4 promptMC histograms!"

    if(sample_prompt.get("fixNegBins", False)):
        for h in (hprompt_A, hprompt_B, hprompt_C, hprompt_D):
            for b in iterate_bins(h):
                if(h.GetBinContent(b) < 0):
                    h.SetBinContent(b, 0.)

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
    
    if(fixNegBins):
        for b in iterate_bins(hFR):
            if(hFR.GetBinContent(b) < 0): hFR.SetBinContent(b, 0.)
    
    with TFileContext("Plot/FakeRate_from_{:s}-{:s}.root".format(sample_data["name"], sample_prompt["name"]), "RECREATE") as tf:
        tf.cd()
        hFR.Write()
        
    to_preserve_FR = beautify(cFR, hFR, logx, logy)
        
    for ext in ["png"]: 
        cFR.SaveAs( "Plot/PhFR_from_{:s}-{:s}.{:s}".format(sample_data["name"], sample_prompt["name"], ext) )
    del cFR, to_preserve_FR
    
    ## Estimate = B*C/D ##
    cES = ROOT.TCanvas( "cES", "Estimate {} - {}".format(sample_data["name"], sample_prompt["name"]) , 200, 10, 1200, 900 )
    cES.cd()

    hES = hFR  # rename for clarity
    hES.Multiply(hdata_B)
    hES.SetTitle( "Fake Photon estimate: B*C/D (from {:s} - {:s})".format(sample_data["title"], sample_prompt["title"]) )

    if(fixNegBins):
        for b in iterate_bins(hES):
            if(hES.GetBinContent(b) < 0): hES.SetBinContent(b, 0.)

    with TFileContext("Plot/Estimate_from_{:s}-{:s}.root".format(sample_data["name"], sample_prompt["name"]), "RECREATE") as tf:
        tf.cd()
        hES.Write()
    
    to_preserve_ESdata = beautify(cES, hES, logx, logy, logz=True)
    
    for ext in ["png"]:
        cES.SaveAs( "Plot/PhEstimate_from_{:s}-{:s}.{:s}".format(sample_data["name"], sample_prompt["name"], ext) )
    del cES, to_preserve_ESdata, hES

    

def kFactor(sample, logx=False, logy=False):
    print("kFactor: MC sample =", sample)
    
    with TFileContext("Plot/Estimate_from_{:s}.root".format(sample["name"]), "READ") as tf_MC, TFileContext("Plot/Estimate_from_data.root", "READ") as tf_data:
        hdata = tf_data.Get("fake_photon_estimate")
        hMC   = tf_MC  .Get("fake_photon_estimate")
        if(not hMC or not hdata):
            print("no histograms:", hMC, hdata)
            return 1
        hKF = copy.deepcopy(hdata)
        # hKF.Sumw2()
        hKF.Divide(hMC)

    hKF.SetTitle("kFactor: data / "+sample["name"])
    c = ROOT.TCanvas("c", "kFactor", 50, 100, 1200, 900)
    if(logx):
        c.SetLogx()
        hKF.GetXaxis().SetMoreLogLabels()
        hKF.GetXaxis().SetNoExponent()
    if(logy):
        c.SetLogy()
        hKF.GetYaxis().SetNoExponent()
        hKF.GetYaxis().SetMoreLogLabels()
    # if(logz):
    #     c.SetLogz()
    #     hKF.GetZaxis().SetNoExponent()
    #     hKF.GetZaxis().SetMoreLogLabels()

    hKF.Draw("colz texte")
    lines, ticks = linesAndTicks(hKF, logx, logy)
    to_preserve = lines + ticks
    for l in to_preserve:
        l.Draw("same")
    
    for ext in ["png"]:
        c.SaveAs( "Plot/kFactor_from_{:s}.{:s}".format(sample["name"], ext) )
    
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

    fakeRate(sampleList["data"], sampleList["ZGToLLG"], analyzer="VVGammaAnalyzer", logx=True, region="CRLFR", fixNegBins=True)
    fakeRate(sampleList["Drell-Yan"], None            , analyzer="VVGammaAnalyzer", logx=True, region="CRLFR", fixNegBins=True)

    # kFactor(sampleList["Drell-Yan"], logx=True, logy=False)
