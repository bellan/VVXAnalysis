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
import itertools
import re

from plotUtils import TFileContext, makedirs_ok

ROOT.gStyle.SetPaintTextFormat(".2f")

def getPlots(inputdir, sample, plots, verbose=True):
    fname = path.join(inputdir, sample+".root")
    retrieved = []
    with TFileContext(fname) as rFile:
        for plot in plots:
            h = rFile.Get(plot)
            if(not h):
                if(verbose): print('WARN: Could not get "%s" from "%s"' % (plot, fname))
                retrieved.append(None)
            else:
                retrieved.append( copy.deepcopy(h) )
    return retrieved


def addIfExisting(*args):
    result = None
    for a in [ a for a in args if a is not None ]:
        if result is None:
            result = a
        else:
            result.Add(a)
    return result


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


def get_path_results_base():
    if(environ.get("CMSSW_BASE")): path_base = path.join(environ['CMSSW_BASE'], "src/VVXAnalysis/TreeAnalysis/results")
    else:                          path_base = "results" # "~/Work/Analysys/rsync_results/Updated/EXT"
    assert path.exists(path_base), 'Please call this script from "VVXAnalysis/TreeAnalysis" or do `cmsenv`'
    return path_base


def get_path_results(year, analyzer, region):
    # return "{:s}/{}/{:s}_{:s}".format(get_path_results_base(), year, analyzer, region)
    return path.join(get_path_results_base(), str(year), '_'.join([analyzer, region]))


def fakeRateABCD(sample_data, sample_prompt, analyzer="VVGammaAnalyzer", year=2016, region="CRLFR", logx=False, logy=False, fixNegBins=False):
    print("Fake Rate ABCD: data     =", sample_data  )
    print("Fake Rate ABCD: promptMC =", sample_prompt)
    path_in = get_path_results(year=year, analyzer=analyzer, region=region)
    outfname = "data/ABCD_FR_{dataname:s}-{promptname:s}_{year}.root".format(dataname=sample_data["name"], promptname=sample_prompt["name"], year=year)
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
        cFR.SaveAs( "Plot/PhFR/ABCD_FR_{dataname:s}-{promptname:s}_{year}.{ext:s}".format(dataname=sample_data["name"], promptname=sample_prompt["name"], year=year, ext=ext) )
    del cFR, to_preserve_FR
    
    ## Estimate = B*C/D ##
    cES = ROOT.TCanvas( "cES", "Estimate {} - {}".format(sample_data["name"], sample_prompt["name"]) , 200, 10, 1200, 900 )
    cES.cd()

    hES = copy.deepcopy(hFR)
    hES.Multiply(hdata_B)
    hES.SetTitle( "Fake Photon estimate: B*C/D (from {:s} - {:s})".format(sample_data["title"], sample_prompt["title"]) )
    hES.SetName( "FakeEstimate" )

    if(fixNegBins):
        fix_neg_bins(hES)

    with TFileContext(outfname, "UPDATE") as tf:
        tf.cd()
        hES.Write(hES.GetName(), ROOT.TObject.kOverwrite)
    
    to_preserve_ESdata = beautify(cES, hES, logx, logy, logz=True)
    
    for ext in ["png"]:
        cES.SaveAs( "Plot/PhFR/ABCD_Estimate_{dataname:s}-{promptname:s}_{year}.{ext:s}".format(dataname=sample_data["name"], promptname=sample_prompt["name"], year=year, ext=ext) )
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
    cKF.SaveAs( 'Plot/PhFR/ABCD_kFactor_{dataname:s}-{promptname:s}_{year}.{ext:s}'.format(dataname=sample_data['name'], promptname=sample_prompt['name'], year=year, ext='png') )
    del cKF, to_preserve_KF
    
    return hFR, hES, hKF


    return hFR, hES


def fakeRateABCD_noSubtract(sample, analyzer="VVGammaAnalyzer", year=2016, region="CRLFR", logx=False, logy=False, fixNegBins=False):
    print("Fake Rate ABCD: sample   =", sample  )
    path_in = get_path_results(year=year, analyzer=analyzer, region=region)
    outfname = "data/ABCD_FR_{samplename:s}_{year}.root".format(samplename=sample["name"], year=year)
    print('Output (without MC subtraction) in: "{:s}"'.format(outfname))
    
    hAp, hBp, hCp, hDp = getPlots(path_in, sample["file"], [ "PhFR_%s_nonprompt" % (c) for c in 'ABCD' ] )
    hAn, hBn, hCn, hDn = getPlots(path_in, sample["file"], [ "PhFR_%s_prompt"    % (c) for c in 'ABCD' ] )
    hA = addIfExisting(hAp, hAn)
    hB = addIfExisting(hAp, hBn)
    hC = addIfExisting(hAp, hCn)
    hD = addIfExisting(hAp, hDn)
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
        cFR.SaveAs( "Plot/PhFR/ABCD_FR_{samplename:s}_{year}.{ext:s}".format(samplename=sample["name"], year=year, ext=ext) )
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
        cES.SaveAs( "Plot/PhFR/ABCD_Estimate_{samplename:s}_{year}.{ext:s}".format(samplename=sample["name"], year=year, ext=ext) )
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
    cKF.SaveAs( 'Plot/PhFR/ABCD_kFactor_{samplename:s}_{year}.{ext:s}'.format(samplename=sample['name'], year=year, ext='png') )
    del cKF, to_preserve_KF
    
    return hFR, hES, hKF


def getPassFailLtoT_noSubtract_regex(sample, method, regex, year, analyzer, region, fixNegBins=False):
    regex_compiled = re.compile(regex)
    path_in = get_path_results(year=year, analyzer=analyzer, region=region)
    channels = findChannelsInFile(path.join(path_in, sample['file']+'.root'), method)

    print('\tDEBUG: regex    =', regex_compiled.pattern)
    selected = [ c  for c in channels if regex_compiled.search(c) ]
    print('\tDEBUG: selected =', selected)

    if(method.startswith('LtoT')):
        raise NotImplementedError('Use VLtoL')

    listPASS = []
    listFAIL = []
    for channel in selected:
        hPASSp, hFAILp = getPlots(path_in, sample['file'], [ 'PhFR_%s_%s_prompt_%s'    % (method, channel, s) for s in ['PASS', 'FAIL'] ] , verbose=False)
        hPASSn, hFAILn = getPlots(path_in, sample['file'], [ 'PhFR_%s_%s_nonprompt_%s' % (method, channel, s) for s in ['PASS', 'FAIL'] ] , verbose=False)

        listPASS.append( addIfExisting(hPASSp, hPASSn) )
        listFAIL.append( addIfExisting(hFAILp, hFAILp) )

    hPASS = addIfExisting(*listPASS)
    hFAIL = addIfExisting(*listFAIL)

    if(sample.get('fixNegBins') and fixNegBins):
        fix_neg_bins(hPASS)
        fix_neg_bins(hFAIL)

    return hPASS, hFAIL


def getPassFailLtoT_regex(sample_data, sample_prompt, method, regex, year, analyzer, region, fixNegBins=False):
    regex_compiled = re.compile(regex)
    path_in = get_path_results(year=year, analyzer=analyzer, region=region)
    channels = findChannelsInFile(path.join(path_in, sample_data['file']+'.root'), method)

    print('\tDEBUG: regex    =', regex_compiled.pattern)
    selected = [ c  for c in channels if regex_compiled.search(c) ]
    print('\tDEBUG: selected =', selected)

    if(method.startswith('LtoT')):
        raise NotImplementedError('Use VLtoL')

    listDataPASS   = []
    listDataFAIL   = []
    listPromptPASS = []
    listPromptFAIL = []
    for channel in selected:
        hdata_PASS  , hdata_FAIL   = getPlots(path_in, sample_data['file']  , [ 'PhFR_%s_%s_%s'        % (method, channel, s) for s in ['PASS', 'FAIL'] ] , verbose=True)
        hprompt_PASS, hprompt_FAIL = getPlots(path_in, sample_prompt['file'], [ 'PhFR_%s_%s_prompt_%s' % (method, channel, s) for s in ['PASS', 'FAIL'] ] , verbose=True)

        listDataPASS  .append(hdata_PASS)
        listDataFAIL  .append(hdata_FAIL)
        listPromptPASS.append(hprompt_PASS)
        listPromptFAIL.append(hprompt_FAIL)

    hDataPASS   = addIfExisting(*listDataPASS)
    hDataFAIL   = addIfExisting(*listDataFAIL)
    hPromptPASS = addIfExisting(*listPromptPASS)
    hPromptFAIL = addIfExisting(*listPromptFAIL)

    assert hDataPASS   and hDataFAIL  , 'Could not get the 2 data sums of histograms'
    assert hPromptPASS and hPromptFAIL, 'Could not get the 2 MC sums of histograms'

    if(sample_data.get('fixNegBins') and fixNegBins):
        fix_neg_bins(hDataPASS)
        fix_neg_bins(hDataFAIL)
    if(sample_prompt.get('fixNegBins') and fixNegBins):
        fix_neg_bins(hPromptPASS)
        fix_neg_bins(hPromptFAIL)

    hPASS = hDataPASS
    hPASS.Add(hPromptPASS, -1)
    hFAIL = hDataFAIL
    hFAIL.Add(hPromptFAIL, -1)

    return hPASS, hFAIL


def fakeRateLtoT_regex(sample_data, sample_prompt, method, regex, year=2016, analyzer='VVGammaAnalyzer', region='CRLFR', pattern_printable=None, logx=False, logy=False, fixNegBins=False):
    if(pattern_printable is None):
        pattern_printable = regex.replace('\\','').replace('"','').replace("'","")

    if(sample_prompt is None):
        print('Fake Rate {}: sample   ='.format(method), sample_data)
        samplename  = sample_data['name']
        sampletitle = sample_data['title']

        hPASS, hFAIL = getPassFailLtoT_noSubtract_regex(sample_data, method=method, regex=regex, year=year, analyzer=analyzer, region=region, fixNegBins=fixNegBins)
    else:
        print('Fake Rate {}_{}: data     ='.format(method, pattern_printable), sample_data  )
        print('Fake Rate {}_{}: promptMC ='.format(method, pattern_printable), sample_prompt)
        samplename  = sample_data['name']  +  '-'  + sample_prompt['name']
        sampletitle = sample_data['title'] + ' - ' + sample_prompt['title']

        hPASS, hFAIL = getPassFailLtoT_regex(sample_data, sample_prompt, method=method, regex=regex, year=year, analyzer=analyzer, region=region, fixNegBins=fixNegBins)

    outname = 'FR_{method}_{regex}_{samplename}{region}_{year}'.format(method=method, samplename=samplename, region='' if region=='CRLFR' else '_'+region, year=year, regex=pattern_printable)
    title = 'Photon Fake Rate: {method} {regex} (from {sampletitle} in {region})'.format(method=method, sampletitle=sampletitle, regex=pattern_printable, region=region)

    hTOTAL = hFAIL
    hTOTAL.Add(hPASS)

    ## Fake rate = PASS/TOTAL ##
    hFR = hPASS
    hFR.Divide(hTOTAL)

    plotFR_LtoT(hFR, outname, title, logx=logx, logy=logy)

    return hFR


def getPassFailLtoT_noSubtract(sample, analyzer, year, region, method, fixNegBins=False):
    path_in = get_path_results(year=year, analyzer=analyzer, region=region)
    if(method.startswith('LtoT')):
        if(sample['name'] == 'ZZTo4l' and region=='CRLFR'):
            h5p, h4p, h3p = getPlots(path_in, sample['file'], [ 'PhFR_%s_%d_prompt'    % (method, s) for s in [5, 4, 3] ] )
            h5n, h4n, h3n = getPlots(path_in, sample['file'], [ 'PhFR_%s_%d_nonprompt' % (method, s) for s in [5, 4, 3] ] )
        else:
            h5p, h4p, h3p = getPlots(path_in, sample['file'], [ 'PhFR_%s_%d_prompt'    % (method, s) for s in [5, 4, 3] ] )
            h5n, h4n, h3n = getPlots(path_in, sample['file'], [ 'PhFR_%s_%d_nonprompt' % (method, s) for s in [5, 4, 3] ] )
        
        hPASS = addIfExisting(h5p, h5n)
        hFAIL = addIfExisting(h3p, h3n, h4p, h4n)
        if(not hPASS): print('Could not get the PASS LtoT histogram')
        if(not hFAIL): print('Could not get any of the the FAIL LtoT histogram')
        assert (hPASS and hFAIL), "Could't get the 2*3 LtoT MC histograms!"
    else:
        hPASSp, hFAILp = getPlots(path_in, sample['file'], [ 'PhFR_%s_prompt_%s'    % (method, s) for s in ['PASS', 'FAIL'] ] )
        hPASSn, hFAILn = getPlots(path_in, sample['file'], [ 'PhFR_%s_nonprompt_%s' % (method, s) for s in ['PASS', 'FAIL'] ] )
        
        hPASS = addIfExisting(hPASSp, hPASSn)
        hFAIL = addIfExisting(hFAILp, hFAILn)
        if(not hPASS): print('Could not get the PASS histogram')
        if(not hFAIL): print('Could not get the FAIL histogram')
        assert (hPASS and hFAIL), "Could't get the 2*2 MC histograms!"

    if(sample.get('fixNegBins', False) and fixNegBins):
        for h in [hPASS, hFAIL]:
            fix_neg_bins(h)

    return hPASS, hFAIL


def getPassFailLtoT(sample_data, sample_prompt, analyzer, year, region, method, fixNegBins):
    path_in = get_path_results(year=year, analyzer=analyzer, region=region)
    
    if(method.startswith('LtoT')):
        hdata_5  , hdata_4  , hdata_3   = getPlots(path_in, sample_data['file']  , [ 'PhFR_LtoT_%d'        % (s) for s in [5, 4, 3] ] )
        hprompt_5, hprompt_4, hprompt_3 = getPlots(path_in, sample_prompt['file'], [ 'PhFR_LtoT_%d_prompt' % (s) for s in [5, 4, 3] ] )
        assert (hdata_5   and (hdata_4   or  hdata_3  )), "Could't get the 3 data histograms!"
        assert (hprompt_5 and (hprompt_4 or  hprompt_3)), "Could't get the 3 prompt MC histograms"

        hdata_PASS   = hdata_5
        hdata_FAIL   = addIfExisting(hdata_3  , hdata_4  )
        hprompt_PASS = hprompt_5
        hprompt_FAIL = addIfExisting(hprompt_3, hprompt_4)

    else:
        hdata_PASS  , hdata_FAIL   = getPlots(path_in, sample_data['file']  , [ 'PhFR_%s_%s'        % (method, s) for s in ['PASS', 'FAIL'] ] )
        hprompt_PASS, hprompt_FAIL = getPlots(path_in, sample_prompt['file'], [ 'PhFR_%s_prompt_%s' % (method, s) for s in ['PASS', 'FAIL'] ] )
        
        assert (hdata_PASS   and hdata_FAIL  ), "Could't get the 2 data histograms!"
        assert (hprompt_PASS and hprompt_FAIL), "Could't get the 2 prompt MC histograms!"

    if(sample_data.get('fixNegBins', False) and fixNegBins):
        for h in [hdata_PASS, hdata_FAIL]:
            fix_neg_bins(h)
    if(sample_prompt.get('fixNegBins', False) and fixNegBins):
        for h in (hprompt_PASS, hprompt_FAIL):
            fix_neg_bins(h)

    hdata_PASS.Add(hprompt_PASS, -1)
    hdata_FAIL.Add(hprompt_FAIL, -1)

    return hdata_PASS, hdata_FAIL


def fakeRateLtoT(sample_data, sample_prompt, analyzer='VVGammaAnalyzer', year=2016, region='CRLFR', method='LtoT', logx=False, logy=False, fixNegBins=False):
    if(sample_prompt is None):
        print('Fake Rate {}: sample   ='.format(method), sample)
        samplename  = sample_data['name']
        sampletitle = sample_data['title']

        hPASS, hFAIL = getPassFailLtoT_noSubtract(sample_data, analyzer=analyzer, year=year, region=region, method=method, fixNegBins=fixNegBins)
    else:
        print('Fake Rate {}: data     ='.format(method), sample_data  )
        print('Fake Rate {}: promptMC ='.format(method), sample_prompt)
        samplename  = sample_data['name']  +  '-'  + sample_prompt['name']
        sampletitle = sample_data['title'] + ' - ' + sample_prompt['title']

        hPASS, hFAIL = getPassFailLtoT(sample_data=sample_data, sample_prompt=sample_prompt, analyzer=analyzer, year=year, region=region, method=method, fixNegBins=fixNegBins)

    outname = '{method}_FR_{samplename}{region}_{year}'.format(method=method, samplename=samplename, region='' if region=='CRLFR' else '_'+region, year=year)
    if(method.startswith('LtoT')):
        title = 'Photon Fake Rate: 5/(3+4+5) (from {sampletitle} in {region})'.format(sampletitle=sampletitle, region=region)
    else:
        title = 'Photon Fake Rate: {method} (from {sampletitle} in {region})'.format(method=method, sampletitle=sampletitle, region=region)

    hTOTAL = hFAIL
    hTOTAL.Add(hPASS)

    ## Fake rate = PASS/TOTAL ##
    hFR = hPASS
    hFR.Divide(hTOTAL)

    plotFR_LtoT(hFR, outname, title, logx=logx, logy=logy)

    return hFR


def plotFR_LtoT(hFR, outname, title, logx=False, logy=False):
    outfname = 'data/{}.root'.format(outname)
    picname = 'Plot/PhFR/{}'.format(outname)

    cFR = ROOT.TCanvas( 'cFR_{}'.format(outname), title, 1200, 900 )
    cFR.cd()
    
    hFR.SetTitle(title)
    hFR.SetName( 'PhFR' )
    
    with TFileContext(outfname, 'UPDATE') as tf:
        tf.cd()
        hFR.Write(hFR.GetName(), ROOT.TObject.kOverwrite)
    print('Output (with prompt MC subtraction) in: "{:s}"'.format(outfname))

    to_preserve_FR = beautify(cFR, hFR, logx, logy)

    for ext in ['png']:
        cFR.SaveAs( picname+'.'+ext )
    del cFR, to_preserve_FR


def plotRatio(h1, h2, picture_name="ratio", title="ratio"):
    assert h1, "h1 missing"
    assert h2, "h2 missing"

    ratio = copy.deepcopy(h1)
    ratio.Divide(h2)
    ratio.SetName('PhFRSF')
    ratio.SetTitle(title)
    
    c = ROOT.TCanvas("cratio_{:s}".format(picture_name), "ratio", 1200, 900)
    c.cd()
    
    ratio.Draw("colz texte")
    to_preserve = beautify(c, ratio, True, False)
    # ratio.GetYaxis().SetRange(1, ratio.GetXaxis().GetNbins())
    c.SaveAs('Plot/PhFR/{:s}.png'.format(picture_name))
    with TFileContext(path.join('data', picture_name+'.root'), 'RECREATE') as tf:
        tf.cd()
        ratio.Write()
        print('Output (ratio) in "{:s}"'.format(tf.GetName()))
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
        def get_bin_2d(h, x, y):
            return h.GetBin(x, y)
    elif(direction=='X'):
        axis_2D_proj = h2.GetXaxis()
        axis_2D_draw = h2.GetYaxis()
        nbins_proj = h2.GetNbinsX()
        nbins_draw = h2.GetNbinsY()
        def get_bin_2d(h, x, y):
            return h.GetBin(y, x)
    
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
            bin_2d = get_bin_2d(h2, i, j)
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
    h1s[0].Draw("P0 E1")
    for h1 in h1s[1:]:
        h1.Draw("P0 E1 same")
    
    legend.Draw("SAME")
    
    cprof.SaveAs( "Plot/PhFR/{:s}.png".format(picture_name) )


def findChannelsInFile(fname, method):
    regex_compiled = re.compile('PhFR_{method}_([^_]+)_((non)?prompt_)?(PASS|FAIL)'.format(method=method))

    with TFileContext(fname) as tf:
        keys = [ k.GetName() for k in tf.GetListOfKeys() ]

    matches = []
    for key in keys:
        m = regex_compiled.match(key)
        if(m):
            matches.append(m.group(1))  # group(0) is the whole expression

    if(len(matches) == 0):
        print("\tERROR: no key matching the following regex:", regex_compiled.pattern,
              "\n\t       Keys starting with PhFR:", [k for k in keys if k.startswith('PhFR')])
    return matches


if __name__ == "__main__":
    sampleList = {
        "data"     : {"file": 'data'},
        "ZGToLLG"  : {"file": 'ZGToLLG', "title": "Z#gamma_{MC}", "fixNegBins": True},
        "Drell-Yan": {"file": 'DYJetsToLL_M50', "fixNegBins": True},
        "ZZTo4l"   : {"file": 'ZZTo4l', "fixNegBins": True},
        "ggTo4l"   : {"file": 'ggTo4l', "fixNegBins": True},
        "WZ"       : {"file": "WZTo3LNu", "fixNegBins":True}
    }
    for key, sample in sampleList.items():
        sample.setdefault("name" , key)
        sample.setdefault("title", sample["name"])

    possible_methods = {"LtoT", "VLtoL", "KtoVL", "KtoVLexcl", "ABCD"}

    from argparse import ArgumentParser
    parser = ArgumentParser()
    parser.add_argument("-y", "--year"   , type=int, default=2016, choices=[2016, 2017, 2018])
    parser.add_argument("-m", "--method" , nargs="+", choices=possible_methods, default=["LtoT"])
    parser.add_argument(      "--channels", action='store_true')
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

    for method in set(args.method) & possible_methods:
        if(args.do_data):
            print('>>>', method, 'data')
            if(args.channels):
                hDe  = fakeRateLtoT_regex(sampleList['data'], sampleList["ZGToLLG"], year=args.year, method=method, regex='2[me]\+e[PF]', pattern_printable='2x+e', logx=True,logy=False,fixNegBins=False)
                hDm  = fakeRateLtoT_regex(sampleList['data'], sampleList["ZGToLLG"], year=args.year, method=method, regex='2[me]\+m[PF]', pattern_printable='2x+m', logx=True,logy=False,fixNegBins=False)
                hDP  = fakeRateLtoT_regex(sampleList['data'], sampleList["ZGToLLG"], year=args.year, method=method, regex='2[me]\+[em]P', pattern_printable='2x+P', logx=True,logy=False,fixNegBins=False)
                hDF  = fakeRateLtoT_regex(sampleList['data'], sampleList["ZGToLLG"], year=args.year, method=method, regex='2[me]\+[em]F', pattern_printable='2x+F', logx=True,logy=False,fixNegBins=False)

                hDeP = fakeRateLtoT_regex(sampleList['data'], sampleList["ZGToLLG"], year=args.year, method=method, regex='2[me]\+eP'   , pattern_printable='2x+eP', logx=True, fixNegBins=True)
                hDeF = fakeRateLtoT_regex(sampleList['data'], sampleList["ZGToLLG"], year=args.year, method=method, regex='2[me]\+eF'   , pattern_printable='2x+eF', logx=True, fixNegBins=True)
                hDmP = fakeRateLtoT_regex(sampleList['data'], sampleList["ZGToLLG"], year=args.year, method=method, regex='2[me]\+mP'   , pattern_printable='2x+mP', logx=True, fixNegBins=True)
                hDmF = fakeRateLtoT_regex(sampleList['data'], sampleList["ZGToLLG"], year=args.year, method=method, regex='2[me]\+mF'   , pattern_printable='2x+mF', logx=True, fixNegBins=True)

                hD2e = fakeRateLtoT_regex(sampleList['data'], sampleList["ZGToLLG"], year=args.year, method=method, regex='2e\+[em][PF]', pattern_printable='2e+x' , logx=True, fixNegBins=True)
                hD2m = fakeRateLtoT_regex(sampleList['data'], sampleList["ZGToLLG"], year=args.year, method=method, regex='2m\+[em][PF]', pattern_printable='2m+x' , logx=True, fixNegBins=True)

                plotRatio(hDe , hDm , picture_name="ratio_{}_data-ZG_e_over_m_{}".format(  method, args.year), title="Ratio FR(2x+e)/FR(2x+m) in data-Z#gamma"  )
                plotRatio(hDF , hDP , picture_name="ratio_{}_data-ZG_F_over_P_{}".format(  method, args.year), title="Ratio FR(2x+F)/FR(2x+P) in data-Z#gamma"  )
                plotRatio(hDeF, hDeP, picture_name="ratio_{}_data-ZG_eF_over_eP_{}".format(method, args.year), title="Ratio FR(2x+eF)/FR(2x+eP) in data-Z#gamma")
                plotRatio(hDmF, hDmP, picture_name="ratio_{}_data-ZG_mF_over_mP_{}".format(method, args.year), title="Ratio FR(2x+mF)/FR(2x+mP) in data-Z#gamma")
                plotRatio(hD2e, hD2m, picture_name="ratio_{}_data-ZG_2e_over_2m_{}".format(method, args.year), title="Ratio FR(2e+x)/FR(2m+x) in data-Z#gamma"  )
            else:
                hFR_LtoT_data = fakeRateLtoT(sampleList["data"], sampleList["ZGToLLG"], year=args.year, method=method, logx=True, fixNegBins=False)
                plotProfiled(hFR_LtoT_data, picture_name='LtoT_FR_profiledX_data-ZGToLLG', title='FR(#gamma) vs #eta' , direction='X')
                plotProfiled(hFR_LtoT_data, picture_name='LtoT_FR_profiledY_data-ZGToLLG', title='FR(#gamma) vs p_{T}', direction='Y')
            print()

        if(args.do_mc):
            print('>>>', method, 'MC')
            if(args.channels):
                he  = fakeRateLtoT_regex(sampleList['ZZTo4l'], None, method=method, regex='2[me]\+e[PF]', year=args.year, pattern_printable='2x+e' , logx=True, fixNegBins=True)
                hm  = fakeRateLtoT_regex(sampleList['ZZTo4l'], None, method=method, regex='2[me]\+m[PF]', year=args.year, pattern_printable='2x+m' , logx=True, fixNegBins=True)
                hP  = fakeRateLtoT_regex(sampleList['ZZTo4l'], None, method=method, regex='2[me]\+[em]P', year=args.year, pattern_printable='2x+P' , logx=True, fixNegBins=True)
                hF  = fakeRateLtoT_regex(sampleList['ZZTo4l'], None, method=method, regex='2[me]\+[em]F', year=args.year, pattern_printable='2x+F' , logx=True, fixNegBins=True)
                heP = fakeRateLtoT_regex(sampleList['ZZTo4l'], None, method=method, regex='2[me]\+eP'   , year=args.year, pattern_printable='2x+eP', logx=True, fixNegBins=True)
                heF = fakeRateLtoT_regex(sampleList['ZZTo4l'], None, method=method, regex='2[me]\+eF'   , year=args.year, pattern_printable='2x+eF', logx=True, fixNegBins=True)
                hmP = fakeRateLtoT_regex(sampleList['ZZTo4l'], None, method=method, regex='2[me]\+mP'   , year=args.year, pattern_printable='2x+mP', logx=True, fixNegBins=True)
                hmF = fakeRateLtoT_regex(sampleList['ZZTo4l'], None, method=method, regex='2[me]\+mF'   , year=args.year, pattern_printable='2x+mF', logx=True, fixNegBins=True)
                h2e = fakeRateLtoT_regex(sampleList['ZZTo4l'], None, method=method, regex='2e\+[em][PF]', year=args.year, pattern_printable='2e+x' , logx=True, fixNegBins=True)
                h2m = fakeRateLtoT_regex(sampleList['ZZTo4l'], None, method=method, regex='2m\+[em][PF]', year=args.year, pattern_printable='2m+x' , logx=True, fixNegBins=True)

                plotRatio(h2e, h2m, picture_name="ratio_{}_ZZ_2e_over_2m_{}".format(method, args.year), title="Ratio FR(2e+x)/FR(2m+x) in ZZTo4l")
                plotRatio(he , hm , picture_name="ratio_{}_ZZ_e_over_m_{}".format(  method, args.year), title="Ratio FR(2x+e)/FR(2x+m) in ZZTo4l")
                plotRatio(hF , hP , picture_name="ratio_{}_ZZ_F_over_P_{}".format(  method, args.year), title="Ratio FR(2x+F)/FR(2x+P) in ZZTo4l")
                plotRatio(heF, heP, picture_name="ratio_{}_ZZ_eF_over_eP_{}".format(method, args.year), title="Ratio FR(2x+eF)/FR(2x+eP) in ZZTo4l")
                plotRatio(hmF, hmP, picture_name="ratio_{}_ZZ_mF_over_mP_{}".format(method, args.year), title="Ratio FR(2x+mF)/FR(2x+mP) in ZZTo4l")
            else:
                hFR_LtoT_DY   = fakeRateLtoT(sampleList["Drell-Yan"]       , None, year=args.year, method=method, logx=True, fixNegBins=True)
                hFR_LtoT_ZG   = fakeRateLtoT(sampleList["ZGToLLG"]         , None, year=args.year, method=method, logx=True, fixNegBins=True)
                hFR_LtoT_ZZ   = fakeRateLtoT(sampleList["ZZTo4l"]          , None, year=args.year, method=method, logx=True, fixNegBins=True)
                # hFR_LtoT_ZZ_4P= fakeRateLtoT(sampleList["ZZTo4l"]          , None, year=args.year, method=method, logx=True, fixNegBins=True, region='SR4P')
                hFR_LtoT_gg   = fakeRateLtoT(sampleList["ggTo4l"]          , None, year=args.year, method=method, logx=True, fixNegBins=True)

                # plotRatio(hFR_LtoT_ZZ_4P, hFR_LtoT_ZZ, picture_name="{}_ratio_ZZ4P_over_ZZCRLFR_{}".format(method, args.year), title="Ratio FR(ZZ_{4P})/FR(ZZ_{CRLFR})")
            print()

        if(args.do_data and args.do_mc and not args.channels):
            print('>>>', method, 'ratio')
            ratioLtoT = plotRatio(hFR_LtoT_data, hFR_LtoT_DY, picture_name="{}_ratio_data-ZG_over_DY_{}".format(method, args.year), title="Ratio FR(data-Z#gamma)/FR(DY)")
            plotRatio(hFR_LtoT_data, hFR_LtoT_ZZ, picture_name="{}_ratio_data-ZG_over_ZZ_{}".format(method, args.year), title="Ratio FR(data-Z#gamma)/FR(ZZ)")
            # ratioLtoT_4P = plotRatio(hFR_LtoT_data, hFR_LtoT_ZZ_4P, picture_name="{}_ratio_data-ZG_over_ZZ_4P_{}".format(method, args.year), title="Ratio FR(data-Z#gamma)_{CRLFR}/FR(ZZ_{4P})")
            # ratioLtoT = plotRatio(hFR_LtoT_data, hFR_LtoT_gg, picture_name="{}_ratio_data-ZG_over_gg_{}".format(method, args.year), title="Ratio FR(data-Z#gamma)/FR(gg)")
            print()
