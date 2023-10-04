#!/usr/bin/env python

##############################################
#                                         
# Creates files for photon fake rate
# 
##############################################

from __future__ import print_function
import sys
from os import path, environ, makedirs
import copy
import ROOT
# from readSampleInfo import *
# from Colours import *
from math import log10
from ctypes import c_int, c_double
from array import array
import itertools
import re
from plotUtils23 import TFileContext, addIfExisting, rebin2D
import Colours

if(sys.version_info.major == 2):
    from plotUtils import makedirs_ok
else:
    def makedirs_ok(*args):
        makedirs(*args, exist_ok=True)


_path_base = None  # Module-wide variable
_outdir_data = "data"
_outdir_plot = path.join("Plot","PhFR")

def getPlots(inputdir, sample, plots, verbose=0):
    fname = path.join(inputdir, sample+".root")
    retrieved = []
    with TFileContext(fname) as rFile:
        for plot in plots:
            h = rFile.Get(plot)
            if(not h):
                if(verbose > 0): print(Colours.Warn('WARN')+': Could not get "%s" from "%s"' % (plot, fname))
                retrieved.append(None)
            else:
                retrieved.append( copy.deepcopy(h) )
                if(verbose > 1):
                    ignore = c_double(0)
                    n = h.IntegralAndError(0, -1, 0, -1, ignore)
                    print('\t\t{:s} - entries: {:6.0f}'.format(plot, n))
    return retrieved


def joinIfNotNone(strings, connection='_'):
    return connection.join([s for s in strings if s is not None])


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
    # hist.GetXaxis().SetRange(0, hist.GetXaxis().GetNbins()+1)

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
    path_base = _path_base if _path_base is not None else get_path_results_base()
    # return "{:s}/{}/{:s}_{:s}".format(get_path_results_base(), year, analyzer, region)
    return path.join(path_base, str(year), '_'.join([analyzer, region]))


def fakeRateABCD(sample_data, sample_prompt, analyzer="VVGammaAnalyzer", year=2016, region="CRLFR", logx=False, logy=False, fixNegBins=False):
    print("Fake Rate ABCD: data     =", sample_data  )
    print("Fake Rate ABCD: promptMC =", sample_prompt)
    path_in = get_path_results(year=year, analyzer=analyzer, region=region)
    outfname = "{outdir}/ABCD_FR_{dataname:s}-{promptname:s}_{year}.root".format(outdir=_outdir_data, dataname=sample_data["name"], promptname=sample_prompt["name"], year=year)
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
        cFR.SaveAs( "{outdir}/ABCD_FR_{dataname:s}-{promptname:s}_{year}.{ext:s}".format(outdir=_outdir_plot, dataname=sample_data["name"], promptname=sample_prompt["name"], year=year, ext=ext) )
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
        cES.SaveAs( "{outdir}/ABCD_Estimate_{dataname:s}-{promptname:s}_{year}.{ext:s}".format(outdir=_outdir_plot, dataname=sample_data["name"], promptname=sample_prompt["name"], year=year, ext=ext) )
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
    cKF.SaveAs( '{outdir}/ABCD_kFactor_{dataname:s}-{promptname:s}_{year}.{ext:s}'.format(outdir=_outdir_plot, dataname=sample_data['name'], promptname=sample_prompt['name'], year=year, ext='png') )
    del cKF, to_preserve_KF
    
    return hFR, hES, hKF


    return hFR, hES


def fakeRateABCD_noSubtract(sample, analyzer="VVGammaAnalyzer", year=2016, region="CRLFR", logx=False, logy=False, fixNegBins=False):
    print("Fake Rate ABCD: sample   =", sample  )
    path_in = get_path_results(year=year, analyzer=analyzer, region=region)
    outfname = "{outdir:s}/ABCD_FR_{samplename:s}_{year}.root".format(outdir=_outdir_data, samplename=sample["name"], year=year)
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
        cFR.SaveAs( "{outdir}/ABCD_FR_{samplename:s}_{year}.{ext:s}".format(outdir=_outdir_plot, samplename=sample["name"], year=year, ext=ext) )
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
        cES.SaveAs( "{outdir}/ABCD_Estimate_{samplename:s}_{year}.{ext:s}".format(outdir=_outdir_plot, samplename=sample["name"], year=year, ext=ext) )
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
    cKF.SaveAs( '{outdir}/ABCD_kFactor_{samplename:s}_{year}.{ext:s}'.format(outdir=_outdir_plot, samplename=sample['name'], year=year, ext='png') )
    del cKF, to_preserve_KF
    
    return hFR, hES, hKF


def getPassFailLtoT_noSubtract_regex(sample, method, variable, regex, year, analyzer, region, fixNegBins=False):
    regex_compiled = re.compile(regex)
    path_in = get_path_results(year=year, analyzer=analyzer, region=region)
    channels = findChannelsInFile(path.join(path_in, sample['file']+'.root'), method, variable)

    print('\tDEBUG: channels =', channels)
    print('\tDEBUG: regex    =', regex_compiled.pattern)
    selected = [ c  for c in channels if regex_compiled.search(c) ]
    print('\tDEBUG: selected =', selected)
    assert len(selected) > 0, "No matches for regex:{regex}, method:{method}, variable:{variable} in {fname}".format(regex=regex, method=method, variable=variable, fname=path.join(path_in, sample['file']+'.root'))

    if(method.startswith('LtoT')):
        raise NotImplementedError('Use VLtoL')

    listPASS = []
    listFAIL = []
    for channel in selected:
        if(sample['name'] == 'data'):
            hPASS, hFAIL = getPlots(path_in, sample['file'], [ 'PhFR_%s_%s_%s_data_%s'      % (method, variable, channel, s) for s in ['PASS', 'FAIL'] ] , verbose=1)

            listPASS.append(hPASS)
            listFAIL.append(hFAIL)

            # hTOTAL=copy.deepcopy(hPASS)
            # hTOTAL.Add(hFAIL)
            # ignore = c_double(0)
            # nP = hPASS .IntegralAndError(0, -1, 0, -1, ignore)
            # nT = hTOTAL.IntegralAndError(0, -1, 0, -1, ignore)
            # print('\tPASS: {:6.0f} - TOTAL: {:6.0f} - <FR>: {:.2f}'.format(nP, nT, nP/nT))
            del hPASS, hFAIL
        else:
            hPASSp, hFAILp = getPlots(path_in, sample['file'], [ 'PhFR_%s_%s_%s_prompt_%s'    % (method, variable, channel, s) for s in ['PASS', 'FAIL'] ] , verbose=1)
            hPASSn, hFAILn = getPlots(path_in, sample['file'], [ 'PhFR_%s_%s_%s_nonprompt_%s' % (method, variable, channel, s) for s in ['PASS', 'FAIL'] ] , verbose=1)

            listPASS.append( addIfExisting(hPASSp, hPASSn) )
            listFAIL.append( addIfExisting(hFAILp, hFAILn) )

    hPASS = addIfExisting(*listPASS)
    hFAIL = addIfExisting(*listFAIL)

    assert hPASS and hFAIL

    if(sample.get('fixNegBins') and fixNegBins):
        fix_neg_bins(hPASS)
        fix_neg_bins(hFAIL)

    return hPASS, hFAIL


def getPassFailLtoT_regex(sample_data, sample_prompt, method, variable, regex, year, analyzer, region, fixNegBins=False):
    regex_compiled = re.compile(regex)
    path_in = get_path_results(year=year, analyzer=analyzer, region=region)
    channels = findChannelsInFile(path.join(path_in, sample_data['file']+'.root'), method, variable)

    # print('\tDEBUG: regex    =', regex_compiled.pattern)
    selected = [ c  for c in channels if regex_compiled.search(c) ]
    # print('\tDEBUG: selected =', selected)

    if(method.startswith('LtoT')):
        raise NotImplementedError('Use VLtoL')

    listDataPASS   = []
    listDataFAIL   = []
    listPromptPASS = []
    listPromptFAIL = []
    for channel in selected:
        hdata_PASS  , hdata_FAIL   = getPlots(path_in, sample_data['file']  , [ 'PhFR_%s_%s_%s_data_%s'   % (method, variable, channel, s) for s in ['PASS', 'FAIL'] ] , verbose=1)
        hprompt_PASS, hprompt_FAIL = getPlots(path_in, sample_prompt['file'], [ 'PhFR_%s_%s_%s_prompt_%s' % (method, variable, channel, s) for s in ['PASS', 'FAIL'] ] , verbose=1)

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


def fakeRateLtoT_regex(sample_data, sample_prompt, method, variable, regex, year=2016, analyzer='VVGammaAnalyzer', region='CRLFR', pattern_printable=None, fixNegBins=False, rebin_eta=False, rebin_pt=False, **kwargs):
    if(pattern_printable is None):
        pattern_printable = regex.replace('\\','').replace('"','').replace("'","")

    if(sample_prompt is None):
        print('Fake Rate {}: sample   ='.format(method), sample_data)
        samplename  = sample_data['name']
        sampletitle = sample_data['title']

        hPASS, hFAIL = getPassFailLtoT_noSubtract_regex(sample_data, method=method, variable=variable, regex=regex, year=year, analyzer=analyzer, region=region, fixNegBins=fixNegBins)
    else:
        print('Fake Rate {}_{}: data     ='.format(method, pattern_printable), sample_data  )
        print('Fake Rate {}_{}: promptMC ='.format(method, pattern_printable), sample_prompt)
        samplename  = sample_data['name']  +  '-'  + sample_prompt['name']
        sampletitle = sample_data['title'] + ' - ' + sample_prompt['title']

        hPASS, hFAIL = getPassFailLtoT_regex(sample_data, sample_prompt, method=method, variable=variable, regex=regex, year=year, analyzer=analyzer, region=region, fixNegBins=fixNegBins)

    outname = 'FR_{method}_{variable}_{regex}_{samplename}{region}_{year}'.format(method=method, variable=variable, samplename=samplename, region='' if region=='CRLFR' else '_'+region, year=year, regex=pattern_printable)
    title = 'Photon Fake Rate: {method} {regex} (from {sampletitle} in {region})'.format(method=method, sampletitle=sampletitle, regex=pattern_printable, region=region)

    hTOTAL = hFAIL
    hTOTAL.Add(hPASS)

    ignore = c_double(0)
    nP = hPASS .IntegralAndError(0, -1, 0, -1, ignore)
    nT = hTOTAL.IntegralAndError(0, -1, 0, -1, ignore)
    print('\tPASS: {:6.0f} - TOTAL: {:6.0f} - <FR>: {:.2f}'.format(nP, nT, nP/nT))

    ## Fake rate = PASS/TOTAL ##
    if(rebin_eta):
        y_bins_orig = array('d', hPASS.GetYaxis().GetXbins())
        y_bins      = array('d', (y_bins_orig[0], y_bins_orig[2], y_bins_orig[3], y_bins_orig[5]))
        hPASS       = rebin2D(hPASS , y_bins=y_bins)
        hTOTAL      = rebin2D(hTOTAL, y_bins=y_bins)
    if(rebin_pt):
        x_bins_orig = array('d', hPASS.GetXaxis().GetXbins())
        x_bins      = array('d', (x_bins_orig[0], x_bins_orig[1], x_bins_orig[2], x_bins_orig[4]))
        hPASS       = rebin2D(hPASS , x_bins=x_bins)
        hTOTAL      = rebin2D(hTOTAL, x_bins=x_bins)

    hFR = hPASS
    hFR.Divide(hTOTAL)

    x_name, y_name = variable.split('-')
    hFR.GetXaxis().SetTitle(varname_to_title(x_name))
    hFR.GetYaxis().SetTitle(varname_to_title(y_name))
    plotFR_LtoT(hFR, outname, title, **kwargs)
    return hFR


def getPassFailLtoT(sample_main, sample_subtr, analyzer, year, region, method, variable, fixNegBins):
    if(method.startswith('LtoT')):
        raise NotImplementedException('LtoT is deprecated. Use VLtoL')

    path_in = get_path_results(year=year, analyzer=analyzer, region=region)

    if(sample_main['name'] == 'data'):
        hmain_PASS, hmain_FAIL   = getPlots(path_in, sample_main['file'] , [ 'PhFR_%s_%s_data_%s'    % (method, variable, s) for s in ['PASS', 'FAIL'] ] )
        assert hmain_PASS, 'Could not get the PASS histogram for data'
        assert hmain_FAIL, 'Could not get the FAIL histogram for data'
    else:
        hmain_PASSp, hmain_FAILp = getPlots(path_in, sample_main['file'], [ 'PhFR_%s_%s_prompt_%s'    % (method, variable, s) for s in ['PASS', 'FAIL'] ] )
        hmain_PASSn, hmain_FAILn = getPlots(path_in, sample_main['file'], [ 'PhFR_%s_%s_nonprompt_%s' % (method, variable, s) for s in ['PASS', 'FAIL'] ] )

        hmain_PASS = addIfExisting(hmain_PASSp, hmain_PASSn)
        hmain_FAIL = addIfExisting(hmain_FAILp, hmain_FAILn)
        assert hmain_PASS, 'Could not get the PASS main (data) histogram'
        assert hmain_FAIL, 'Could not get the FAIL main (data) histogram'

    assert (hmain_PASS   and hmain_FAIL  ), "Could't get the 2 main (data) histograms!"

    if(sample_subtr is not None):
        hsubtr_PASS, hsubtr_FAIL = getPlots(path_in, sample_subtr['file'], [ 'PhFR_%s_%s_prompt_%s' % (method, variable, s) for s in ['PASS', 'FAIL'] ] )
        assert (hsubtr_PASS and hsubtr_FAIL), "Could't get the 2 prompt MC histograms!"

    if(variable.startswith('pt-dRl')):
        to_rebin = [hmain_PASS, hmain_FAIL]
        if(sample_subtr is not None):
            to_rebin += [hsubtr_PASS, hsubtr_FAIL]
        for h in to_rebin:
            h.RebinY(5) #4

    # if(sample_main.get('fixNegBins', False) and fixNegBins):
    #     for h in [hmain_PASS, hmain_FAIL]:
    #         fix_neg_bins(h)
    # if(sample_subtr is not None and sample_subtr.get('fixNegBins', False) and fixNegBins):
    #     for h in [hsubtr_PASS, hsubtr_FAIL]:
    #         fix_neg_bins(h)

    if(sample_subtr is not None):
        hmain_PASS.Add(hsubtr_PASS, -1)
        hmain_FAIL.Add(hsubtr_FAIL, -1)

    if(fixNegBins):
        for h in [hmain_PASS, hmain_FAIL]:
            fix_neg_bins(h)

    return hmain_PASS, hmain_FAIL


def varname_to_title(name):
    if  (name == 'pt'  ): return 'p_{T} [GeV/c]'
    elif(name == 'aeta'): return '|#eta|'
    if  (name == 'dRl' ): return '#DeltaR(l, #gamma)'

def fakeRateLtoT(sample_data, sample_prompt, analyzer='VVGammaAnalyzer', year=2016, region='CRLFR', method='LtoT', variable='pt-aeta', fixNegBins=False, rebin_eta=False, rebin_pt=False, **kwargs):
    if(sample_prompt is None):
        print('Fake Rate {}: sample   ='.format(method), sample_data)
        samplename  = sample_data['name']
        sampletitle = sample_data['title']
    else:
        print('Fake Rate {}: data     ='.format(method), sample_data  )
        print('Fake Rate {}: promptMC ='.format(method), sample_prompt)
        samplename  = sample_data['name']  +  '-'  + sample_prompt['name']
        sampletitle = sample_data['title'] + ' - ' + sample_prompt['title']

    hPASS, hFAIL = getPassFailLtoT(sample_data, sample_prompt, analyzer=analyzer, year=year, region=region, method=method, variable=variable, fixNegBins=fixNegBins)

    outname = 'FR_{method}_{variable}_{samplename}{region}_{year}'.format(method=method, variable=variable, samplename=samplename, region='' if region=='CRLFR' else '_'+region, year=year)
    if(method.startswith('LtoT')):
        title = 'Photon Fake Rate: 5/(3+4+5) (from {sampletitle} in {region})'.format(sampletitle=sampletitle, region=region)
    else:
        title = 'Photon Fake Rate: {method} (from {sampletitle} in {region})'.format(method=method, sampletitle=sampletitle, region=region)

    hTOTAL = hFAIL
    hTOTAL.Add(hPASS)

    ignore = c_double(0)
    nP = hPASS .IntegralAndError(0, -1, 0, -1, ignore)
    nT = hTOTAL.IntegralAndError(0, -1, 0, -1, ignore)
    print('\tPASS: {:6.0f} - TOTAL: {:6.0f} - <FR>: {:.2f}'.format(nP, nT, nP/nT))

    ## Fake rate = PASS/TOTAL ##
    if(rebin_eta):
        y_bins_orig = array('d', hPASS.GetYaxis().GetXbins())
        y_bins      = array('d', (y_bins_orig[0], y_bins_orig[2], y_bins_orig[3], y_bins_orig[5]))
        hPASS       = rebin2D(hPASS , y_bins=y_bins)
        hTOTAL      = rebin2D(hTOTAL, y_bins=y_bins)
    if(rebin_pt):
        x_bins_orig = array('d', hPASS.GetXaxis().GetXbins())
        x_bins      = array('d', (x_bins_orig[0], x_bins_orig[1], x_bins_orig[2], x_bins_orig[4]))
        hPASS       = rebin2D(hPASS , x_bins=x_bins)
        hTOTAL      = rebin2D(hTOTAL, x_bins=x_bins)

    hFR = hPASS
    hFR.Divide(hTOTAL)

    x_name, y_name = variable.split('-')
    hFR.GetXaxis().SetTitle(varname_to_title(x_name))
    hFR.GetYaxis().SetTitle(varname_to_title(y_name))
    plotFR_LtoT(hFR, outname, title, **kwargs)

    return hFR


def plotFR_LtoT(hFR, outname, title, logx=False, logy=False, do_title=True, range_FR_z=[0.,1.], **kwargs):
    outfname = path.join(_outdir_data, outname+'.root')
    picname  = path.join(_outdir_plot, outname)

    cFR = ROOT.TCanvas( 'cFR_{}'.format(outname), title, 1200, 900 )
    cFR.cd()

    min_draw , max_draw  = range_FR_z
    min_value = hFR.GetBinContent(hFR.GetMinimumBin())
    max_value = hFR.GetBinContent(hFR.GetMaximumBin())
    if(min_value < min_draw): print(Colours.Warn('WARN')+': MIN value (%f) smaller than Z draw limit (%f) for %s' % (min_value, min_draw, outname))
    if(max_value > max_draw): print(Colours.Warn('WARN')+': MAX value (%f) greater than Z draw limit (%f) for %s' % (max_value, max_draw, outname))

    hFR.SetTitle(title if do_title else '')
    hFR.SetName( 'PhFR' )
    hFR.SetMinimum(min_draw)
    hFR.SetMaximum(max_draw)
    
    with TFileContext(outfname, 'UPDATE') as tf:
        tf.cd()
        hFR.Write(hFR.GetName(), ROOT.TObject.kOverwrite)
    print('Output (with prompt MC subtraction) in: "{:s}"'.format(outfname))

    to_preserve_FR = beautify(cFR, hFR, logx, logy)

    for ext in ['png', 'pdf']:
        cFR.SaveAs( picname+'.'+ext )
    del cFR, to_preserve_FR


def plotRatio(h1, h2, name="ratio", title="ratio", do_title=True):
    assert h1, "h1 missing"
    assert h2, "h2 missing"

    ratio = copy.deepcopy(h1)
    ratio.Divide(h2)
    ratio.SetName('PhFRSF')
    ratio.SetTitle(title if do_title else '')

    ratio.SetMaximum(2.)
    ratio.SetMinimum(0.)
    ratio.SetContour(21)  # Set to an odd number so that 1 is in the middle of a bin in the color gradient
    
    c = ROOT.TCanvas("cratio_{:s}".format(name), "ratio", 1200, 900)
    c.cd()
    
    ratio.Draw("colz texte")
    to_preserve = beautify(c, ratio, True, False)
    # ratio.GetYaxis().SetRange(1, ratio.GetXaxis().GetNbins())
    c.SaveAs('{:s}/{:s}.png'.format(_outdir_plot, name))
    with TFileContext(path.join(_outdir_data, name+'.root'), 'RECREATE') as tf:
        tf.cd()
        ratio.Write()
        print('Output (ratio) in "{:s}"'.format(tf.GetName()))
    del to_preserve, c

    return ratio


def plotProfiled(h2, name=None, title='profile', direction='X', do_title=True, **kwargs):
    if(not h2.Class().InheritsFrom("TH2")):
        print(Important('ERROR')+': "{:s}" does not inherit from TH2'.format(h2.GetName()))
        return
    if(name is None):
        name = h2.GetName()

    _style = [
        {"color":ROOT.kRed     , "marker":ROOT.kFullCircle},
        {"color":ROOT.kBlue    , "marker":ROOT.kFullTriangleUp},
        {"color":ROOT.kGreen   , "marker":ROOT.kFullTriangleDown},
        {"color":ROOT.kBlack   , "marker":ROOT.kFullSquare},
        {"color":ROOT.kMagenta , "marker":ROOT.kFullCrossX},
        {"color":ROOT.kCyan    , "marker":ROOT.kFullStar},
        {"color":ROOT.kYellow+3, "marker":ROOT.kFullCross}
    ]

    cprof = ROOT.TCanvas("c_"+name, h2.GetName(), 1200, 900)
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

    if(nbins_proj > len(_style)):
        to_divide = nbins_proj // len(_style)
        print(Colours.Warn("WARN")+": style(%d) insufficient for %d. Rebinning(%d)" %(len(_style), nbins_proj, to_divide))

        if  (direction=='Y'):
            h2.RebinY(to_divide)
            nbins_proj = h2.GetNbinsY()
        elif(direction=='X'):
            h2.RebinX(to_divide)
            nbins_proj = h2.GetNbinsX()

        print("INFO: resulting in:", nbins_proj)

    h1s = []
    draw_edges = array('d')
    for i in range(1, axis_2D_draw.GetNbins()+2):
        draw_edges.append(axis_2D_draw.GetBinLowEdge(i))
    draw_edges.append(axis_2D_draw.GetBinUpEdge(axis_2D_draw.GetNbins()+1))
    # print('edges:',draw_edges)

    zmax, zmin = 0, 0
    legend = ROOT.TLegend(.7,.7,.9,.9)
    
    for j in range(1, axis_2D_proj.GetNbins()+1):
        h1 = ROOT.TH1D("{:s}_bin{:d}".format(h2.GetName(), j), title if do_title else '', nbins_draw+1, draw_edges)#, axis_2D_draw.GetNbins(), axis_2D_draw.GetXbins().GetArray())
        h1.GetXaxis().SetTitle(axis_2D_draw.GetTitle())
        h1.GetYaxis().SetTitle('FR^{#gamma}')
        allzeroes = True
        for i in range(0, h1.GetNbinsX()+2):
            bin_2d = get_bin_2d(h2, i, j)
            content = h2.GetBinContent(bin_2d)
            error = h2.GetBinError  (bin_2d)
            zmax = max(zmax, content+error)
            zmin = min(zmin, content-error)
            # print("\t>>> bin:", i, "bin_2d:", bin_2d, "content: {:.2f} +- {:.2f}".format(content, error))
            h1.SetBinContent(i, content)
            h1.SetBinError  (i, error)
            if(content > 0): allzeroes = False

        varTitle = axis_2D_proj.GetTitle()
        match = re.search('[^[]+', varTitle)
        if(match): varTitle = match.group().strip()
        legendTitle = "{}<{:s}<{}".format(axis_2D_proj.GetBinLowEdge(j), varTitle, axis_2D_proj.GetBinUpEdge(j))

        if(allzeroes):
            del h1
            print('INFO: skipping slice "%s" because all bins are empty' % (legendTitle))
        else:
            h1s.append(h1)
            legend.AddEntry(h1, legendTitle)

    for i,h in enumerate(h1s):
        h.SetLineColor(  _style[i]["color"] )
        h.SetMarkerColor(_style[i]["color"] )
        h.SetMarkerStyle(_style[i]["marker"])
    
    h1s[0].GetYaxis().SetRangeUser(max(0, zmin - (zmax-zmin)/10), zmax + (zmax-zmin)/10)
    h1s[0].Draw("P0 E1")
    for h1 in h1s[1:]:
        h1.Draw("P0 E1 same")
    
    legend.Draw("SAME")

    for ext in ('png','pdf'):
        cprof.SaveAs( "{:s}/{:s}.{:s}".format(_outdir_plot, name, ext) )
    del legend, cprof


def findChannelsInFile(fname, method, variable):
    regex_compiled = re.compile('PhFR_{method}_{variable}_([^_]+)_(data|(non)?prompt)?_(PASS|FAIL)'.format(method=method, variable=variable))

    with TFileContext(fname) as tf:
        keys = [ k.GetName() for k in tf.GetListOfKeys() ]

    matches = []
    for key in keys:
        m = regex_compiled.match(key)
        if(m):
            matches.append(m.group(1))  # group(0) is the whole expression

    if(len(matches) == 0):
        start = "PhFR_{}:".format(method)
        print("\tERROR: no key matching the following regex:", regex_compiled.pattern,
              "\n\t       Keys starting with", start, [k for k in keys if k.startswith(start)])
    return matches


if __name__ == "__main__":
    ROOT.gStyle.SetPaintTextFormat(".2f")

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

    possible_methods = {"VLtoL", "KtoVL", "KtoVLexcl", "90to80", "ABCD"}

    from argparse import ArgumentParser
    parser = ArgumentParser()
    parser.add_argument("-y", "--year"   , type=int, default=2016, choices=[2016, 2017, 2018], help='Default: %(default)s')
    parser.add_argument("-m", "--method" , choices=possible_methods, default="VLtoL", help='Default: %(default)s')
    parser.add_argument("-t", "--variable", default = "pt-aeta", help='Default: %(default)s')
    parser.add_argument("-s", "--final-state", default=None, help='Default: %(default)s')
    parser.add_argument("-i", "--inputdir", default="results", help='Top directory containing input (default: %(default)s)')
    parser.add_argument("-o", "--outputdir", default=None, help='Subdirectory name for output. If None will use the name of the input dir')
    parser.add_argument(      "--channels", action='store_true', help='Divide by channel (e.g. 2e1m, lepton fails, etc.)')
    parser.add_argument(      "--no-mc"  , dest="do_mc"  , action="store_false", help="Skip MC plots" )
    parser.add_argument(      "--no-data", dest="do_data", action="store_false", help="Skip data plots")
    parser.add_argument(      "--no-title",dest="do_title",action="store_false", help="Do not paint the title on the canvas (default: False)")
    parser.add_argument(      "--linearx", dest='logx', action="store_false", default=True )
    parser.add_argument(      "--logx"   , dest='logx', action="store_true" , default=True )
    parser.add_argument(      "--lineary", dest='logy', action="store_false", default=False)
    parser.add_argument(      "--logy"   , dest='logy', action="store_true" , default=False)
    parser.add_argument(      "--range-FR-z"   , type=float, nargs=2, default=[0., 1.], help='Manually set the z range for FR plots    (default: %(default)s)')
    parser.add_argument(      "--range-ratio-z", type=float, nargs=2, default=[0., 2.], help='Manually set the z range for ratio plots (default: %(default)s)')
    parser.add_argument(      "--rebin-eta", action='store_true', help='Rebin the y axis (eta) so that it has only 3 bins: EB, gap, EE')
    parser.add_argument(      "--rebin-pt" , action='store_true', help='Rebin the x axis (pt) so that it has only 3 bins: [20-35], [35,50], [50,120]')
    parser.add_argument(      "--fix-negative"   , dest='fixNegBins', action='store_true' , help='Set bins with negative content to zero (default: %(default)s)')
    parser.add_argument(      "--no-fix-negative", dest='fixNegBins', action='store_false', help='Set bins with negative content to zero')

    args = parser.parse_args()

    # Set paths for output
    args.inputdir = args.inputdir.rstrip('/')
    _path_base = args.inputdir
    if(args.outputdir is not None):
        _outdir_data = path.join(_outdir_data, args.outputdir)
        _outdir_plot = path.join(_outdir_plot, args.outputdir)
    elif(args.inputdir != parser.get_default('inputdir')):
        _outdir_data = path.join(_outdir_data, args.inputdir.split('results_')[-1])
        _outdir_plot = path.join(_outdir_plot, args.inputdir.split('results_')[-1])
    # print(f'{_outdir_data=}\n{_outdir_plot=}')

    # Set up ROOT options and create dirs
    ROOT.gStyle.SetOptStat(0)
    ROOT.gROOT.SetBatch(True)

    makedirs_ok(_outdir_data)
    makedirs_ok(_outdir_plot)

    # ROOT.gStyle.SetPalette(ROOT.kBird)

    # Start FR plots
    if(args.method == "ABCD"):
        if(args.do_data):
            print('##### ABCD data #####')
            hFR_ABCD_data, hES_data, hKF_data= fakeRateABCD(sampleList["data"], sampleList["ZGToLLG"], year=args.year, logx=True, fixNegBins=False)
            print()
        if(args.do_mc):
            print('##### ABCD MC #####')
            hFR_ABCD_MC  , hES_MC  , hKF_MC  = fakeRateABCD_noSubtract(sampleList["Drell-Yan"]       , year=args.year, logx=True, fixNegBins=True)
            print()
        if(args.do_data and args.do_mc):
            print('##### ABCD ratio #####')
            plotRatio(hFR_ABCD_data, hFR_ABCD_MC, name="ABCD_ratio_data-ZG_over_DY_{}".format(args.year), title="Ratio FR(data-Z#gamma)/FR(DY)")
            print()
        exit(0)


    method   = args.method
    varState = joinIfNotNone([args.variable, args.final_state])
    argsdict = vars(args)

    if(args.channels):
        print("########## CHANNELS   method:", args.method, " variable:", args.variable, " final_state:", args.final_state, "##########")
        if(args.variable.startswith('pt-dRl')):
            if(args.do_data):
                hd_e  = fakeRateLtoT(sampleList['data'], None, **dict(argsdict, variable='pt-dRl_a-e-all-a'))
                hd_m  = fakeRateLtoT(sampleList['data'], None, **dict(argsdict, variable='pt-dRl_a-m-all-a'))
                hd_P  = fakeRateLtoT(sampleList['data'], None, **dict(argsdict, variable='pt-dRl_P-a-all-a'))
                hd_F  = fakeRateLtoT(sampleList['data'], None, **dict(argsdict, variable='pt-dRl_F-a-all-a'))

                hd_eP = fakeRateLtoT(sampleList['data'], None, **dict(argsdict, variable='pt-dRl_P-e-all-a', fixNegBins=True))
                hd_eF = fakeRateLtoT(sampleList['data'], None, **dict(argsdict, variable='pt-dRl_F-m-all-a', fixNegBins=True))
                hd_mP = fakeRateLtoT(sampleList['data'], None, **dict(argsdict, variable='pt-dRl_P-m-all-a', fixNegBins=True))
                hd_mF = fakeRateLtoT(sampleList['data'], None, **dict(argsdict, variable='pt-dRl_F-m-all-a', fixNegBins=True))

                plotRatio(hd_e ,hd_m , name="ratio_{}_{}_data_e_over_m_{}".format(  method, args.variable, args.year), title="Ratio FR(l=e)/FR(l=m) in data"  )
                plotRatio(hd_F ,hd_P , name="ratio_{}_{}_data_F_over_P_{}".format(  method, args.variable, args.year), title="Ratio FR(l=F)/FR(l=P) in data"  )
                plotRatio(hd_eF,hd_eP, name="ratio_{}_{}_data_eF_over_eP_{}".format(method, args.variable, args.year), title="Ratio FR(l=eF)/FR(l=eP) in data")
                plotRatio(hd_mF,hd_mP, name="ratio_{}_{}_data_mF_over_mP_{}".format(method, args.variable, args.year), title="Ratio FR(l=mF)/FR(l=mP) in data")
                print()

            if(args.do_data and args.do_mc):
                hdZG_e  = fakeRateLtoT(sampleList['data'], sampleList['ZGToLLG'], **dict(argsdict, variable='pt-dRl_a-e-all-a'))
                hdZG_m  = fakeRateLtoT(sampleList['data'], sampleList['ZGToLLG'], **dict(argsdict, variable='pt-dRl_a-m-all-a'))
                hdZG_P  = fakeRateLtoT(sampleList['data'], sampleList['ZGToLLG'], **dict(argsdict, variable='pt-dRl_P-a-all-a'))
                hdZG_F  = fakeRateLtoT(sampleList['data'], sampleList['ZGToLLG'], **dict(argsdict, variable='pt-dRl_F-a-all-a'))

                hdZG_eP = fakeRateLtoT(sampleList['data'], sampleList['ZGToLLG'], **dict(argsdict, variable='pt-dRl_P-e-all-a', fixNegBins=True))
                hdZG_eF = fakeRateLtoT(sampleList['data'], sampleList['ZGToLLG'], **dict(argsdict, variable='pt-dRl_F-m-all-a', fixNegBins=True))
                hdZG_mP = fakeRateLtoT(sampleList['data'], sampleList['ZGToLLG'], **dict(argsdict, variable='pt-dRl_P-m-all-a', fixNegBins=True))
                hdZG_mF = fakeRateLtoT(sampleList['data'], sampleList['ZGToLLG'], **dict(argsdict, variable='pt-dRl_F-m-all-a', fixNegBins=True))

                plotRatio(hdZG_e ,hdZG_m , name="ratio_{}_{}_data-ZG_e_over_m_{}".format(  method,args.variable,args.year), title="Ratio FR(l=e)/FR(l=m) in data-Z#gamma"  )
                plotRatio(hdZG_F ,hdZG_P , name="ratio_{}_{}_data-ZG_F_over_P_{}".format(  method,args.variable,args.year), title="Ratio FR(l=F)/FR(l=P) in data-Z#gamma"  )
                plotRatio(hdZG_eF,hdZG_eP, name="ratio_{}_{}_data-ZG_eF_over_eP_{}".format(method,args.variable,args.year), title="Ratio FR(l=eF)/FR(l=eP) in data-Z#gamma")
                plotRatio(hdZG_mF,hdZG_mP, name="ratio_{}_{}_data-ZG_mF_over_mP_{}".format(method,args.variable,args.year), title="Ratio FR(l=mF)/FR(l=mP) in data-Z#gamma")
                print()

        elif(args.variable.startswith('pt-aeta')):
            if(args.do_data):
                print("### data ###")
                hd_e  = fakeRateLtoT_regex(sampleList['data'], None, regex='2[me]\+e[PF]-[YN]', pattern_printable='2x+e', **dict(argsdict))
                hd_m  = fakeRateLtoT_regex(sampleList['data'], None, regex='2[me]\+m[PF]-[YN]', pattern_printable='2x+m', **dict(argsdict))
                hd_P  = fakeRateLtoT_regex(sampleList['data'], None, regex='2[me]\+[em]P-[YN]', pattern_printable='2x+P', **dict(argsdict))
                hd_F  = fakeRateLtoT_regex(sampleList['data'], None, regex='2[me]\+[em]F-[YN]', pattern_printable='2x+F', **dict(argsdict))

                hd_eP = fakeRateLtoT_regex(sampleList['data'], None, regex='2[me]\+eP-[YN]'   , pattern_printable='2x+eP', **dict(argsdict, fixNegBins=True))
                hd_eF = fakeRateLtoT_regex(sampleList['data'], None, regex='2[me]\+eF-[YN]'   , pattern_printable='2x+eF', **dict(argsdict, fixNegBins=True))
                hd_mP = fakeRateLtoT_regex(sampleList['data'], None, regex='2[me]\+mP-[YN]'   , pattern_printable='2x+mP', **dict(argsdict, fixNegBins=True))
                hd_mF = fakeRateLtoT_regex(sampleList['data'], None, regex='2[me]\+mF-[YN]'   , pattern_printable='2x+mF', **dict(argsdict, fixNegBins=True))

                hd_2e = fakeRateLtoT_regex(sampleList['data'], None, regex='2e\+[em][PF]-[YN]', pattern_printable='2e+x' , **dict(argsdict, fixNegBins=True))
                hd_2m = fakeRateLtoT_regex(sampleList['data'], None, regex='2m\+[em][PF]-[YN]', pattern_printable='2m+x' , **dict(argsdict, fixNegBins=True))

                plotRatio(hd_e ,hd_m , name="ratio_{}_{}_data_e_over_m_{}".format(  method,args.variable,args.year),title="Ratio FR(2x+e)/FR(2x+m) in data"  )
                plotRatio(hd_F ,hd_P , name="ratio_{}_{}_data_F_over_P_{}".format(  method,args.variable,args.year),title="Ratio FR(2x+F)/FR(2x+P) in data"  )
                plotRatio(hd_eF,hd_eP, name="ratio_{}_{}_data_eF_over_eP_{}".format(method,args.variable,args.year),title="Ratio FR(2x+eF)/FR(2x+eP) in data")
                plotRatio(hd_mF,hd_mP, name="ratio_{}_{}_data_mF_over_mP_{}".format(method,args.variable,args.year),title="Ratio FR(2x+mF)/FR(2x+mP) in data")
                plotRatio(hd_2e,hd_2m, name="ratio_{}_{}_data_2e_over_2m_{}".format(method,args.variable,args.year),title="Ratio FR(2e+x)/FR(2m+x) in data"  )
                print()

            if(args.do_mc):
                print("###  MC  ###")
                he  = fakeRateLtoT_regex(sampleList['ZZTo4l'], None, regex='2[me]\+e[PF]', pattern_printable='2x+e' , **dict(argsdict, fixNegBins=True))
                hm  = fakeRateLtoT_regex(sampleList['ZZTo4l'], None, regex='2[me]\+m[PF]', pattern_printable='2x+m' , **dict(argsdict, fixNegBins=True))
                hP  = fakeRateLtoT_regex(sampleList['ZZTo4l'], None, regex='2[me]\+[em]P', pattern_printable='2x+P' , **dict(argsdict, fixNegBins=True))
                hF  = fakeRateLtoT_regex(sampleList['ZZTo4l'], None, regex='2[me]\+[em]F', pattern_printable='2x+F' , **dict(argsdict, fixNegBins=True))
                heP = fakeRateLtoT_regex(sampleList['ZZTo4l'], None, regex='2[me]\+eP'   , pattern_printable='2x+eP', **dict(argsdict, fixNegBins=True))
                heF = fakeRateLtoT_regex(sampleList['ZZTo4l'], None, regex='2[me]\+eF'   , pattern_printable='2x+eF', **dict(argsdict, fixNegBins=True))
                hmP = fakeRateLtoT_regex(sampleList['ZZTo4l'], None, regex='2[me]\+mP'   , pattern_printable='2x+mP', **dict(argsdict, fixNegBins=True))
                hmF = fakeRateLtoT_regex(sampleList['ZZTo4l'], None, regex='2[me]\+mF'   , pattern_printable='2x+mF', **dict(argsdict, fixNegBins=True))
                h2e = fakeRateLtoT_regex(sampleList['ZZTo4l'], None, regex='2e\+[em][PF]', pattern_printable='2e+x' , **dict(argsdict, fixNegBins=True))
                h2m = fakeRateLtoT_regex(sampleList['ZZTo4l'], None, regex='2m\+[em][PF]', pattern_printable='2m+x' , **dict(argsdict, fixNegBins=True))

                plotRatio(h2e, h2m, name="ratio_{}_ZZ_2e_over_2m_{}".format(method, args.year), title="Ratio FR(2e+x)/FR(2m+x) in ZZTo4l")
                plotRatio(he , hm , name="ratio_{}_ZZ_e_over_m_{}".format(  method, args.year), title="Ratio FR(2x+e)/FR(2x+m) in ZZTo4l")
                plotRatio(hF , hP , name="ratio_{}_ZZ_F_over_P_{}".format(  method, args.year), title="Ratio FR(2x+F)/FR(2x+P) in ZZTo4l")
                plotRatio(heF, heP, name="ratio_{}_ZZ_eF_over_eP_{}".format(method, args.year), title="Ratio FR(2x+eF)/FR(2x+eP) in ZZTo4l")
                plotRatio(hmF, hmP, name="ratio_{}_ZZ_mF_over_mP_{}".format(method, args.year), title="Ratio FR(2x+mF)/FR(2x+mP) in ZZTo4l")
                print()
            
            if(args.do_data and args.do_mc):
                print("### both ###")
                hdZG_e  = fakeRateLtoT_regex(sampleList['data'], sampleList["ZGToLLG"], regex='2[me]\+e[PF]', pattern_printable='2x+e' , **dict(argsdict))
                hdZG_m  = fakeRateLtoT_regex(sampleList['data'], sampleList["ZGToLLG"], regex='2[me]\+m[PF]', pattern_printable='2x+m' , **dict(argsdict))
                hdZG_P  = fakeRateLtoT_regex(sampleList['data'], sampleList["ZGToLLG"], regex='2[me]\+[em]P', pattern_printable='2x+P' , **dict(argsdict))
                hdZG_F  = fakeRateLtoT_regex(sampleList['data'], sampleList["ZGToLLG"], regex='2[me]\+[em]F', pattern_printable='2x+F' , **dict(argsdict))

                hdZG_eP = fakeRateLtoT_regex(sampleList['data'], sampleList["ZGToLLG"], regex='2[me]\+eP'   , pattern_printable='2x+eP', **dict(argsdict, fixNegBins=True))
                hdZG_eF = fakeRateLtoT_regex(sampleList['data'], sampleList["ZGToLLG"], regex='2[me]\+eF'   , pattern_printable='2x+eF', **dict(argsdict, fixNegBins=True))
                hdZG_mP = fakeRateLtoT_regex(sampleList['data'], sampleList["ZGToLLG"], regex='2[me]\+mP'   , pattern_printable='2x+mP', **dict(argsdict, fixNegBins=True))
                hdZG_mF = fakeRateLtoT_regex(sampleList['data'], sampleList["ZGToLLG"], regex='2[me]\+mF'   , pattern_printable='2x+mF', **dict(argsdict, fixNegBins=True))

                hdZG_2e = fakeRateLtoT_regex(sampleList['data'], sampleList["ZGToLLG"], regex='2e\+[em][PF]', pattern_printable='2e+x' , **dict(argsdict, fixNegBins=True))
                hdZG_2m = fakeRateLtoT_regex(sampleList['data'], sampleList["ZGToLLG"], regex='2m\+[em][PF]', pattern_printable='2m+x' , **dict(argsdict, fixNegBins=True))

                plotRatio(hdZG_e ,hdZG_m , name="ratio_{}_{}_data-ZG_e_over_m_{}".format(  method,args.variable,args.year),title="Ratio FR(2x+e)/FR(2x+m) in data-Z#gamma"  )
                plotRatio(hdZG_F ,hdZG_P , name="ratio_{}_{}_data-ZG_F_over_P_{}".format(  method,args.variable,args.year),title="Ratio FR(2x+F)/FR(2x+P) in data-Z#gamma"  )
                plotRatio(hdZG_eF,hdZG_eP, name="ratio_{}_{}_data-ZG_eF_over_eP_{}".format(method,args.variable,args.year),title="Ratio FR(2x+eF)/FR(2x+eP) in data-Z#gamma")
                plotRatio(hdZG_mF,hdZG_mP, name="ratio_{}_{}_data-ZG_mF_over_mP_{}".format(method,args.variable,args.year),title="Ratio FR(2x+mF)/FR(2x+mP) in data-Z#gamma")
                plotRatio(hdZG_2e,hdZG_2m, name="ratio_{}_{}_data-ZG_2e_over_2m_{}".format(method,args.variable,args.year),title="Ratio FR(2e+x)/FR(2m+x) in data-Z#gamma"  )
                print()

    else:  # args.channels is false
        print("########## INCLUSIVE   method:", args.method, " variable:", args.variable, " final_state:", args.final_state, "##########")
        if(args.do_data):
            hFR_data    = fakeRateLtoT(sampleList["data"], None                 , **dict(argsdict, variable=varState))
            plotProfiled(hFR_data   , name='FR_profiledX_{}_{}_{}_data'        .format(method, varState, args.year), title='FR(#gamma) vs #eta' , direction='X', **argsdict)
            plotProfiled(hFR_data   , name='FR_profiledY_{}_{}_{}_data'        .format(method, varState, args.year), title='FR(#gamma) vs p_{T}', direction='Y', **argsdict)
            print()

        if(args.do_mc):
            hFR_DY   = fakeRateLtoT(sampleList["Drell-Yan"]       , None, **dict(argsdict, variable=varState, fixNegBins=True))
            hFR_ZG   = fakeRateLtoT(sampleList["ZGToLLG"]         , None, **dict(argsdict, variable=varState, fixNegBins=True))
            hFR_ZZ   = fakeRateLtoT(sampleList["ZZTo4l"]          , None, **dict(argsdict, variable=varState, fixNegBins=True))
            hFR_ZZ_4P= fakeRateLtoT(sampleList["ZZTo4l"]          , None, **dict(argsdict, variable=varState, fixNegBins=True, region='SR4P'))
            # hFR_gg   = fakeRateLtoT(sampleList["ggTo4l"]          , None, year=args.year, method=method, variable=varState, logx=True, fixNegBins=True)

            plotRatio(hFR_ZZ_4P, hFR_ZZ, name="ratio_{}_{}_ZZ4P_over_ZZCRLFR_{}".format(method, varState, args.year), title="Ratio FR(ZZ_{4P})/FR(ZZ_{CRLFR}) with "+joinIfNotNone([method, args.final_state], " "))
            print()
        if(args.do_data and args.do_mc):
            hFR_data_ZG = fakeRateLtoT(sampleList["data"], sampleList["ZGToLLG"], **dict(argsdict, variable=varState))
            plotProfiled(hFR_data_ZG, name='FR_profiledX_{}_{}_{}_data-ZGToLLG'.format(method, varState, args.year), title='FR(#gamma) vs #eta' , direction='X', **argsdict)
            plotProfiled(hFR_data_ZG, name='FR_profiledY_{}_{}_{}_data-ZGToLLG'.format(method, varState, args.year), title='FR(#gamma) vs p_{T}', direction='Y', **argsdict)

            plotRatio(hFR_data   , hFR_DY, name="ratio_{}_{}_data_over_DY_{}"   .format(method, varState, args.year), title="Ratio FR(data)/FR(DY) with "        +joinIfNotNone([method, args.final_state], " "))
            plotRatio(hFR_data   , hFR_ZZ, name="ratio_{}_{}_data_over_ZZ_{}"   .format(method, varState, args.year), title="Ratio FR(data)/FR(ZZ) with "        +joinIfNotNone([method, args.final_state], " "))
            plotRatio(hFR_data   , hFR_ZZ_4P, name="ratio_{}_{}_data_over_ZZ4P_{}".format(method, varState, args.year), title="Ratio FR(data)/FR(ZZ_{SR4P}) with "+joinIfNotNone([method, args.final_state], " "))
            plotRatio(hFR_data_ZG, hFR_DY, name="ratio_{}_{}_data-ZG_over_DY_{}".format(method, varState, args.year), title="Ratio FR(data-Z#gamma)/FR(DY) with "+joinIfNotNone([method, args.final_state], " "))
            plotRatio(hFR_data_ZG, hFR_ZZ, name="ratio_{}_{}_data-ZG_over_ZZ_{}".format(method, varState, args.year), title="Ratio FR(data-Z#gamma)/FR(ZZ) with "+joinIfNotNone([method, args.final_state], " "))
            plotRatio(hFR_data_ZG, hFR_ZZ_4P, name="ratio_{}_{}_data-ZG_over_ZZ4P_{}".format(method, varState, args.year), title="Ratio FR(data-Z#gamma)/FR(ZZ_{SR4P}) with "+joinIfNotNone([method, args.final_state], " "))
            # plotRatio(hFR_LtoT_data, hFR_LtoT_ZZ_4P, name="ratio_{}_{}_data-ZG_over_ZZ_4P_{}".format(method, varState, args.year), title="Ratio FR(data-Z#gamma)_{CRLFR}/FR(ZZ_{4P}) with "+joinIfNotNone([method, args.final_state], " "))
            # plotRatio(hFR_LtoT_data, hFR_LtoT_gg, name="ratio_{}_{}_data-ZG_over_gg_{}".format(method, varState, args.year), title="Ratio FR(data-Z#gamma)/FR(gg) with "+joinIfNotNone([method, args.final_state], " "))
            print()

