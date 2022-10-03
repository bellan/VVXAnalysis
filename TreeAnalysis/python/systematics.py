#!/usr/bin/env python2

from __future__ import print_function
from sys import argv
from os import path, environ, getcwd
import ROOT
from ctypes import c_double
from math import log10, ceil
from json import load, dump
from plotUtils import TFileContext, makedirs_ok


ROOT.gStyle.SetOptStat('0000')
if '-b' in argv:
    ROOT.gROOT.SetBatch(True)


def plotSystematics(hCentral, hUp, hDn, syst_values, var='[var]', syst='[syst]', sample='[sample]'):  # <TH1>, <TH1>, <TH1>, <dict> (modified), <str>, <str>, <str>
    formatInfo = dict(var=var, syst=syst, sample=sample)
    c = ROOT.TCanvas('c_{sample}_{var}_{syst}'.format(**formatInfo), '{var}: {syst}'.format(**formatInfo), 1600, 900)
    c.cd()
    pad_histo = ROOT.TPad ('pad_histo', '', 0., 0.34, 1., 1.)
    pad_ratio = ROOT.TPad ('pad_ratio', '', 0., 0.0 , 1., 0.34)
    pad_histo.SetTopMargin   (0.1)
    pad_histo.SetBottomMargin(0.013)
    pad_ratio.SetTopMargin   (0.05)
    pad_ratio.SetBottomMargin(0.20)
    for pad in (pad_histo, pad_ratio):
        pad.SetRightMargin   (0.06)
        pad.SetLeftMargin    (0.10)
        c.cd()
        pad.Draw()
    
    # histograms
    pad_histo.cd()
    legend = ROOT.TLegend(0.70, 0.68, 0.90, 0.88)
    # pad_histo.SetLogy()
    
    hCentral.SetLineWidth(2)
    # hCentral.SetFillStyle(3305)
    # hCentral.SetFillColor(ROOT.kBlack)
    # hUp.SetFillColor(ROOT.kRed)
    # hUp.SetFillStyle(3354)
    # hDn.SetFillColor(ROOT.kBlue)
    # hDn.SetFillStyle(3345)
    hCentral.SetLineColor(ROOT.kBlack)
    hUp     .SetLineColor(ROOT.kRed  )
    hDn     .SetLineColor(ROOT.kBlue )
    hCentral.GetYaxis().SetLabelSize(0.045)

    errStat = c_double(0.)
    integralCentral = hCentral.IntegralAndError(0, hCentral.GetNbinsX()+1, errStat)
    digitsBeforeDecimal = ceil(log10(errStat.value))

    if(digitsBeforeDecimal <= 2):
        theFormat = '{:.%df}' % (2 - digitsBeforeDecimal)
    else:
        theFormat = '{:.3e}'
    
    legend.AddEntry(hCentral, ('centr: %s +- %s' %(theFormat,theFormat)).format(integralCentral, errStat.value))
    integralUp = hUp.IntegralAndError(0, hUp.GetNbinsX()+1, errStat)
    legend.AddEntry(hUp     , ('Up   : %s +- %s' %(theFormat,theFormat)).format(integralUp, errStat.value))
    integralDn = hDn.IntegralAndError(0, hDn.GetNbinsX()+1, errStat)
    legend.AddEntry(hDn     , ('Dn   : %s +- %s' %(theFormat,theFormat)).format(integralDn, errStat.value))

    upVar = (integralUp-integralCentral)/integralCentral
    dnVar = (integralDn-integralCentral)/integralCentral
    # print('\t{var}_{syst}'.format(**formatInfo), ' Up: {:.1f} %  Dn: {:.1f} %'.format(100*upVar, 100*dnVar))
    syst_values.setdefault(sample, {}).setdefault(var, {})[syst] = {'up':1+upVar, 'dn':1+dnVar}

    hCentral.GetXaxis().SetLabelSize(0)  # remove x axis tick labels
    hCentral.SetTitle('{var} {syst} --> {:+.1f}/{:+.1f} %'.format(100*upVar, 100*dnVar, **formatInfo))
    hCentral.Draw('hist')
    hUp.Draw('hist same')
    hDn.Draw('hist same')
    legend.Draw('same')
    
    # ratio
    pad_ratio.cd()
    hRatioUp = ROOT.TGraphAsymmErrors(hUp, hCentral, 'pois')
    hRatioDn = ROOT.TGraphAsymmErrors(hDn, hCentral, 'pois')
    
    hRatioUp.SetTitle(';'+hCentral.GetXaxis().GetTitle())
    hRatioUp.SetMarkerStyle(ROOT.kFullTriangleUp)
    hRatioDn.SetMarkerStyle(ROOT.kFullTriangleDown)
    hRatioUp.SetLineColor(ROOT.kRed )
    hRatioDn.SetLineColor(ROOT.kBlue)
    hRatioUp.SetMarkerColor(ROOT.kRed )
    hRatioDn.SetMarkerColor(ROOT.kBlue)
    
    Line = ROOT.TLine(hCentral.GetXaxis().GetXmin(), 1, hCentral.GetXaxis().GetXmax(), 1)
    Line.SetLineWidth(1)
    Line.SetLineStyle(ROOT.kDashed)

    frame = pad_ratio.DrawFrame(
        hCentral.GetXaxis().GetXmin(),
        0.9,
        hCentral.GetXaxis().GetXmax(),
        1.1)
    frame.GetXaxis().SetLabelSize(0.1)
    frame.GetYaxis().SetLabelSize(0.08)
    
    Line.Draw('same')
    hRatioUp.Draw('P')
    hRatioDn.Draw('P')
    
    c.SaveAs('Plot/SYS/{sample}_{var}_{syst}.png'.format(**formatInfo))
    del c
    return upVar, dnVar


def doSystematics(tf, var, syst, syst_values):  # <TFile>, <str>, <str>, <dict> (is modified)
    formatInfo = dict(var=var, syst=syst, file=tf.GetName())
    hCentral = tf.Get('SYS_{var}_central'.format(**formatInfo))
    hUp = tf.Get('SYS_{var}_{syst}_Up'  .format(**formatInfo))
    hDn = tf.Get('SYS_{var}_{syst}_Down'.format(**formatInfo))
    if((not hUp) or (not hDn)):
        print('WARN: var={var}, syst={syst} not found in file={file}'.format(**formatInfo))
        return
    sample = tf.GetName().split('/')[-1].split('.')[0]
    upVar, dnVar = plotSystematics(hCentral, hUp, hDn, syst_values, var=var, syst=syst, sample=sample)
    # syst_values.setdefault(sample, {}).setdefault(var, {})[syst] = {'up':dnVar, 'dn':dnVar}


def doSystOnFile(path, syst_values):  # <str>, <dict> (to be passed to doSystematics(...) in order to be filled
    with TFileContext(path, 'READ') as tf:
        names = set()
        for key in tf.GetListOfKeys():
            name = key.GetName()
            if(name[:3] == 'SYS'):
                names.add(name)
    
        variables   = set([n.split('_')[1] for n in names])
        systematics = set([n.split('_')[2] for n in names])
    
        print('variables =', variables)
        print('systematics =', systematics)
    
        for var in variables:
            for syst in systematics:
                if('central' in syst):
                    continue
                elif('QCD' in syst):
                    hCentral  = tf.Get('SYS_{}_central'.format(var))
                    hmuR0p5F1 = tf.Get('SYS_{}_QCDscalemuR_Down'.format(var))
                    hmuR2F1   = tf.Get('SYS_{}_QCDscalemuR_Up'  .format(var))
                    hmuR1F0p5 = tf.Get('SYS_{}_QCDscaleF_Down'  .format(var))
                    hmuR1F2   = tf.Get('SYS_{}_QCDscaleF_Up'    .format(var))
                    plotSystematics(hCentral, hmuR2F1, hmuR0p5F1, syst_values, var=var, syst='QCD-muR', sample=path.split('/')[-1].split('.')[0])
                    plotSystematics(hCentral, hmuR1F2, hmuR1F0p5, syst_values, var=var, syst='QCD-F'  , sample=path.split('/')[-1].split('.')[0])
                
                doSystematics(tf, var, syst, syst_values)

if __name__ == '__main__':
    syst_values = {}
    results_folder = 'results/2016/VVGammaAnalyzer_{region}'
    
    doSystOnFile(path.join(results_folder.format(region='SR3P'), 'WZTo3LNu.root'      ), syst_values)
    doSystOnFile(path.join(results_folder.format(region='SR3P'), 'WZGTo3LNuG.root'    ), syst_values)
    # doSystOnFile(path.join(results_folder.format(region='SR3P'), 'ZGToLLG.root'       ), syst_values)
    # doSystOnFile(path.join(results_folder.format(region='SR3P'), 'DYJetsToLL_M50.root'), syst_values)
    
    # doSystOnFile(path.join(results_folder.format(region='SR4P'), 'fake_leptons.root'  ), syst_values)
    # doSystOnFile(path.join(results_folder.format(region='SR4P'), 'WZTo3LNu.root'      ), syst_values)
    doSystOnFile(path.join(results_folder.format(region='SR4P'), 'ZZTo4l.root'        ), syst_values)
    doSystOnFile(path.join(results_folder.format(region='SR4P'), 'ggTo4l.root'        ), syst_values)
    # doSystOnFile(path.join(results_folder.format(region='SR4P'), 'ggTo4e_Contin_MCFM701.root'), syst_values)
    doSystOnFile(path.join(results_folder.format(region='SR4P'), 'ZZGTo4LG.root'      ), syst_values)

    if(environ.get('CMSSW_BASE', False)):
        basepath = path.join(environ['CMSSW_BASE'], 'src', 'VVXAnalysis', 'TreeAnalysis', 'data')
    elif('VVXAnalysis' in getcwd()):
        basepath = path.join(getcwd().split('VVXAnalysis')[0], 'VVXAnalysis', 'TreeAnalysis', 'data')
    else:
        basepath = '.'
    
    if(not path.isdir(basepath)):
        makedirs_ok(basepath)
    
    sysJSON = path.join(basepath, 'systematics.json')
    original = {}
    try:
        with open(sysJSON, 'r') as f:
            original = load(f)
    except (IOError, ValueError):
        pass
    
    original.update(syst_values)
    
    with open(sysJSON, 'w') as fout:
        dump(original, fout, indent=2)
    print('INFO: wrote systematics to "{}"'.format(sysJSON))
