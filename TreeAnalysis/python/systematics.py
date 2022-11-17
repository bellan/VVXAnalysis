#!/usr/bin/env python2

from __future__ import print_function
from os import path, environ, getcwd
import ROOT
from ctypes import c_double
from copy import deepcopy
from math import log10, ceil
from json import load, dump
from plotUtils import TFileContext, makedirs_ok
from argparse import ArgumentParser

ROOT.gStyle.SetOptStat('0000')
ROOT.gROOT.SetBatch(True)


def getYrange(*graphs, **kwargs):  # <TGraphAsymmErrors>
    extremes = [0.9, 1.1]
    for graph in graphs:
        if(not graph or not graph.GetY()):
            continue
        for y,eyh,eyl in zip(graph.GetY(), graph.GetEYhigh(), graph.GetEYlow()):
            if(kwargs.get('use_errors', False)):
               extremes.append( y + eyh )
               extremes.append( y - eyh )
               extremes.append( y + eyl )
               extremes.append( y - eyl )
            else:
               extremes.append( y )
    return min(extremes), max(extremes)

def plotSystematics(hCentral, hUp, hDn, syst_values, var='[var]', syst='[syst]', sample='[sample]', region='[region]'):  # <TH1>, <TH1>, <TH1>, <dict> (modified), <str>, <str>, <str>, <str>
    formatInfo = dict(var=var, syst=syst, sample=sample, region=region)

    errorCentral = c_double(0.)
    errorUp = c_double(0.)
    errorDn = c_double(0.)
    integralCentral = hCentral.IntegralAndError(0, hCentral.GetNbinsX()+1, errorCentral)
    integralUp = hUp.IntegralAndError(0, hUp.GetNbinsX()+1, errorUp)
    integralDn = hDn.IntegralAndError(0, hDn.GetNbinsX()+1, errorDn)

    upVar = (integralUp-integralCentral)/integralCentral
    dnVar = (integralDn-integralCentral)/integralCentral
    # print('\t{var}_{syst}'.format(**formatInfo), ' Up: {:.1f} %  Dn: {:.1f} %'.format(100*upVar, 100*dnVar))

    ################################### Definition of the schema for the dict ####################################
    syst_values.setdefault(region, {}).setdefault(var, {}).setdefault(sample, {})[syst] = {'up':upVar, 'dn':dnVar}
    ##############################################################################################################

    
    if(_do_plots):
        c = ROOT.TCanvas('c_{region}_{sample}_{var}_{syst}'.format(**formatInfo), '{var}: {syst} ({region})'.format(**formatInfo), 1600, 900)
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

        digitsBeforeDecimal = ceil(log10(errorCentral.value))
        if(digitsBeforeDecimal <= 2):
            theFormat = '{:.%df}' % (2 - digitsBeforeDecimal)
        else:
            theFormat = '{:.3e}'
        legend.AddEntry(hUp     , ('Up   : %s +- %s' %(theFormat,theFormat)).format(integralUp, errorCentral.value))
        legend.AddEntry(hCentral, ('centr: %s +- %s' %(theFormat,theFormat)).format(integralCentral, errorUp.value))
        legend.AddEntry(hDn     , ('Dn   : %s +- %s' %(theFormat,theFormat)).format(integralDn, errorDn.value))

        
        hCentral.GetXaxis().SetLabelSize(0)  # remove x axis tick labels
        hCentral.SetTitle('{var} {syst} --> {:+.1f}/{:+.1f} %'.format(100*upVar, 100*dnVar, **formatInfo))
        hCentral.Draw('hist')
        hUp.Draw('hist same')
        hDn.Draw('hist same')
        legend.Draw('same')
    
        # ratio
        pad_ratio.cd()
        stat_method = 'pois'
        verbose = False #'muoFake' in syst
        if(verbose):
            stat_method = 'pois v'
        if(integralCentral < 0):
            hCe_forRatio = deepcopy(hCentral)
            hUp_forRatio = deepcopy(hUp)
            hDn_forRatio = deepcopy(hDn)
            hCe_forRatio.Scale(-1)
            hUp_forRatio.Scale(-1)
            hDn_forRatio.Scale(-1)
        else:
            hCe_forRatio = hCentral
            hUp_forRatio = hUp
            hDn_forRatio = hDn
            
        hRatioUp = ROOT.TGraphAsymmErrors(hUp_forRatio, hCe_forRatio, stat_method)
        hRatioDn = ROOT.TGraphAsymmErrors(hDn_forRatio, hCe_forRatio, stat_method)
            
        if(verbose):
            for i in range(hCentral.GetNbinsX()):
                nc = hCentral.GetBinContent(i)
                xc = hCentral.GetBinCenter(i)
                nu = hUp.GetBinContent(i)
                nd = hDn.GetBinContent(i)
                if(nc != 0):
                    print('>>> bin:{:d} ({:.0f}), central:{:+.3f}, up:{:+.3f} ({:+.3f}), dn:{:+.3f} ({:+.3f})'.format(i, xc, nc, nu, nu/nc, nd, nd/nc))
                else:
                    print('>>> bin:{:d} ({:.0f}), central:{:+.3f}, up:{:+.3f} (  nan ), dn:{:.3f} (  nan )'.format(i, xc, nc, nu, nd))
    
        hRatioUp.SetTitle(';'+hCentral.GetXaxis().GetTitle())
        hRatioUp.SetMarkerStyle(ROOT.kFullTriangleUp)
        hRatioDn.SetMarkerStyle(ROOT.kFullTriangleDown)
        hRatioUp.SetLineColor(ROOT.kRed )
        hRatioDn.SetLineColor(ROOT.kBlue)
        hRatioUp.SetMarkerColor(ROOT.kRed )
        hRatioDn.SetMarkerColor(ROOT.kBlue)
        for h in (hRatioUp, hRatioDn):
            h.SetMarkerSize(1.6)
            # h.SetLineWidth(2)
    
        Line = ROOT.TLine(hCentral.GetXaxis().GetXmin(), 1, hCentral.GetXaxis().GetXmax(), 1)
        Line.SetLineWidth(1)
        Line.SetLineStyle(ROOT.kDashed)
        
        ymin, ymax = getYrange(hRatioUp, hRatioDn, use_errors=False)
        range_dn = max(1 + (ymin - 1)*1.2, 0.)
        range_up = min(1 + (ymax - 1)*1.2, 2.)
        frame = pad_ratio.DrawFrame(
            hCentral.GetXaxis().GetXmin(),
            range_dn,
            hCentral.GetXaxis().GetXmax(),
            range_up
        )
        frame.GetXaxis().SetLabelSize(0.1)
        frame.GetYaxis().SetLabelSize(0.08)
    
        Line.Draw('same')
        hRatioUp.Draw('P')
        hRatioDn.Draw('P')
        
        c.SaveAs('Plot/SYS/{region}/{sample}_{var}_{syst}.png'.format(**formatInfo))
        del c
    return upVar, dnVar


def doSystematics(tf, var, syst, syst_values, **kwargs):  # <TFile>, <str>, <str>, <dict> (is modified)
    formatInfo = dict(var=var, syst=syst, file=tf.GetName())
    hCentral = tf.Get('SYS_{var}_central'.format(**formatInfo))
    hUp = tf.Get('SYS_{var}_{syst}_Up'  .format(**formatInfo))
    hDn = tf.Get('SYS_{var}_{syst}_Down'.format(**formatInfo))
    if((not hUp) or (not hDn)):
        print('WARN: var={var}, syst={syst} not found in file={file}'.format(**formatInfo))
        return
    
    upVar, dnVar = plotSystematics(hCentral, hUp, hDn, syst_values, var=var, syst=syst, **kwargs)
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

        sample = path.split('/')[-1].split('.')[0]
        region = path.split('/')[-2].split('_')[-1]
        makedirs_ok('Plot/SYS/{}'.format(region))
    
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
                    plotSystematics(hCentral, hmuR2F1, hmuR0p5F1, syst_values, var=var, syst='QCD-muR', sample=sample, region=region)
                    plotSystematics(hCentral, hmuR1F2, hmuR1F0p5, syst_values, var=var, syst='QCD-F'  , sample=sample, region=region)
                
                doSystematics(tf, var, syst, syst_values, sample=sample, region=region)

if __name__ == '__main__':
    parser = ArgumentParser(description='Calculate systematic variations from rootfiles produced by VVGammaAnalyzer')
    parser.add_argument('-p', '--plots', dest='do_plots', action='store_true')
    args = parser.parse_args()
    
    syst_values = {}
    results_folder = 'results/2016/VVGammaAnalyzer_{region}'

    if(environ.get('CMSSW_BASE', False)):
        basepath = path.join(environ['CMSSW_BASE'], 'src', 'VVXAnalysis', 'TreeAnalysis')
    elif('VVXAnalysis' in getcwd()):
        basepath = path.join(getcwd().split('VVXAnalysis')[0], 'VVXAnalysis', 'TreeAnalysis')
    else:
        basepath = '.'
    
    datapath = path.join(basepath, 'data')
    if(not path.isdir(datapath)):
        makedirs_ok(datapath)
    
    sysJSON = path.join(datapath, 'systematics.json')
    try:
        with open(sysJSON, 'r') as f:
            syst_values = load(f)
    except (IOError, OSError, ValueError):
        syst_values = {}

    _do_plots = args.do_plots
    
    doSystOnFile(path.join(results_folder.format(region='SR3P'), 'fake_leptons.root'  ), syst_values)
    doSystOnFile(path.join(results_folder.format(region='SR3P'), 'WZTo3LNu.root'      ), syst_values)
    doSystOnFile(path.join(results_folder.format(region='SR3P'), 'WZGTo3LNuG.root'    ), syst_values)
    doSystOnFile(path.join(results_folder.format(region='SR3P'), 'ZZTo4l.root'        ), syst_values)
    doSystOnFile(path.join(results_folder.format(region='SR3P'), 'ZGToLLG.root'       ), syst_values)
    doSystOnFile(path.join(results_folder.format(region='SR3P'), 'DYJetsToLL_M50.root'), syst_values)
    doSystOnFile(path.join(results_folder.format(region='SR3P'), 'ZZGTo4LG.root'      ), syst_values)
    
    doSystOnFile(path.join(results_folder.format(region='CR3P1F'), 'data.root'  ), syst_values)
    doSystOnFile(path.join(results_folder.format(region='CR2P2F'), 'data.root'  ), syst_values)
    doSystOnFile(path.join(results_folder.format(region='SR4P'), 'fake_leptons.root'  ), syst_values)
    doSystOnFile(path.join(results_folder.format(region='SR4P'), 'WZTo3LNu.root'      ), syst_values)
    doSystOnFile(path.join(results_folder.format(region='SR4P'), 'ZZTo4l.root'        ), syst_values)
    doSystOnFile(path.join(results_folder.format(region='SR4P'), 'ggTo4l.root'        ), syst_values)
    doSystOnFile(path.join(results_folder.format(region='SR4P'), 'ZZGTo4LG.root'      ), syst_values)
    
    with open(sysJSON, 'w') as fout:
        dump(syst_values, fout, indent=2)
    print('INFO: wrote systematics to "{}"'.format(sysJSON))
