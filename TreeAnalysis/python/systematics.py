#!/usr/bin/python3

from __future__ import print_function
import sys
import ROOT
from ctypes import c_double
from math import log10, ceil

ROOT.gStyle.SetOptStat('0000')
if '-b' in sys.argv:
    ROOT.gROOT.SetBatch(True)

# from plotUtils import TFileContext  # plotUtils is in python2
class TFileContext(object):
    def __init__(self, *args):
        # print('>>>Opening with args:', args)
        self.tfile = ROOT.TFile(*args)
    def __enter__(self):
        return self.tfile
    def __exit__(self, type, value, traceback):
        # print('<<<Closing TFile '%s'' % (self.tfile.GetName()))
        self.tfile.Close()



def plotSystematics(hCentral, hUp, hDn, var='[var]', syst='[syst]', sample='[sample]'):
    formatInfo = dict(var=var, syst=syst)
    c = ROOT.TCanvas('c_{var}_{syst}'.format(**formatInfo), '{var}: {syst}'.format(**formatInfo), 1600, 900)
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
    print('\t{var}_{syst}'.format(**formatInfo), ' Up: {:.1f} %  Dn: {:.1f} %'.format(100*upVar, 100*dnVar))

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
    
    c.SaveAs('Plots/SYS/{sample}_{var}_{syst}.png'.format(**formatInfo, sample=sample))


def doSystematics(tf, var, syst):  # TFile, str, str
    formatInfo = dict(var=var, syst=syst, file=tf.GetName())
    hCentral = tf.Get('SYS_{var}_central'.format(**formatInfo))
    hUp = tf.Get('SYS_{var}_{syst}_Up'.format(**formatInfo))
    hDn = tf.Get('SYS_{var}_{syst}_Dn'.format(**formatInfo))
    if(not hUp or not hDn):
        print('Warning: var={var}, syst={syst} not found in file={file}'.format(**formatInfo))
        return
    plotSystematics(hCentral, hUp, hDn, var=var, syst=syst, sample=tf.GetName().split('/')[-1].split('.')[0])


def doSystOnFile(path):
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
                    hmuR0p5F1 = tf.Get('SYS_{}_QCDscale_muR0p5F1'.format(var))
                    hmuR2F1   = tf.Get('SYS_{}_QCDscale_muR2F1'  .format(var))
                    hmuR1F0p5 = tf.Get('SYS_{}_QCDscale_muR1F0p5'.format(var))
                    hmuR1F2   = tf.Get('SYS_{}_QCDscale_muR1F2'  .format(var))
                    plotSystematics(hCentral, hmuR2F1, hmuR0p5F1, var=var, syst='QCD-muR', sample=path.split('/')[-1].split('.')[0])
                    plotSystematics(hCentral, hmuR1F2, hmuR1F0p5, var=var, syst='QCD-F'  , sample=path.split('/')[-1].split('.')[0])
                
                doSystematics(tf, var, syst)

results_folder = 'results/2016/VVGammaAnalyzer_SR3P'
doSystOnFile(results_folder+'/WZTo3LNu.root')
doSystOnFile(results_folder+'/ZZTo4l.root')
