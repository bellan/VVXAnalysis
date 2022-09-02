#!/usr/bin/python3

from __future__ import print_function
import sys
import ROOT
if '-b' in sys.argv:
    ROOT.gROOT.SetBatch(True)

from plotUtils import * 

def plotSystematics(tf, var, syst):  # TFile, str, str
    hcentral = tf.Get('SYS_{}_central'.format(var))
    hUp = tf.Get('SYS_{}_{}_Up'.format(var, syst))
    hDn = tf.Get('SYS_{}_{}_Dn'.format(var, syst))
    if(not hUp or not hDn):
        print(f'not found: {var=}, {syst=} in', tf.GetName())
        return
    
    c = ROOT.TCanvas(f"c_{var}_{syst}", f"{var}: {syst}", 200, 100, 1200, 900)
    c.cd()
    hcentral.SetTitle(f'{var} {syst}')
    
    hcentral.SetLineColor(ROOT.kBlack)
    hUp     .SetLineColor(ROOT.kRed  )
    hDn     .SetLineColor(ROOT.kBlue )
    
    hcentral.Draw("hist")
    hUp.Draw("hist same")
    hDn.Draw("hist same")

    c.BuildLegend()
    
    c.SaveAs(f'Plots/SYS/{var}_{syst}.png')



with TFileContext("results/2016/VVGammaAnalyzer_SR3P/WZTo3LNu.root", "READ") as tf:
    names = set()
    for key in tf.GetListOfKeys():
        name = key.GetName()
        if(name[:3] == 'SYS'):
            names.add(name)

    variables   = set([n.split('_')[1] for n in names])
    systematics = set([n.split('_')[2] for n in names])

    print(f'{variables = }')
    print(f'{systematics = }')

    for var in variables:
        for syst in systematics:
            if('QCD' in syst or 'central' in syst):
                continue
            
            plotSystematics(tf, var, syst)
