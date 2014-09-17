#! /usr/bin/env python

import ROOT, copy, sys

from ROOT import TH1F

from readSampleInfo import *
from collections import OrderedDict

##### Extract MC plot #####
def getMCPlot(inputdir, category, plot):
    signal = OrderedDict([('ZZZJets',ROOT.kAzure-4) , 
                          ('WZZJets',ROOT.kAzure-5),
                          ('H',ROOT.kAzure-6),
                          ('qq_4l2j',ROOT.kAzure-7),
                          ('gg_4l',ROOT.kAzure-8),
                          ('qq_ZZ',ROOT.kAzure-9)])

    background = OrderedDict([('WWZJets', ROOT.kOrange+2),
                              ('WZ',ROOT.kOrange+1),
                              ('tt',ROOT.kOrange), 
                              ('Z',ROOT.kYellow)])

    typeofsamples = []
    if category == 'all':
        allsamples = background
        allsamples.update(signal)
        typeofsamples = allsamples
    elif category == 'signal':
        typeofsamples = signal
    elif category == 'background':
        typeofsamples = background
    

    files = {}
    for sample in typeofsamples:
        files[sample] = ROOT.TFile(inputdir+sample+".root")
    
    totalMC = 0

    stack = ROOT.THStack("stack",plot+"_stack")
    for sample in typeofsamples:
        h2d = files[sample].Get(plot)
        h = h2d.ProjectionY("_signal",2,-1).Clone()
        if not h: continue
        print sample, ", total number of events: ", h.Integral()
        totalMC += h.Integral()
        h.SetLineColor(typeofsamples[sample])
        h.SetFillColor(typeofsamples[sample])
        h.SetMarkerStyle(21)
        h.SetMarkerColor(typeofsamples[sample])
        stack.Add(h)

    print "Total MC = {0:.2f}".format(totalMC)
    return copy.deepcopy(stack)

######

def getDataPlot(inputdir, plot):
    fdata = ROOT.TFile(inputdir+"data.root")
    hdata = fdata.Get(plot)
    hdata.SetMarkerColor(ROOT.kBlack)
    hdata.SetLineColor(ROOT.kBlack)
    hdata.SetMarkerStyle(21)
    print "Total data = {0:0.1f}".format(hdata.Integral())
    return copy.deepcopy(hdata)


if __name__ == '__main__':

    plot = sys.argv[1] #"ZZTo4l_nCentralJets"

    analysis = "VVXAnalyzer"
    cregion_1 = "SR"
    cregion_2 = "SR"
    type_1 = 'MC'
    type_2 = 'datasubtracted'
    category = 'signal'

    c1 = ROOT.TCanvas(plot,plot)
    c1.cd()

    if type_1 == 'MC':
        hmc = getMCPlot("results/{0:s}_{1:s}/".format(analysis,cregion_1), category, plot)
        hmc.Draw('hist')

    if type_1 == 'data':
        hdata = getDataPlot("results/{0:s}_{1:s}/".format(analysis,cregion_1),plot)
        hdata.Draw()

    if type_2 == 'data':
        hdata = getDataPlot("results/{0:s}_{1:s}/".format(analysis,cregion_2),plot)
        hdata.Draw("same")
        
    if type_2 == 'MC':
        hmc = getMCPlot("results/{0:s}_{1:s}/".format(analysis,cregion_2),category,plot)
        hmc.Draw('histsame')
    
    if type_2 == 'datasubtracted':
        hdata_SR = getDataPlot("results/{0:s}_{1:s}/".format(analysis,'SR'),plot)
        hdata_CR = getDataPlot("results/{0:s}_{1:s}/".format(analysis,'CR'),plot)
        hsub = hdata_SR.Clone("BackgroundSubtracted")
        hsub.Add(hdata_CR, -1)
        hsub.Draw('same')

    c1.Update()

    
    c1.SaveAs(plot+"_stack.png")
    input()
