#! /usr/bin/env python

import ROOT, copy, sys

from ROOT import TH1F

from readSampleInfo import *
from collections import OrderedDict

analysis = "VVXAnalyzer"
cregion = "baseline"
sample = "ZZ"
plot = sys.argv[1] #"ZZTo4l_nCentralJets"

results = "results/{0:s}_{1:s}/".format(analysis,cregion)

typeofsamples = OrderedDict([('WWZJets', ROOT.kOrange+2),
                             ('WZ',ROOT.kOrange+1),
                             ('tt',ROOT.kOrange), 
                             #('Z',ROOT.kYellow),
                             ('H',ROOT.kAzure-6), 
                             ('ZZZJets',ROOT.kAzure-7), 
                             ('WZZJets',ROOT.kAzure-8),
                             ('ZZ',ROOT.kAzure-9)])
#typeofsamples = OrderedDict([('triboson',ROOT.kAzure-7) , ('H',ROOT.kRed) , ('WZ',ROOT.kOrange) , ('WZZ',ROOT.kAzure-8) , ('tt',ROOT.kYellow) , ('Z',ROOT.kGreen-5) , ('ZZ',ROOT.kAzure-9)])
#typeofsamples = OrderedDict([('ZZZJets',ROOT.kAzure-7) , ('WZZJetsSignal',ROOT.kAzure-8) ,  ('ZZTo2e2tau',ROOT.kAzure-9)])

#typeofsamples = OrderedDict([('ZZZJets',ROOT.kAzure) , ('WZZJetsSignal',ROOT.kAzure-1),
#                             ('ZZTo4tau',ROOT.kAzure-2), 
#                             ('ZZTo2mu2tau',ROOT.kAzure-3),  
#                             ('ZZTo2e2tau',ROOT.kAzure-4), 
#                             ('ggZZ4l',ROOT.kAzure-5),
#                             ('ggZZ2l2l',ROOT.kAzure-6),
#                             ('ZZTo4e',ROOT.kAzure-7),
#                             ('ZZTo4mu',ROOT.kAzure-8), 
#                             ('ZZTo2e2mu',ROOT.kAzure-9)])


print typeofsamples

files = {}
for sample in typeofsamples:
    files[sample] = ROOT.TFile(results+sample+".root")
    
totalMC = 0

stack = ROOT.THStack("stack",plot+"_stack")
for sample in typeofsamples:
    h = files[sample].Get(plot)
    print sample, ", total number of events: ", h.Integral()
    totalMC += h.Integral()
    h.SetLineColor(typeofsamples[sample])
    h.SetFillColor(typeofsamples[sample])
    h.SetMarkerStyle(21)
    h.SetMarkerColor(typeofsamples[sample])
    stack.Add(h)

c1 = ROOT.TCanvas(plot+"_stack",plot+"_stack")
c1.cd()

stack.Draw('hist')

fdata = ROOT.TFile(results+"data.root")
hdata = fdata.Get(plot)
print "Total MC = {0:.2f}, data = {1:0.1f}".format(totalMC, hdata.Integral())
hdata.SetMarkerColor(ROOT.kBlack)
hdata.SetLineColor(ROOT.kBlack)
hdata.SetMarkerStyle(21)
hdata.Draw("same")
c1.Update()

c1.SaveAs(plot+"_stack.png")
input()

