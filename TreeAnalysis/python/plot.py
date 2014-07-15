#! /usr/bin/env python

import ROOT, copy

from ROOT import TH1F

from readSampleInfo import *
from collections import OrderedDict

analysis = "VVXAnalyzer"
cregion = "baseline"
sample = "ZZ"
plot = "nCentralJets"

results = "results/{0:s}_{1:s}/".format(analysis,cregion)

typeofsamples_uno = {'triboson':2,'H':3,'WZ':4,'WZZ':5,'tt':6,'Z':5,'ZZ':1}
typeofsamples = OrderedDict([('triboson',2) , ('H',3) , ('WZ',4) , ('WZZ',5) , ('tt',6) , ('Z',5) , ('ZZ',1)])

print typeofsamples

files = {}
for sample in typeofsamples:
    files[sample] = ROOT.TFile(results+sample+".root")
    

stack = ROOT.THStack("stack",plot+"_stack")
for sample in typeofsamples:
    h = files[sample].Get(plot)
    print sample, ", total number of events: ", h.Integral()
    h.SetLineColor(typeofsamples[sample])
    h.SetFillColor(typeofsamples[sample])
    h.SetMarkerStyle(21)
    h.SetMarkerColor(typeofsamples[sample])
    stack.Add(h)

c1 = ROOT.TCanvas(plot+"_stack",plot+"_stack")
c1.cd()
stack.Draw('hist')
c1.SaveAs(plot+"_stack.png")
input()

