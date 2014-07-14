#! /usr/bin/env python

import ROOT, copy

from ROOT import TH1F

from readSampleInfo import *
csvfile = '../Producers/python/samples_8TeV.csv'

analysis = "VVXAnalyzer"
cregion = "baseline"
sample = "ZZ"
plot = "nCentralJets"

results = "results/{0:s}_{1:s}/".format(analysis,cregion)

typeofsamples = typeOfSamples(csvfile)

print typeofsamples

stack = ROOT.THStack("stack",plot+"_stack")
for sample in reversed(typeofsamples):
    if not sample == "Phantom_2e2mu2j" and not sample == "ZZJetsTo4L" and not sample == "W" and not sample == "WW" \
    and not (sample == 'DoubleMu' or sample == 'DoubleEle' or sample == 'MuEG'): 
        f = ROOT.TFile(results+sample+".root")
        h = copy.deepcopy(f.Get(plot))
        print sample, h.Integral()
        h.SetLineColor(2)
        stack.Add(h)

c1 = ROOT.TCanvas(plot+"_stack",plot+"_stack")
c1.cd()
stack.Draw()
c1.SaveAs(plot+"_stack.png")
input()

