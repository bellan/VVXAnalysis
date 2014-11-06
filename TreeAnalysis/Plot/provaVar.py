#! /usr/bin/env python

import ROOT,copy
from ROOT import gSystem, TCanvas, TH1,  TPad, gStyle, TLegend
from collections import OrderedDict
from optparse import OptionParser
import sys
import math
import operator

f = ROOT.TFile("../results/ZZjAnalyzer_SR/H.root")
#f2 = ROOT.TFile("../results/ZZjAnalyzer_SR/DoubleEle.root")
#f3 = ROOT.TFile("../results/ZZjAnalyzer_SR/MuEG.root")

#files=[f1,f2,f3]

hsum=ROOT.TH1F()
hsumvar=ROOT.TH1F()

#for f in files:
h=f.Get("ZZTo4e_nCentralJets")
if h == None:
    print f, " has no entries"
  #  continue
#isFirst =1
#if isFirst:
hsum=copy.deepcopy(h)
#isFirst=0
#hsum.Add(h)
    
#for f in files:
hvar=f.Get("ZZTo4e_nCentralJets_FRVar")
if h == None:
    print f, " has no entries"
    #    continue
    #isFirst =1
    #if isFirst:
hsumvar=copy.deepcopy(hvar)
#        isFirst=0
#    hsumvar.Add(hvar)

print "Integral ", hsum.Integral()

print "Var ", hsumvar.Integral() 

CloseVar=input("digit anything you want to end the script \n")
