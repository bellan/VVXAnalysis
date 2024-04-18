#! /usr/bin/env python

from __future__ import print_function
from optparse import OptionParser
import ROOT,copy
from ROOT import gSystem, TCanvas, TH1,  TPad, gStyle, TLegend
import sys,ast
from array import*
import math
import operator
import CrossInfo
from CrossInfo import*
#Type = sys.argv[1]
#StringFileIn  = sys.argv[2]


for var in VarList:
    print(var)
    FileOut = ROOT.TFile("macros/UnfoldingMacros/"+var+"_test/DataToUnfoldFake.root","recreate") 
    FileIn = ROOT.TFile("macros/UnfoldingMacros/"+var+"_test/DataToUnfold.root") 
    for i in ["2e2m","4e","4m"]:
        
        FileIn.cd()        
        hGen = FileIn.Get("DataminusBkg_"+var+"_ZZTo"+i)
        
        Integral = (hGen.Integral())
        value  = Integral/hGen.GetNbinsX()
        # print "Int",Integral[0],"value",value
        hFake = copy.deepcopy(hGen) 
        hFake.SetName("DataminusBkg_"+var+"_ZZTo"+i)   
        for b in range(1,hFake.GetNbinsX()+1):
            hFake.SetBinContent(b,value)
        FileOut.cd()
        hFake.Write()
