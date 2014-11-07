#! /usr/bin/env python

import ROOT,copy
from ROOT import gSystem, TCanvas, TH1,  TPad, gStyle, TLegend
from collections import OrderedDict
from optparse import OptionParser
from Stack import GetFakeRate
import sys
import math
import operator


region = sys.argv[1]


c2 = TCanvas( 'c2', 'c2', 200, 10, 900, 700 )

leg = TLegend(0.61,0.56,0.85,0.81)
leg.SetBorderSize(0)
leg.SetTextSize(0.025)


hfake = GetFakeRate(region,"MC")
hfakedata = GetFakeRate(region,"data")

hfake.SetFillColor(ROOT.kGray)
hfake.SetLineColor(ROOT.kBlack)
hfakedata.SetLineColor(ROOT.kRed)


if region=="SR_compare":
    hfake.Draw("hist")
    hfakedata.Draw("same")
    leg.AddEntry(hfake,"Red background from FR method, MC driven","f") # change "data" with MC for the MC driven method
    leg.AddEntry(hfakedata,"Red background from FR method, data driven","lep")

hfake.GetXaxis().SetTitle("M_{4l}")
hfake.GetYaxis().SetTitle("Events/10GeV")

leg.Draw("same")
    
ROOT.gStyle.SetOptStat(0);
    
c2.Update()


c2.SaveAs("Stack_"+region+"_"+"MCvsDATA"+".png")
c2.SaveAs("Stack_"+region+"_"+"MCvsDATA"+".root")
CloseVar=input("digit anything you want to end the script \n")        
    
