#! /usr/bin/env python

import ROOT,copy
from ROOT import gSystem, TCanvas, TH1,  TPad, gStyle, TLegend
from collections import OrderedDict
from optparse import OptionParser
from Stack import GetFakeRate, GetStackPlot, GetDataPlot, GetStackPlot_fstate, GetSignalDefPlot
import sys
import math
import operator


region = sys.argv[1]
category = sys.argv[2]
#method = sys.argv[3]

c2 = TCanvas( 'c2', 'c2', 200, 10, 900, 700 )

hfake = GetFakeRate(region,"data") ## Set the method!
(hstack,leg) = GetStackPlot(region, category)      #decomment the interesting one
#(hstack,leg) = GetStackPlot_fstate(region, category) #decomment the interesting one
#(hstack, leg) = GetSignalDefPlot(category) #decomment the interesting one
#hdata = GetDataPlot(region)

hstack.Draw("hist")

if region=="SR_compare":
    hfake.Draw("same")
    leg.AddEntry(hfake,"Red background from FR method, data driven","lep")

hstack.GetXaxis().SetTitle("M_{4l}")
hstack.GetYaxis().SetTitle("Events/10GeV")
#hstack.SetMaximum(40)


#hdata.Draw("same")
#hdata.GetYaxis().SetRangeUser(0.,40.)
#leg.AddEntry(hdata,"Data","lep")
leg.Draw("same")
    
ROOT.gStyle.SetOptStat(0);
    
c2.Update()
#c2.SaveAs("Stack_"+region+"_"+category+"_FinState"+".png")
#c2.SaveAs("Stack_"+region+"_"+category+"_FinState"+".root")
c2.SaveAs("Stack_"+region+"_"+category+".png")
c2.SaveAs("Stack_"+region+"_"+category+".root")
CloseVar=input("digit anything you want to end the script \n")        
    
