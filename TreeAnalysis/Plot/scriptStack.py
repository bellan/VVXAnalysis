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

#------------------Set the method!---------

hfake = GetFakeRate(region,"data")


#------------------decomment the interesting one: ---------

#(hstack,leg) = GetStackPlot(region, category)            # ----> Normal stack plot with samples selected by "category" arg
#(hstack,leg) = GetStackPlot_fstate(region, category)    # ----> Stack plot with signal divided in the different final states
#(hstack, leg) = GetSignalDefPlot(category)              # ----> Stack Plot with the signal samples division based on the signal definition (at gen level)

#------------------data plot: ---------
#hdata = GetDataPlot(region)
#--------------------------------------

#hstack.Draw("hist")

if region=="SR_compare":
    hfake.Draw("hist")
    leg.AddEntry(hfake,"Red background from FR method, data driven","lep")  # change "data" with MC for the MC driven method

#hstack.GetXaxis().SetTitle("M_{4l}")
#hstack.GetYaxis().SetTitle("Events/10GeV")

#------------------ To be used for the stack plot in the signal region, and to be changed for different plots ------------------
#hstack.SetMaximum(0.9) 

#------------------- To be used only with data plot ---------
#hdata.Draw("same")
#leg.AddEntry(hdata,"Data","lep")
#--------------------------------------

leg.Draw("same")
    
ROOT.gStyle.SetOptStat(0);
    
c2.Update()

#------------------Change the name just using "GetStackPlot_fstate" ------------------

#c2.SaveAs("Stack_"+region+"_"+category+"_FinState"+".png")
#c2.SaveAs("Stack_"+region+"_"+category+"_FinState"+".root")

c2.SaveAs("Stack_"+region+"_"+category+".png")
c2.SaveAs("Stack_"+region+"_"+category+".root")
CloseVar=input("digit anything you want to end the script \n")        
    
