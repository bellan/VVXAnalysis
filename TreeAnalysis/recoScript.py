#! /usr/bin/env python
from optparse import OptionParser
import ROOT
ROOT.PyConfig.IgnoreCommandLineOptions = True
from ROOT import TH1F,TCanvas, TLegend
import sys,ast
import math
import operator
import os


parser = OptionParser(usage="usage: %prog <final state> [options]")

parser.add_option("-D", "--Dir", dest="Dir",
                  default="test",
                  help="Save folder default is test.")

(options, args) = parser.parse_args()

Dir        = options.Dir

Variables = ["Mass","Mjj","Mjj_Central","Z1Mass","Z2Mass","Z1lep0_sip","Z1lep0_iso","Z0lep0_pt","nJets","nJetsOrig","PtJet1","EtaJet1","PtJet2","EtaJet2","Z1pt","Z2pt","Z2z","Z1z","ptJRatio","ptRatio","PtZZ" ,"Deta","Deta_Central","nJets","nJets_central","Dphi"]
#Variables = ["z","Z2z","Z1z","PtZZ"]
VariablesStudy = ["Deta2Jet","Deta3Jet","Deta_noCentral"]

Variables_CR = ["Mass","nJets","PtZZ","PtJet1"]

cmd = "./python/recoPlot.py  -r SR -c All -S -D " + Dir + " -t "
for var in (Variables+VariablesStudy):
    print cmd+var
    os.system(cmd+var)

    print cmd+var+" -s Fs"
    os.system(cmd+var+" -s Fs")
    print cmd+var+" -s Fs -m pow"
    os.system(cmd+var+" -s Fs -m pow")

cmd = "./python/recoPlot.py -c IrrBkg -S -D " + Dir + " -t "
for var in Variables:
    print cmd+var
    os.system(cmd+var)

cmd = "./python/recoPlot.py -c RedBkg -S -D " + Dir + " -t "
for var in Variables:
    print cmd+var
    os.system(cmd+var)


cmd = "./python/recoPlot.py -c CR -S True -D " + Dir + " -t "
for var in Variables:
    print cmd+var+" -r CR2P2F"
    os.system(cmd+var+" -r CR2P2F")
    print cmd+var+" -r CR3P1F"
    os.system(cmd+var+" -r CR3P1F")

# ./python/recoPlot.py  -r SR -c All -S True 
# ./python/recoPlot.py  -r SR -c All -s Fs -S True 

# ./python/recoPlot.py  -r SR -c All -S True -m mad 
# ./python/recoPlot.py  -r SR -c All -s Fs -S True -m mad 

# ./python/recoPlot.py  -r CR2P2F -c CR -S True 
# ./python/recoPlot.py  -r CR3P1F -c CR -S True 

# ./python/recoPlot.py  -r SR -c RedBkg -S True 
# ./python/recoPlot.py  -r SR -c IrrBkg -S True 

# ./python/recoPlot.py  -r SR -c RedBkg -S True -m mad
# ./python/recoPlot.py  -r SR -c IrrBkg -S True -m mad

# ./python/recoPlot.py  -r SR -c RedBkg -S True -t nJets 
# ./python/recoPlot.py  -r SR -c IrrBkg -S True -t nJets

# ./python/recoPlot.py  -r SR -c RedBkg -S True -t nJets -m mad
# ./python/recoPlot.py  -r SR -c IrrBkg -S True -t nJets -m mad
