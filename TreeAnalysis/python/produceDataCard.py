#! /usr/bin/env python
##################################
## G. Pinna (UNITO) - Jun 2015 ##
##################################

from optparse import OptionParser
import ROOT,copy
from ROOT import gSystem, TCanvas, TH1,  TPad, gStyle, TLegend,TGraphAsymmErrors,Math,TArrayD
import collections 
import CrossInfo
from CrossInfo import*
from Colours import *
import sys,ast,os
import math
import operator
import textwrap
import shutil

parser = OptionParser(usage="usage: %prog [options]")

parser.add_option("-V", "--VBS", dest="doVBS",
                  action="store_true",
                  default=False,
                  help="compute data card for all")

parser.add_option("-o", "--out", dest="outputName",
                  default="datacard",
                  help="File output name. Default is datacard")

parser.add_option("-c", "--copy", dest="doCopy",
                  action="store_true",
                  default=False,
                  help="Copy output file in CombinedLimit directory. Default is False")

parser.add_option("-A", "--analysis", dest="analysis",
                  default="ZZ",
                  help="Choose the analysis. Default is ZZ, other option is VBS")

parser.add_option("-f", "--finalstate", dest="finState",
                  default="4l",
                  help="Choose the final state. Default is 4l")

parser.add_option("-r", "--right", dest="doRight",
                  action="store_true",
                  default=False,
                  help="do zz bkg treatment")


(options, args) = parser.parse_args()

doVBS        = options.doVBS
doCopy       = options.doCopy
outputName   = options.outputName
analysis     = options.analysis
finState     = options.finState
doRight      = options.doRight


inputdir_SR  = "./results/"+analysis+"RecoAnalyzer_SR/"
inputdir_MC  = "./results/"+analysis+"MCAnalyzer_MC/"
inputdir_CR  = "./results/"+analysis+"RecoAnalyzer_CR/"

'''
Gen level part:
Compute fiducial cross section in VBS search region for VBS alone and for all ZZ prdouction.
'''

if doVBS:
    fileMCSig = ROOT.TFile(inputdir_MC+"qq_4l2j.root") 
    fileMCBkg = ROOT.TFile(inputdir_MC+"VBSbkg.root") 

else:
    fileMCSig = ROOT.TFile(inputdir_MC+"sig_pow.root") 

print "MC fiducial cross section"
xsTot = 0;
for fn in ('2e2m','4m','4e'):

    h    = copy.deepcopy(fileMCSig.Get("ZZTo"+fn+"_MassGen_01_fr"))
    if h == None: sys.exit("ZZTo"+fn+"_MassGen_01_fr Does not exist in file "+fileMCSig.GetName())
    print Red(fn)
#    if doAll:  xs=(1000*(h.Integral(0,-1)+hbkg.Integral(0,-1)))/Lumi
    xs=(1000*(h.Integral(0,-1)))/Lumi
    print "{0:.3f} fb".format(xs)
    xsTot+=xs 

print Red('4l')
print "{0:.3f} fb\n".format(xsTot)

'''

Reco level part:
procuce datacard for fiducial region cross-section

'''

TotObs = 0;
TotExp = 0;

SignalEntries = 0;

if doVBS:
    fileSig = ROOT.TFile(inputdir_SR+"qq_4l2j.root") 
    fileIrrBkg = ROOT.TFile(inputdir_SR+"VBSbkg.root") 

else:
    fileSig = ROOT.TFile(inputdir_SR+"sig_pow.root") 
    fileIrrBkg = ROOT.TFile(inputdir_SR+"Irr.root") 

fileRedBkg = ROOT.TFile(inputdir_CR+"data.root") 
fileData = ROOT.TFile(inputdir_SR+"data.root") 

if finState not in ("4l","4m","2e2m","4e"): sys.exit("Chose a final state between 4l.4m,2e2m and 4e")
if finState == "4l":
    SigDic    = {"2e2m":{"yield":0,"variation":{},"unc":0},"4m":{"yield":0,"variation":{},"unc":0},"4e":{"yield":0,"variation":{},"unc":0}}
    ZZBkgDic = {"2e2m":{"yield":0,"variation":{},"unc":0},"4m":{"yield":0,"variation":{},"unc":0},"4e":{"yield":0,"variation":{},"unc":0}}
    irrBkgDic = {"2e2m":{"yield":0,"variation":{},"unc":0},"4m":{"yield":0,"variation":{},"unc":0},"4e":{"yield":0,"variation":{},"unc":0}}
    redBkgDic = {"2e2m":{"yield":0,"variation":{},"unc":0},"4m":{"yield":0,"variation":{},"unc":0},"4e":{"yield":0,"variation":{},"unc":0}}
    DataDic   = {"2e2m":{"yield":0,"variation":{},"unc":0},"4m":{"yield":0,"variation":{},"unc":0},"4e":{"yield":0,"variation":{},"unc":0}}
else: 
    SigDic    = {finState:{"yield":0,"variation":{},"unc":0}}
    irrBkgDic = {finState:{"yield":0,"variation":{},"unc":0}}
    ZZBkgDic  = {finState:{"yield":0,"variation":{},"unc":0}}
    redBkgDic = {finState:{"yield":0,"variation":{},"unc":0}}
    DataDic   = {finState:{"yield":0,"variation":{},"unc":0}}



Type =  ["EleSFSq","MuSFSq","Pu","PDF","As"]

if analysis == "VBS":  Type  += ["JES","JER"]

SystType = ["state","FR"]+Type  

for i, (key, value) in enumerate(SigDic.items()):   

    if doRight: h = copy.deepcopy(fileSig.Get("ZZTo"+key+"_Mass_01_fr"))
    else: h = copy.deepcopy(fileSig.Get("ZZTo"+key+"_Mass_01"))

    hIrrBkg  = copy.deepcopy(fileIrrBkg.Get("ZZTo"+key+"_Mass_01"))
    hRedBkg  = copy.deepcopy(fileRedBkg.Get("ZZTo"+key+"_Mass_01"))
    hZZBkg = copy.deepcopy(fileSig.Get("ZZTo"+key+"_Mass_01_nofr"))

    Err=ROOT.Double(0.)
    value["yield"]=h.IntegralAndError(0,-1,Err)
    SigDic[key]["variation"]["stat"] = Err

    TotExp+=value['yield']
    irrBkgDic[key]["yield"]=hIrrBkg.IntegralAndError(0,-1,Err)
    irrBkgDic[key]["variation"]["stat"] = Err

    redBkgDic[key]["yield"]=hRedBkg.IntegralAndError(0,-1,Err)
    redBkgDic[key]["variation"]["FR"]=redBkgDic[key]["yield"]*0.3 
    if redBkgDic[key]["yield"]<0:  redBkgDic[key]["yield"]=abs(redBkgDic[key]["yield"])

    ZZBkgDic[key]["yield"]=hZZBkg.IntegralAndError(0,-1,Err)
    ZZBkgDic[key]["variation"]["stat"] = Err

    TotExp+=irrBkgDic[key]["yield"]
    TotExp+=redBkgDic[key]["yield"]
    TotExp+=ZZBkgDic[key]["yield"]
   
    for t in Type:
        print t
        if doRight:
            hup   = copy.deepcopy(fileSig.Get("ZZTo"+key+"_Mass_"+t+"Up_01_fr"))
            hdown = copy.deepcopy(fileSig.Get("ZZTo"+key+"_Mass_"+t+"Dn_01_fr"))
        else:
            hup   = copy.deepcopy(fileSig.Get("ZZTo"+key+"_Mass_"+t+"Up_01"))
            hdown = copy.deepcopy(fileSig.Get("ZZTo"+key+"_Mass_"+t+"Dn_01"))

        hupZZbkg   = copy.deepcopy(fileSig.Get("ZZTo"+key+"_Mass_"+t+"Up_01_nofr"))
        hdownZZbkg = copy.deepcopy(fileSig.Get("ZZTo"+key+"_Mass_"+t+"Dn_01_nofr"))

        value["variation"][t]     =abs(hup.Integral(0,-1)-hdown.Integral(0,-1))/2.
        ZZBkgDic[key]["variation"][t]=abs(hupZZbkg.Integral(0,-1)-hdownZZbkg.Integral(0,-1))/2.
        irrBkgDic[key]["variation"][t]=abs(hupZZbkg.Integral(0,-1)-hdownZZbkg.Integral(0,-1))/2.


# irrBkgDic[key]["variation"]["stat"]=(hIrrBkg.Integral(0,-1)-irrBkgDic[key]["yield"])       
#        print h.Integral(0,-1),hup.Integral(0,-1),hdown.Integral(0,-1),value["variation"][t]






for i, (key, value) in enumerate(DataDic.items()):   

    h = copy.deepcopy(fileData.Get("ZZTo"+key+"_Mass_01"))
    if h==None:
        print "No events for",key,"in data"
        value["yield"]=0
    else:
        value["yield"]=h.Integral(0,-1)

    TotObs+=value["yield"] 
for i, (key, value) in enumerate(SigDic.items()):        

    for i,(k,v) in enumerate(value["variation"].items()):
        a=(value["yield"]+v)/value["yield"]
        v=a
       

# Irreducible background

for i, (key, value) in enumerate(irrBkgDic.items()):        
    for i,(k,v) in enumerate(value["variation"].items()):
        v=(value["yield"]+v)/value["yield"]


# Reducible background

for i, (key, value) in enumerate(redBkgDic.items()):        
    for i,(k,v) in enumerate(value["variation"].items()):
        v=(value["yield"]+v)/value["yield"]

# ZZ background

for i, (key, value) in enumerate(ZZBkgDic.items()):        
    for i,(k,v) in enumerate(value["variation"].items()):
        v=(value["yield"]+v)/value["yield"]


out_file = open(outputName+"_"+finState+".txt","w")    

if finState == "4l": out_file.write("imax 3 number of channels\n")
else:                out_file.write("imax 1 number of channels\n")

out_file.write("jmax*  number of backgrounds\n")
out_file.write("kmax* number of nuisance parameters (sources of systematical uncertainties)\n")
out_file.write("\n------------\n")
out_file.write("2e2m   4m   4e")
out_file.write("\nobservation ")

for i, (key, value) in enumerate(DataDic.items()):   
    out_file.write(str(value["yield"])+" ")

out_file.write("\n------------\n")

out_file.write("\nbin  ")
for key, value in SigDic.items():   
    out_file.write(key+" ")
for key, value in irrBkgDic.items():   
    out_file.write(key+" ")
for key, value in redBkgDic.items():   
    out_file.write(key+" ")

if doRight:
    for key, value in ZZBkgDic.items():   
        out_file.write(key+" ")

if doRight:
    if finState == "4l": 
        out_file.write("\nprocess  qq_4l2j  qq_4l2j  qq_4l2j  irrBkg  irrBkg  irrBkg redBkg  redBkg  redBkg  ZZbkg  ZZbkg  ZZbkg \n")
        out_file.write("\nprocess  0  0  0  1  1  1  2  2  2  3  3  3\n")
    else: 
        out_file.write("\nprocess  qq_4l2j irrBkg redBkg  ZZbkg \n")
        out_file.write("\nprocess  0  1  2  3\n")
else:
    if finState == "4l": 
        out_file.write("\nprocess  qq_4l2j  qq_4l2j  qq_4l2j  irrBkg  irrBkg  irrBkg redBkg  redBkg  redBkg\n")
        out_file.write("\nprocess  0  0  0  1  1  1  2  2  2\n")
    else: 
        out_file.write("\nprocess  qq_4l2j irrBkg redBkg \n")
        out_file.write("\nprocess  0  1  2 \n")

out_file.write("\nrate  ")
#signal channel

for key, value in SigDic.items():   
    out_file.write("{0:.3f} ".format(value["yield"]))

for key, value in irrBkgDic.items():   
    out_file.write("{0:.3f} ".format(value["yield"]))

for key, value in redBkgDic.items():   
    out_file.write("{0:.3f} ".format(value["yield"]))

if doRight:
    for key, value in ZZBkgDic.items():   
        out_file.write("{0:.3f} ".format(value["yield"]))
 

out_file.write("\n------------\n")

#systematics

for T in SystType:
    out_file.write("\n"+T+" lnN ")
    for key, value in SigDic.items():   
        if value["variation"].has_key(T): out_file.write("{0:.3f} ".format((value["variation"][T]+value["yield"])/value["yield"]))
        else:   out_file.write("   -   ")

    for key, value in irrBkgDic.items():   
        if value["variation"].has_key(T): out_file.write("{0:.3f} ".format((value["variation"][T]+value["yield"])/value["yield"]))
        else:   out_file.write("   -   ")

    for key, value in redBkgDic.items():   
        if value["variation"].has_key(T):  out_file.write("{0:.3f} ".format((value["variation"][T]+value["yield"])/value["yield"]))
        else:   out_file.write("   -   ")

    if doRight:
        for key, value in ZZBkgDic.items():   
            if value["variation"].has_key(T):  out_file.write("{0:.3f} ".format((value["variation"][T]+value["yield"])/value["yield"]))
            else:   out_file.write("   -   ")

GlobSystList.pop()

for syst in GlobSystList:
    out_file.write("\n"+syst["name"]+" lnN ")
    for key, value in SigDic.items():   
        out_file.write("{0:.3f} ".format(1+syst["value"]))
                       
    for key, value in irrBkgDic.items():   
        out_file.write("{0:.3f} ".format(1+syst["value"]))

    for key, value in redBkgDic.items():   
        out_file.write("  -  ")

    if doRight:
        for key, value in ZZBkgDic.items():   
            out_file.write("  -  ")


out_file.write("\nscale lnN")
for key, value in SigDic.items():   
    out_file.write("  1.01  ")
                       
for key, value in irrBkgDic.items():   
    out_file.write("  1.01  ")

for key, value in redBkgDic.items():   
    out_file.write("  1.01  ")

if doRight:
    for key, value in ZZBkgDic.items():   
        out_file.write("  1.01  ")



out_file.write("\nmet lnN")
for key, value in SigDic.items():   
    out_file.write("  -  ")
                       
for key, value in irrBkgDic.items():   
    out_file.write("  -  ")

for key, value in redBkgDic.items():   
    out_file.write("  1.01  ")

if doRight:
    for key, value in ZZBkgDic.items():   
        out_file.write("  -  ")


print "Tot expected {0:.3f} Tot observed {1:.3f} \n".format(TotExp,TotObs)

if doCopy:
    #shutil.copy2(outputName+".txt", "/afs/cern.ch/user/g/gpinnaan/Work/VVX/ZZ/CMSSW_7_4_15_patch1/src/HiggsAnalysis/CombinedLimit/test/") #it doesn't work
    print "cp "+outputName+"_"+finState+".txt ~/Work/VVX/ZZ/CMSSW_7_4_15_patch1/src/HiggsAnalysis/CombinedLimit/test/"
    print "cd  ~/Work/VVX/ZZ/CMSSW_7_4_15_patch1/src/HiggsAnalysis/CombinedLimit/test/ \ncombine -M MaxLikelihoodFit --forceRecreateNLL",outputName+"_"+finState+".txt" 

