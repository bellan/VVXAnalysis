#! /usr/bin/env python

##################################
## G. Pinna (UNITO) - Jun 2015 ##
##################################

import ROOT,copy
from ROOT import gSystem, TCanvas, TH1,  TPad, gStyle, TLegend,TGraphAsymmErrors,Math,TArrayD,TTree
import collections 
import CrossInfo
from CrossInfo import*
from optparse import OptionParser
import sys,ast,os
import math
import operator
import textwrap
import collections                                                                                    


parser = OptionParser(usage="usage: %prog <final state> [options]")

parser.add_option("-S", "--Set", dest="Set",
                  default="Pow",
                  help="MC samples Set, default is Pow (Powheg) the other one is Amc (Amcatnlo)")

parser.add_option("-A", "--Analysis", dest="Analysis",
                  default="ZZ",
                  help="Analysis, default is  ZZ, others are ZZFull and HZZ")


(options, args) = parser.parse_args()

Set = options.Set
Analysis  = options.Analysis


list2e2m = {"fs":"2e2m","name":"2e2\mu","mu":1.00276  ,"up":+0.0987631 , "down": -0.0915201 ,"stup":+0.0758923,"stdown":-0.0724288 }  

list4m = {"fs":"4m","name":"4\mu","mu":0.963087  ,"up":+0.0973039 , "down": -0.0910539 ,"stup":+0.0925716,"stdown":-0.0864597 }  

list4e = {"fs":"4e","name":"4e","mu":1.0214  ,"up":+0.166146 , "down": -0.147676 ,"stup":+0.126615,"stdown":-0.117369 }  

list4l = {"fs":"4l","name":"4\ell","mu":0.9855   ,"up":+0.0713972  , "down": -0.0683678 ,"stup":+0.0526103,"stdown":-0.0508797 }  

listFin = {"2e2m":list2e2m,"4m":list4m,"4e":list4e,"4l":list4l}

finstate = ("2e2m","4e","4m","4l")

skipFit  =True
if not skipFit:


    for fin in finstate:
        print "combine -M MaxLikelihoodFit --forceRecreateNLL datacard_"+fin+".txt -n " +fin+"\n\n"
        os.system("combine -M MaxLikelihoodFit --forceRecreateNLL datacard_"+fin+".txt -n " +fin)
        print "combine -M MaxLikelihoodFit --forceRecreateNLL datacard_"+fin+".txt -n -S 0 " +fin+"_stat\n\n"
        os.system("combine -M MaxLikelihoodFit --forceRecreateNLL datacard_"+fin+".txt -S 0 -n " +fin+"_stat")



for fin in finstate:

    fStat = ROOT.TFile("mlfit"+fin+"_stat.root")
    tree = fStat.Get("tree_fit_sb")
    tree.GetEntry(0)
    listFin[fin]["stup"]    =  fStat.tree_fit_sb.muHiErr
    listFin[fin]["stdown"]  =  fStat.tree_fit_sb.muLoErr


    fTot = ROOT.TFile("mlfit"+fin+".root")
    tree = fTot.Get("tree_fit_sb")
    tree.GetEntry(0)
    listFin[fin]["mu"]    =  fTot.tree_fit_sb.mu
    listFin[fin]["up"]    =  fTot.tree_fit_sb.muHiErr
    listFin[fin]["down"]  =  fTot.tree_fit_sb.muLoErr
#    listFin[fin]["up"]    =  math.sqrt(pow(fTot.tree_fit_sb.muHiErr,2) + pow(listFin[fin]["stup"],2))
 #   listFin[fin]["down"]  =  math.sqrt(pow(fTot.tree_fit_sb.muLoErr,2) + pow(listFin[fin]["stdown"],2))


for gb in GlobSystList:
    if gb["name"]=="Lumi": LumiUnc = gb["value"]

print "Signal strength\n\n\n"

print "\\begin{tabular}{|c|c|}"
print "\\hline Process &  Signal strength \\\\"                                                       
print "\\hline"

for l in (list4m,list4e,list2e2m,list4l):
    if l["fs"]=="4l": print "\\hline"
    print "pp~$\\to\Z\Z \\to {0}$ & $ {1:.3f} {2} {3:.3f} {4} {5} +{6:.3f} {4}{7} {2} -{8:.3f} {4}{5} +{9:.3f} {4}{10} \\pm {11:.3f} {12} $\\\\".format(l["name"],l["mu"],"_{",l["stdown"],"}","^{", l["stup"],"~\mathrm{(stat.)}", math.sqrt(math.pow(l["down"],2) - math.pow(l["stdown"],2)),math.sqrt(math.pow(l["up"],2) - math.pow(l["stup"],2)),"~\mathrm{(syst.)}",LumiUnc,"~\mathrm{(lumi.)}")

print " \\hline \n \\end{tabular} \n\n"


print "\\begin{tabular}{lc}"
print "\\hline Process & Fiducial cross section [fb] \\\\ "
print "\\hline"
for l in (list4m,list4e,list2e2m,list4l):
    if l["fs"]=="4l": print "\\hline"
    print "pp~$\\to\Z\Z \\to {0}$ & $ {1:.2f} {2} {3:.2f} {4} {5} +{6:.2f} {4}{7} {2} -{8:.2f} {4}{5} +{9:.2f} {4}{10} \\pm {11:.2f} {12} $\\\\".format(l["name"],l["mu"]*xs_tight[l["fs"]],"_{",l["stdown"]*xs_tight[l["fs"]],"}","^{", l["stup"]*xs_tight[l["fs"]],"~\mathrm{(stat.)}", math.sqrt(math.pow(l["down"],2) - math.pow(l["stdown"],2))*xs_tight[l["fs"]],math.sqrt(math.pow(l["up"],2) - math.pow(l["stup"],2))*xs_tight[l["fs"]],"~\mathrm{(syst.)}",LumiUnc*xs_tight[l["fs"]] ,"~\mathrm{(lumi.)}")

print " \\hline \n \\end{tabular} \n"


BR_4m   = BRmu*BRmu
BR_2e2m = 2*BRmu*BRele
BR_4e   = BRele*BRele    
BR      = BR_4m+BR_2e2m+BR_4e


AccFile     = ROOT.TFile("./Acceptance/Acceptance_Mad_Mass.root")  
Acc     = AccFile.Get("TotAcc4l_Acc").GetVal() 
print "Acc",Acc,"BR",BR
#Acc = 0.5388

#print "\sigma_{pp~$\\to\Z\Z \\to {0}} {1:.2f} {2} {3:.2f} {4} {5} +{6:.2f} {4}{7} {2} -{8:.2f} {4}{5} +{9:.2f} {4}{10} \\pm {11:.2f} {12} $".format(list4l["name"],list4l["mu"]*xs_tight[list4l["fs"]]/(1000*BR*Acc),"_{",list4l["stdown"]*xs_tight[list4l["fs"]]/(1000*BR*Acc),"}","^{", list4l["stup"]*xs_tight[list4l["fs"]]/(1000*BR*Acc),"~\mathrm{(stat.)}", math.sqrt(math.pow(list4l["down"],2) - math.pow(list4l["stdown"],2))*xs_tight[list4l["fs"]]/(1000*BR*Acc),math.sqrt(math.pow(list4l["up"],2) - math.pow(list4l["stup"],2))*xs_tight[list4l["fs"]]/(1000*BR*Acc),"~\mathrm{(syst.)}",0.062*xs_tight[list4l["fs"]]/(1000*BR*Acc) ,"~\mathrm{(lumi.)}")

print "$\sigma {2} pp~\\to\Z\Z \\to {0} {4} {1:.2f} {2} {3:.2f} {4} {5} +{6:.2f} {4}{7} {2} -{8:.2f} {4}{5} +{9:.2f} {4}{10} \\pm {11:.2f} {12} $".format(list4l["name"],list4l["mu"]*xs_tight[list4l["fs"]]/(1000*BR*Acc),"_{",list4l["stdown"]*xs_tight[list4l["fs"]]/(1000*BR*Acc),"}","^{", list4l["stup"]*xs_tight[list4l["fs"]]/(1000*BR*Acc),"~\mathrm{(stat.)}", math.sqrt(math.pow(list4l["down"],2) - math.pow(list4l["stdown"],2))*xs_tight[list4l["fs"]]/(1000*BR*Acc),math.sqrt(math.pow(list4l["up"],2) - math.pow(list4l["stup"],2))*xs_tight[list4l["fs"]]/(1000*BR*Acc),"~\mathrm{(syst.)}",LumiUnc*xs_tight[list4l["fs"]]/(1000*BR*Acc) ,"~\mathrm{(lumi.)}")

#print (list4l["mu"]*xs_tight[list4l["fs"]])/(1000*BR*Acc),"\pm",(list4l["stup"]*xs_tight[list4l["fs"]])/(1000*BR*Acc),(list4l["stdown"]*xs_tight[list4l["fs"]])/(1000*BR*Acc)



