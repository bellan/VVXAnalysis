#! /usr/bin/env python
import os
import sys,ast
from os import walk
from python.Colours import *

What= sys.argv[1]
isUnfold = sys.argv[2]
Dir = sys.argv[3]

isUnfold = ast.literal_eval(isUnfold)

InfoType = ["Mass","nJets","nJets_Central","Mjj","Mjj_Central","Deta","Deta_Central","PtJet1","EtaJet1","PtJet2","EtaJet2","PtZZ"]

CmdDir = " -d "+Dir

# Acceptance
if What=="All" or What=="Acc" and not isUnfold: 
    Command = "./python/Acceptance.py -s "+CmdDir+" -t "
    for t in InfoType:
        print Blue("\n"+Command+t+CmdDir)
        os.system(Command+t)
        print Blue("\n"+Command+t+" -S Pow")
        os.system(Command+t+" -S Pow")

if isUnfold: Unfold = " -u"
else: Unfold = ""

#SetHistos
if What=="All" or What=="Set": 
    Command = "./python/SetHistos.py -t "
    for t in InfoType:
        for f in (" -f",""):
            print Blue("\n"+Command+t+f+Unfold)
            os.system(Command+t+f+Unfold)
            print Blue("\n"+Command+t+" -S Pow"+f+Unfold)
            os.system(Command+t+" -S Pow"+f+Unfold)

#Cross
if What=="All" or What=="Cross": 
    Command = "./python/ComputeCross.py -s "+CmdDir+" -t "
    for t in InfoType:
        for f in (" -f",""):
            for n in (" -n",""):
                print Blue("\n"+Command+t+f+n+Unfold)
                os.system(Command+t+f+n+Unfold)
                print Blue("\n"+Command+t+" -S Pow"+f+n+Unfold)
                os.system(Command+t+" -S Pow"+f+n+Unfold)


else: sys.exit("Choose betwewn, All, Acc, Set, Cross")

    


