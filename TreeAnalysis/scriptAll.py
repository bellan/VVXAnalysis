#! /usr/bin/env python
import os
import sys,ast
from os import walk
from python.Colours import *

What= sys.argv[1]
isUnfold = sys.argv[2]

isUnfold = ast.literal_eval(isUnfold)

InfoType = ["Mass","Jets","Jets_Central","Mjj","Mjj_Central","Deta","Deta_Central","PtJet1","EtaJet1","PtJet2","EtaJet2"]

# Acceptance
if What=="All" or What=="Acc" and not isUnfold: 
    Command = "./python/Acceptance.py -s -t "
    for t in InfoType:
        print Blue("\n"+Command+t)
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
            os.system(Command+t+f)
            print Blue("\n"+Command+t+" -S Pow"+f+Unfold)
            os.system(Command+t+" -S Pow"+f+Unfold)

#Cross
if What=="All" or What=="Cross": 
    Command = "./python/ComputeCross.py -s -t "
    for t in InfoType:
        for f in (" -f",""):
            for n in (" -n",""):
                print Blue("\n"+Command+t+f+n+Unfold)
                os.system(Command+t+f+n+Unfold)
                print Blue("\n"+Command+t+" -S Pow"+f+n+Unfold)
                os.system(Command+t+" -S Pow"+f+n+Unfold)




# ./python/Acceptance.py  -s -t Jets
# ./python/Acceptance.py  -s -t Jets  -S Pow

# ./python/Acceptance.py  -s -t Mjj_Central
# ./python/Acceptance.py  -s -t Mjj_Central  -S Pow

# ./python/Acceptance.py  -s -t Deta
# ./python/Acceptance.py  -s -t Deta  -S Pow

# ./python/Acceptance.py  -s -t Deta_Central
# ./python/Acceptance.py  -s -t Deta_Central  -S Pow


# #"Mass","Jets","Jets_Central","PtJet1","EtaJet1","PtJet2","EtaJet2","Deta","Mjj","Deta_Central","Mjj_Central","dRZZ", "PtZZ","DphiZZ"                                                

# ./python/Acceptance.py  -s -t PtJet1  -S Pow
# ./python/Acceptance.py  -s -t PtJet1

# ./python/Acceptance.py  -s -t EtaJet1  -S Pow
# ./python/Acceptance.py  -s -t EtaJet1
