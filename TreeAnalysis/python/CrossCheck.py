#!/usr/bin/python
from __future__ import print_function
import sys
from optparse import OptionParser
from operator import itemgetter
from Colours import *

parser = OptionParser(usage="usage: %prog <final state> [options]")


parser.add_option("-t", "--type", dest="Type",
                  default="",
                  help="Type of cross-check, default is ''. Choose between nJets, m4l,   mZ1,   mZ2,  jet1pt, jet2pt and  choose between m3l, mZ1, lpt, pass for ZL control region (option -l ")

parser.add_option("-l", "--zl", dest="doZL",
                  action= "store_true",
                  default=False,
                  help="Do ZL control region")


(options, args) = parser.parse_args()

Type = options.Type
doZL = options.doZL

[['273158', '62', '98191855', 'mme', '144.03', '91.84', '23.94', '0\n']]

Index = 0;
if doZL:
    if  Type == "m3l":    Index = 4
    elif  Type == "mZ1":    Index = 5
    elif  Type == "lpt":    Index = 6
    elif  Type == "pass":   Index = 7
    elif  Type != "":       sys.exit("wrong type, look option in help -h")
else:
    if    Type == "nJets":  Index = 7
    elif  Type == "m4l":    Index = 4
    elif  Type == "mZ1":    Index = 5
    elif  Type == "mZ2":    Index = 6
    elif  Type == "jet1pt": Index = 8
    elif  Type == "jet2pt": Index = 9
    elif  Type == "mjj":    Index = 10
    elif  Type == "weight": Index = 11
    elif  Type != "":       sys.exit("wrong type, look option in help -h")

def channels(listIn):
    iseeee = False
    ismmmm = False
    iseemm = False
    eeee = []
    eemm = []
    mmmm = []
    llll = []
    
    for i,evlist in enumerate(listIn):
        ev = evlist.split(":")
        if  not ev or ev[0]=="" or ev[0]=="#" or ev[0]=="\n": continue
        if    ev[3] == "eeee":  eeee.append(ev)
        elif  ev[3] == "mmmm":  mmmm.append(ev)
        elif  ev[3] == "eemm":  eemm.append(ev)
        llll.append(ev)

    eeee=sorted(eeee, key=itemgetter(2))
    eeee=sorted(eeee, key=itemgetter(1))
    eeee=sorted(eeee, key=itemgetter(0))
    eemm=sorted(eemm, key=itemgetter(2))
    eemm=sorted(eemm, key=itemgetter(1))
    eemm=sorted(eemm, key=itemgetter(0))
    mmmm=sorted(mmmm, key=itemgetter(2))
    mmmm=sorted(mmmm, key=itemgetter(1))
    mmmm=sorted(mmmm, key=itemgetter(0))
    llll=sorted(llll, key=itemgetter(2))
    llll=sorted(llll, key=itemgetter(1))
    llll=sorted(llll, key=itemgetter(0))

    return eeee,eemm,mmmm,llll

def channelsZL(listIn):
    iseee = False
    iseem = False
    ismmm = False
    ismme = False
    eee = []
    eem = []
    mmm = []
    mme = []
    lll = []
    
    for i,evlist in enumerate(listIn):
        ev = evlist.split(":")
        if  not ev or ev[0]=="" or ev[0]=="#" or ev[0]=="\n": continue
        if    ev[3] == "eee":  eee.append(ev)
        elif  ev[3] == "eem":  eem.append(ev)
        elif  ev[3] == "mmm":  mmm.append(ev)
        elif  ev[3] == "mme":  mme.append(ev)
        lll.append(ev)

    eee=sorted(eee, key=itemgetter(2))
    eee=sorted(eee, key=itemgetter(1))
    eee=sorted(eee, key=itemgetter(0))

    eem=sorted(eem, key=itemgetter(2))
    eem=sorted(eem, key=itemgetter(1))
    eem=sorted(eem, key=itemgetter(0))

    mmm=sorted(mmm, key=itemgetter(2))
    mmm=sorted(mmm, key=itemgetter(1))
    mmm=sorted(mmm, key=itemgetter(0))

    mme=sorted(mme, key=itemgetter(2))
    mme=sorted(mme, key=itemgetter(1))
    mme=sorted(mme, key=itemgetter(0))

    lll=sorted(lll, key=itemgetter(2))
    lll=sorted(lll, key=itemgetter(1))
    lll=sorted(lll, key=itemgetter(0))

    return eee,eem,mmm,mme,lll





def Diff(list1,list2,Type=""):


    isIn = lambda x,c:[y for y in c if (y[0]==x[0] and y[1]==x[1] and y[2]==x[2])]
    diff = lambda c,d:[x for x in c if not len(isIn(x,d))]

    if    Type == "":   return diff(list1,list2)
    else:
        diffList = diff(list1,list2) 
        for i,j in enumerate(diffList):
            list1.remove(j)

    isIn = lambda x,c:[y for y in c if (y[0]==x[0] and y[1]==x[1] and y[2]==x[2] and x[Index]==y[Index])]

    return diff(list1,list2)



def CheckList(lines1,lines2,Type):
    EEEE_1,EEMM_1,MMMM_1,LLLL_1 = channels(lines1)
    EEEE_2,EEMM_2,MMMM_2,LLLL_2 = channels(lines2)
    
    EEEE12_diff = Diff(EEEE_1,EEEE_2,Type)       
    EEEE21_diff = Diff(EEEE_2,EEEE_1,Type)
    EEMM12_diff = Diff(EEMM_1,EEMM_2,Type)
    EEMM21_diff = Diff(EEMM_2,EEMM_1,Type)
    MMMM12_diff = Diff(MMMM_1,MMMM_2,Type)
    MMMM21_diff = Diff(MMMM_2,MMMM_1,Type)
    


    if Type=="":
        print("File 1 eeee",len(EEEE_1),"eemm",len(EEMM_1),"mmmm",len(MMMM_1))
        print("File 2 eeee",len(EEEE_2),"eemm",len(EEMM_2),"mmmm",len(MMMM_2))
        print("event in file 1 not in file2")

        for fin, fin_str in zip((EEEE12_diff,EEMM12_diff,MMMM12_diff),("eeee:","eemm:","mmmm:")):
            print(fin_str)
            if Type=="": 
                for i in fin:
                    for j in i: 
                        print(j, end=' ')
            else: 
                for i in fin: 
                    for l,j in enumerate(i):
                        if l==Index: print(Red(j), end=' ')
                        else: print(j, end=' ')

        
        print("event in file 2 not in file 1")
        for fin, fin_str in zip((EEEE21_diff,EEMM21_diff,MMMM21_diff),("eeee:","eemm:","mmmm:")):
            print(fin_str)
            if Type=="":
                for i in fin:
                    for j in i:
                        print(j, end=' ')
            else:
                for i in fin:
                    for l,j in enumerate(i):
                        if l==Index: print(Red(j), end=' ')
                        else: print(j, end=' ')

    else:

        for fin_1, fin_2, fin_str in zip((EEEE12_diff,EEMM12_diff,MMMM12_diff),(EEEE21_diff,EEMM21_diff,MMMM21_diff),("eeee:","eemm:","mmmm:")):
            print(fin_str)


            for i,j in zip(fin_1,fin_2):
                print("file1", end=' ')
                for l,m in enumerate(i):
                    if l==Index: print(Red(m), end=' ')
                    else: print(m, end=' ')
                print("file2", end=' ')
                for l,m in enumerate(j):
                    if l==Index: print(Red(m), end=' ')
                    else: print(m, end=' ')
                print("\n")
        print("N different events eeee",len(EEEE12_diff),"eemm",len(EEMM12_diff),"mmmm",len(MMMM12_diff))




def CheckListZL(lines1,lines2,Type):
    EEE_1,EEM_1,MMM_1,MME_1,LLL_1 = channelsZL(lines1)
    EEE_2,EEM_2,MMM_2,MME_2,LLL_2 = channelsZL(lines2)

    
    EEE12_diff = Diff(EEE_1,EEE_2,Type)       
    EEE21_diff = Diff(EEE_2,EEE_1,Type)
    EEM12_diff = Diff(EEM_1,EEM_2,Type)
    EEM21_diff = Diff(EEM_2,EEM_1,Type)
    MMM12_diff = Diff(MMM_1,MMM_2,Type)
    MMM21_diff = Diff(MMM_2,MMM_1,Type)
    MME12_diff = Diff(MME_1,MME_2,Type)
    MME21_diff = Diff(MME_2,MME_1,Type)


    if Type=="":
        print("File 1 eee",len(EEE_1),"eem",len(EEM_1),"mmm",len(MMM_1),"mme",len(MME_1))
        print("File 2 eee",len(EEE_2),"eem",len(EEM_2),"mmm",len(MMM_2),"mme",len(MME_2))

        print("event in file 1 not in file2")

        for fin, fin_str in zip((EEE12_diff,EEM12_diff,MMM12_diff,MME12_diff),("eee:","eem:","mmm:","mme:")):
            print(fin_str)
            if Type=="": 
                for i in fin:
                    for j in i: 
                        print(j, end=' ')
            else: 
                for i in fin: 
                    for l,j in enumerate(i):
                        if l==Index: print(Red(j), end=' ')
                        else: print(j, end=' ')

        
        print("event in file 2 not in file 1")
        for fin, fin_str in zip((EEE21_diff,EEM21_diff,MMM21_diff,MME21_diff),("eee:","eem:","mmm:","mme:")):
            print(fin_str)
            if Type=="": 
                for i in fin:
                    for j in i: 
                        print(j, end=' ')
            else: 
                for i in fin: 
                    for l,j in enumerate(i):
                        if l==Index: print(Red(j), end=' ')
                        else: print(j, end=' ')


    else:

        for fin_1, fin_2, fin_str in zip((EEE12_diff,EEM12_diff,MMM12_diff,MME12_diff),(EEE21_diff,EEM21_diff,MMM21_diff,MME21_diff),("eee:","eem:","mmm:","mme:")):
            print(fin_str)


            for i,j in zip(fin_1,fin_2):
                print("file1", end=' ')
                for l,m in enumerate(i):
                    if l==Index:
                        if Index == len(i)-1:  print(Red(m), end=' ')
                        else:  print(Red(m), end=' ')
                    else: print(m, end=' ')
                print("file2", end=' ')
                for l,m in enumerate(j):
                    if l==Index: print(Red(m), end=' ')
                    else: print(m, end=' ')
                print("\n")
        print("N different events eee",len(EEE12_diff),"eem",len(EEM12_diff),"mmm",len(MMM12_diff),"mme",len(MME12_diff))


filename1 = sys.argv[1]
filename2 = sys.argv[2]


try:
    f = open(filename1, "r")
    g = open(filename2, "r")
    
    try:
        lines1 = f.readlines()
        lines2 = g.readlines()
        print("file 1",len(lines1),"events  file2",len(lines2))
        print("Difference in events:") 
        if doZL:  CheckListZL(lines1,lines2,"")
        else:     CheckList(lines1,lines2,"")
        if Type!="":
            print("Difference in events for ",Type) 
            if doZL:  CheckListZL(lines1,lines2,Type)
            else:     CheckList(lines1,lines2,Type)
        f.close()
        g.close()
    finally:
        g.close()
        f.close()
except IOError:
    pass


