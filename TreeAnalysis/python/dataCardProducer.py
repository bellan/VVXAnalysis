#! /usr/bin/env python

import CrossInfo
from CrossInfo import*
import collections

def produceDataCard(SigDic,RedDic,IrrDic,DataDic,SystDicUp,SystDicDown):
#    print SigDic,RedDic,IrrDic,DataDic
#    print SystDicUp
    
    space=8

    del SigDic["Total Signal"]
    del IrrDic["Total Irreducible"]
    SigDic["qq_4l"] = SigDic.pop("ZZTo4lamcatnlo")

    out_file = open("datacard.txt","w")    
    out_file.write("imax 3 number of channels\n")
    out_file.write("jmax*  number of backgrounds\n")
#   out_file.write("jmax "+str(len(IrrDic))+" number of backgrounds\n")
#    out_file.write("jmax "+str(len(SigDic))+" number of signals\n")
#    out_file.write("kmax "+str(len(SystDicUp))+" number of nuisance parameters (sources of systematical uncertainties\n")
    out_file.write("kmax* number of nuisance parameters (sources of systematical uncertainties)\n")
#    out_file.write("imax 3 number of channels\nbin ")

    out_file.write("\n------------\n")
    # data
    for i, (key, value) in enumerate(DataDic["data"].items()):   
        out_file.write(key+" ")
    out_file.write("\nobservation ")


    for i, (key, value) in enumerate(DataDic["data"].items()):   
        out_file.write(str(value["yield"])+" ")

    out_file.write("\n------------\n")
    out_file.write("\nbin     ")

    #signal channel
    for i, (key, value) in enumerate(SigDic.items()):   
        for (k,val) in value.items():
            out_file.write(k+(space-len(k))*" ")
    #irreducible channel
    for i, (key, value) in enumerate(IrrDic.items()):   
        for (k,val) in value.items():
            out_file.write(k+(space-len(k))*" ")

    #reducible channel
    for i, (key, value) in enumerate(RedDic.items()):   
        for (k,val) in value.items():
            out_file.write(k+(space-len(k))*" ")

    out_file.write("\nprocess ")

    #signal process

    for i, (key, value) in enumerate(SigDic.items()):   
        for (k,val) in value.items():
            out_file.write(key+(space-len(key))*" ")

    #irreducible bkg process
    for i, (key, value) in enumerate(IrrDic.items()):   
        for (k,val) in value.items():
            out_file.write(key+(space-len(key))*" ")

    #reducible bkg process
    for i, (key, value) in enumerate(RedDic.items()):   
        for (k,val) in value.items():
            out_file.write("red"+(space-len('red'))*" ")

    out_file.write("\nprocess ")

    for i,val in enumerate(SigDic):
        out_file.write((str(-len(SigDic)+i+1)+(" "*(space-2)))*3)

    for i,val in enumerate(IrrDic):
        print i
        out_file.write((str(1+i)+(" "*(space-1)))*3)

    for i,val in enumerate(RedDic):
        out_file.write((str(len(IrrDic)+1+i)+(" "*(space-1)))*3)

    out_file.write("\nrate    ")

    #signal rate
    for i, (key, value) in enumerate(SigDic.items()):   
        for (k,val) in value.items():
            out_file.write(str(val["yield"])+((space-len(val["yield"]))*" "))

    #irreducible bkg rate
    for i, (key, value) in enumerate(IrrDic.items()):   
        for (k,val) in value.items():
            out_file.write(str(val["yield"])+((space-len(val["yield"]))*" "))

    #reducible bkg rate
    for i, (key, value) in enumerate(RedDic.items()):   
        for (k,val) in value.items():
            out_file.write(str(val["yield"])+((space-len(val["yield"]))*" "))


    out_file.write("\n------------\n")
    out_file.write("\nirrbkg_unc lnN ")

    #irreducible bkg syst

    for i, (key, value) in enumerate(SigDic.items()):   
        for (k,val) in value.items():
            out_file.write("-"+(" "*(space-1)))

    for i, (key, value) in enumerate(IrrDic.items()):   
        for (k,val) in value.items():
            Unc= float(val["Err"])/float(val["yield"])+1
            out_file.write(str(Unc)+" ")

    for i, (key, value) in enumerate(RedDic.items()):   
        for (k,val) in value.items():
            out_file.write("- ")
    out_file.write("\nredbkg_unc lnN ")

    #reducible bkg syst

    for i, (key, value) in enumerate(SigDic.items()):   
        for (k,val) in value.items():
            out_file.write("- ")

    for i, (key, value) in enumerate(IrrDic.items()):   
        for (k,val) in value.items():
            out_file.write("- ")

    for i, (key, value) in enumerate(RedDic.items()):   
        for (k,val) in value.items():
            Unc= float(val["Err"])/float(val["yield"])+1
            out_file.write(str(Unc)+" ")

    out_file.write("\n")

    print SystDicUp

    #scale factor
    out_file.write("\nscalefactor lnN ")
    print SystDicUp["Scale Factor"]
    for i, (key, value) in enumerate(SigDic.items()):   
        for (k,val) in value.items():
            print type(SystDicUp["Scale Factor"][k]["yield"])
            out_file.write(str(float(SystDicUp["Scale Factor"][k]["yield"])/100. +1)+((space-len(k))*" "))

    for i, (key, value) in enumerate(IrrDic.items()):   
        for (k,val) in value.items():
            out_file.write("- ")

    for i, (key, value) in enumerate(RedDic.items()):   
        for (k,val) in value.items():
            out_file.write("- ")

    out_file.write("\n")


    #mc choice
    out_file.write("\nmcchoice lnN ")

    for i, (key, value) in enumerate(SigDic.items()):   
        for (k,val) in value.items():
            out_file.write(str(float(SystDicDown["Monte Carlo choice"][k]["yield"])/100. +1)+((space-len(k))*" "))

    for i, (key, value) in enumerate(IrrDic.items()):   
        for (k,val) in value.items():
            out_file.write("- ")

    for i, (key, value) in enumerate(RedDic.items()):   
        for (k,val) in value.items():
            out_file.write("- ")

    out_file.write("\n")




    #Global systematics 
    for glob in GlobSystList: 
        out_file.write(glob["name"]+" lnN ")
        out_file.write(len(SigDic)*3*(str(glob["value"]+1)+" "))
        out_file.write(len(IrrDic)*3*("- "))
        out_file.write(len(RedDic)*3*("- "))
        out_file.write("\n")

    out_file.close()
