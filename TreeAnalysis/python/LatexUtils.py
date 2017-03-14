#! /usr/bin/env python

import collections

def YieldLatex(Text,SigDic,RedDic,IrrDic,TotDic,DataDic):                                                                                                                                                                                              
    Line = "\\hline"
    Line = ""
    print "\\begin{tabular}{llll}"

    print "\\hline "+Text[0]+" & "+ Text[1] + " & " + Text[2] + " & " + Text[3] + " \\\\ \\hline"   
    for i, (key, value) in enumerate(SigDic.items()):
        print Line, key ," & ", value['2e2m']["yield"] ," & ",value['4m']["yield"], "&",value['4e']["yield"],"\\\\"
    for i, (key, value) in enumerate(RedDic.items()):
        print Line, key ," & ", value['2e2m']["yield"] ," & ",value['4m']["yield"], "&",value['4e']["yield"],"\\\\"
    for i, (key, value) in enumerate(IrrDic.items()):
        print Line, key ," & ", value['2e2m']["yield"] ," & ",value['4m']["yield"], "&",value['4e']["yield"],"\\\\"
    for i, (key, value) in enumerate(TotDic.items()):
        print "\\hline", key ," & ", value['2e2m']["yield"] ," & ",value['4m']["yield"], "&",value['4e']["yield"],"\\\\"
    for i, (key, value) in enumerate(DataDic.items()):
        print "\\hline", key ," & ", value['2e2m']["yield"] ," & ",value['4m']["yield"], "&",value['4e']["yield"]," \\\\ \\hline"
    print "\\end{tabular} \n"



def YieldLatex_2(Text,SigDic,RedDic,IrrDic,TotDic,DataDic):                                                                                                                                                                                              
    Line = "\\hline"
    Line = ""
    print "\\begin{tabular}{llll}"

    print "\\hline "+Text[0]+" & "+ Text[1] + " & " + Text[2] + " & " + Text[3] + " \\\\ \\hline"   
    for i, (key, value) in enumerate(SigDic.items()):
        print " & ", value['2e2m']["yield"] ," & ",value['4m']["yield"], "&",value['4e']["yield"]
    for i, (key, value) in enumerate(RedDic.items()):
        print  " & ", value['2e2m']["yield"] ," & ",value['4m']["yield"], "&",value['4e']["yield"]
    for i, (key, value) in enumerate(IrrDic.items()):
        print  " & ", value['2e2m']["yield"] ," & ",value['4m']["yield"], "&",value['4e']["yield"]
    for i, (key, value) in enumerate(TotDic.items()):
        print " & ", value['2e2m']["yield"] ," & ",value['4m']["yield"], "&",value['4e']["yield"]
    for i, (key, value) in enumerate(DataDic.items()):
        print " & ", value['2e2m']["yield"] ," & ",value['4m']["yield"], "&",value['4e']["yield"]
    print "\\end{tabular} \n"




def SystLatex(Text,SystDicUp,SystDicDown):                                                                                                                                                                                            
    SystDic = SystDicUp

    for (keymed, valuemed),(keyup, valueup),(keydown, valuedown) in  zip(SystDic.items(),SystDicUp.items(),SystDicDown.items()):
        for (k1,valmed),(k2,valup),(k3,valdown) in zip(valuemed.items(),valueup.items(),valuedown.items()):
            if   valup["yield"] == "-" and valdown["yield"] == "-":   continue
            elif valup["yield"] == "-":  valmed["yield"] = valdown["yield"]
            elif valdown["yield"] == "-":  valmed["yield"] = valup["yield"]
            elif valup["yield"] == valdown["yield"]:  valmed["yield"] =  valup["yield"]
            elif valup["yield"] > valdown["yield"]:   valmed["yield"] =   str(valdown["yield"])+"-"+str(valup["yield"])
            else:                                     valmed["yield"] =   str(valup["yield"])+"-"+str(valdown["yield"])

    Line = "\\hline"
    Line = ""
    print "\\begin{tabular}{llll}"
    print "\\hline "+Text[0]+" & "+ Text[1] + " & " + Text[2] + " & " + Text[3] + " \\\\ \\hline"   
    for i, (key, value) in enumerate(SystDic.items()):
        print Line, key ," & ", value['2e2m']["yield"] ," \\%  &",value['4m']["yield"], " \\% &",value['4e']["yield"]," \\% \\\\ \\hline"
    print "\\end{tabular} \n"


def CrossLatex(Text,CrossDic):
    print "\\begin{tabular}{ccccc}"
    print "\\hline "+Text[0]+" & "+ Text[1]  + " & "+ Text[2] + " & "+ Text[3]+ " & "+Text[4]+ " \\\\ \\hline"   
    for i, (key, value) in enumerate(CrossDic.items()):
        if i==3: i="$\ge3$"
        print i," & ",value["cross"]," & $\pm$",value["stat"]," & $^{+",value["up"],"}_{-",value["down"],"}$ & $\pm$",value["lumi"],"  \\\\"
    print " \\hline \n \\end{tabular} \n"


    for i, (key, value) in enumerate(CrossDic.items()):

        if i==3: i="$\ge3$"
 #       print "$$\sigma_{",i,"~\mathrm{jet}}=",value["cross"],"\pm",value["stat"],"~\mathrm{(stat.)}^{+",value["up"],"}_{-",value["down"],"}~\mathrm{(syst.)} \pm",value["lumi"],"\mathrm{(lumi.)} ~\mathrm{fb} $$"
        print i,"&",value["cross"],"$\pm",value["stat"],"~\mathrm{(stat.)}^{+",value["up"],"}_{-",value["down"],"}~\mathrm{(syst.)} \pm",value["lumi"],"\mathrm{(lumi.)} $ \\\\"
