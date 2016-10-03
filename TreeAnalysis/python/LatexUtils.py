#! /usr/bin/env python

import collections



def TabLatex_n(Title):
    
#print '\\begin\n'
#print "\\begin{sidewaystable}\\footnotesize\n \\begin{center}\\scriptsize \n\\topcaption{"+Title+"} \n \\label{tab:listofsamples} \n \\begin{tabular}{lll} \n"
    
    a=3
    b=4
    c=5
    print "\\begin{tabular}{lll} \n"
    print "\\hline  \\multicolumn{3}{l}{"+Title+"}\\\\"
    print "\\hline  Process & Cross Sections [pb] & Sample\\\\"
    print "\\hline", a ," & ", b ," & ",c,"\\\\"
    
    print "\\end{tabular} \n"
    


def TabLatex(SigDic,RedDic,IrrDic,TotDic,DataDic):                                                                                                                                                                                              
    Line = "\\hline"
    Line = ""
    print "\\begin{tabular}{llll}"
    print "\\hline  Sample & 2e2mu & 4mu & 4e \\\\ \\hline"
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





