#! /usr/bin/env python
from __future__ import print_function
import CrossInfo
from   CrossInfo import*
import collections

class Latex:
    def __init__(self, Type, other):
        self.Type = Type
        self.other = other
        self.Dic = collections.OrderedDict()
        self.sampleTypeStatus = ""
        self.sampleStatus = ""

    def printStatus(self):
        print(self.Type, self.other)
        print(self.Dic)
        print("sampleTypeStatus", self.sampleTypeStatus)
        print("sampleStatus"    , self.sampleStatus    )

    def addElement(self,sampleType,sample,channel,obj,value):
        if sampleType not in self.Dic:
            self.Dic[sampleType] = {}
        if sample not in self.Dic[sampleType]:
            self.Dic[sampleType][sample] = {}
            self.Dic[sampleType][sample] = {'2e2m':{},'4m':{},'4e':{}}

        self.Dic[sampleType][sample][channel][obj]=value

    def setSampleTypeStatus(self,status):
        self.sampleTypeStatus = status

    def setSampleStatus(self,status):
        self.sampleStatus = status


    def add(self,channel,obj,value):
        if self.sampleTypeStatus not in self.Dic:
            self.Dic[self.sampleTypeStatus] = {}
        if self.sampleStatus not in self.Dic[self.sampleTypeStatus]:
            self.Dic[self.sampleTypeStatus][self.sampleStatus] = {}
            self.Dic[self.sampleTypeStatus][self.sampleStatus] = {'2e2m':{},'4m':{},'4e':{}}

        self.Dic[self.sampleTypeStatus][self.sampleStatus][channel][obj]=value

    def addBin(self,bin,channel,obj,value):
        if self.sampleTypeStatus not in self.Dic:
            self.Dic[self.sampleTypeStatus] = {}
        if self.sampleStatus not in self.Dic[self.sampleTypeStatus]:
            self.Dic[self.sampleTypeStatus][self.sampleStatus] = {}
        if bin not in self.Dic[self.sampleTypeStatus][self.sampleStatus]:
            self.Dic[self.sampleTypeStatus][self.sampleStatus][bin] = {'2e2m':{},'4m':{},'4e':{}}
        if obj not in  self.Dic[self.sampleTypeStatus][self.sampleStatus][bin][channel]:
            self.Dic[self.sampleTypeStatus][self.sampleStatus][bin][channel][obj]=value
        else:
            self.Dic[self.sampleTypeStatus][self.sampleStatus][bin][channel][obj]+=value

    def setBin(self,bin,channel,obj,value):
        if self.sampleTypeStatus not in self.Dic:
            self.Dic[self.sampleTypeStatus] = {}
        if self.sampleStatus not in self.Dic[self.sampleTypeStatus]:
            self.Dic[self.sampleTypeStatus][self.sampleStatus] = {}
        if bin not in self.Dic[self.sampleTypeStatus][self.sampleStatus]:
            self.Dic[self.sampleTypeStatus][self.sampleStatus][bin] = {'2e2m':{},'4m':{},'4e':{}}
        self.Dic[self.sampleTypeStatus][self.sampleStatus][bin][channel][obj]=value


    def setSystBin(self,bin,channel,obj,value):
        if self.sampleTypeStatus not in self.Dic:
            self.Dic[self.sampleTypeStatus] = {}
        if self.sampleStatus not in self.Dic[self.sampleTypeStatus]:
            self.Dic[self.sampleTypeStatus][self.sampleStatus] = {}
        if bin not in self.Dic[self.sampleTypeStatus][self.sampleStatus]:
            self.Dic[self.sampleTypeStatus][self.sampleStatus][bin] = {'2e2m':{},'4m':{},'4e':{},'4l':{}}
        self.Dic[self.sampleTypeStatus][self.sampleStatus][bin][channel][obj]=value



    def printCode(self):
#def YieldLatex(Text,SigDic,RedDic,IrrDic,TotDic,DataDic):                                                                                                                                                                                    
        Text =("Sample","2e2mu" ,"4mu" ,"4e")
    
        Line = "\\hline"
        Line = ""
        print("\\begin{tabular}{llll}")
        
        print("\\hline "+Text[0]+" & "+ Text[1] + " & " + Text[2] + " & " + Text[3] + " \\\\ \\hline")
        for i, (key, value) in enumerate(self.Dic["Signal"].items()):
            print(Line, key ," & ", value['2e2m']["yield"] ," & ",value['4m']["yield"], "&",value['4e']["yield"],"\\\\")
        for i, (key, value) in enumerate(self.Dic["Irr"].items()):
            print(Line, key ," & ", value['2e2m']["yield"] ," & ",value['4m']["yield"], "&",value['4e']["yield"],"\\\\")
        for i, (key, value) in enumerate(self.Dic["red"].items()):
            print(Line, key ," & ", value['2e2m']["yield"] ," & ",value['4m']["yield"], "&",value['4e']["yield"],"\\\\")
        for i, (key, value) in enumerate(self.Dic["red"].items()):
            print("\\hline", key ," & ", value['2e2m']["yield"] ," & ",value['4m']["yield"], "&",value['4e']["yield"],"\\\\")
        for i, (key, value) in enumerate(self.Dic["data"].items()):
            print("\\hline", key ," & ", value['2e2m']["yield"] ," & ",value['4m']["yield"], "&",value['4e']["yield"]," \\\\ \\hline")

        print("\\end{tabular}\n")


    def printPerBinCode(self):

        self.Dic["Total"] = {}
    
        self.setSampleTypeStatus("Total")
        self.setSampleStatus("Total")
        #    self.addBin(1,"4e","yield",2.2)                                                                                                                                                                                              
        print(self.Dic)
        
        Line = ""
        print(len( self.Dic["Data"]["Data"]))
        print("\\begin{tabular}{"+(len( self.Dic["Data"]["Data"] ))*"lll"+"ll}")
        print("\\hline   ", end=' ')
        
        for bin in range(0,len( self.Dic["Signal"] )+1): print(" & \\multicolumn{3}{c}{$\geq $", bin , " jet}", end="")
        print("\\\\ \\hline")
        print(" sample ", end=' ')
        for bin in range(0,len( self.Dic["Signal"] )+1): print(" & $2e2\mu$ & $4\mu$ & $4e$ ", end="")
        print("\\\\ \\hline")
        for Type in ("Signal","Background","Total","Data"):
            if Type == "Data" or Type =="Total":  print(" \\hline")
            for i, (key, value) in enumerate(self.Dic[Type].items()):
                print(Line, key ," & ", end="")
                for bin in range(1,len(value)+1):
                    if Type != "Data" and Type != "Total":
                        self.addBin(bin,"2e2m","yield",value[bin]['2e2m']["yield"])
                        self.addBin(bin,"4m","yield",value[bin]['4m']["yield"])
                        self.addBin(bin,"4e","yield",value[bin]['4e']["yield"])
            
                    print("{0:.2f} & {1:.2f} & {2:.2f} &".format( value[bin]['2e2m']["yield"],value[bin]['4m']["yield"],value[bin]['4e']["yield"]), end=" ")
                print("\\\\")
        print("\\hline")
        print("\\end{tabular}\n")


    def printSystRanges(self):

        print("\\begin{tabular}{lc}")
        print("Systematic source &   \\\\")
        print("\\hline")
        for i, (key, value) in enumerate(self.Dic["Syst"].items()):
            List1 = []
        
            for bin in range(1,len(value)+1):
            #    print value[bin]["4l"][1],"\n\n"
                List1.append(value[bin]["4l"][1])
                List1.append(value[bin]["4l"][-1])
            #print "{0} {1:.2f} - {2:.2f}".format(key,min(List1),max(List1))
        
            if min(List1) == max(List1):  print("{0} & {1:.2f} \\%  \\\\".format(key,min(List1),max(List1)))
            else:                         print("{0} &  {1:.2f} - {2:.2f} \\% \\\\".format(key,min(List1),max(List1)))
        for glob in GlobSystList:
            print("{0} & {1:.2f} \\% \\\\".format(glob["name"],glob["value"]*100))
#        for i, (key, value) in enumerate(GlobSystList):
        print("\\hline")
        print("\\end{tabular}\n")
