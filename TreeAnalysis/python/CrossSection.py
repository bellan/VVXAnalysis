#! /usr/bin/env python

import ROOT,copy
from ROOT import gSystem, TCanvas, TH1,  TPad, gStyle, TLegend
#from collections import OrderedDict
from optparse import OptionParser
import sys
import math
import operator


def GetNevAndErr(samples,inputdir,FinState):
    for sample in samples:
        print "opening",inputdir+sample+".root"
        files[sample] = ROOT.TFile(inputdir+sample+".root")  

    hdata=ROOT.TH1F()
    isFirst=1
    print 'ZZTo'+FinState+'_Mass'
    for dat in samples:
        h = files[dat].Get('ZZTo'+FinState+'_Mass')
        #h = files[dat].Get('ZZMass')
        if h==None:
            print dat," has no entries or is a zombie" 
            continue
        if isFirst:
        #    print dat, ", total number of events: ", h.Integral()," entries ",h.GetEntries()
            hdata=copy.deepcopy(h)
            isFirst=0
            continue
       # print dat, ", total number of events: ", h.Integral()," entries ",h.GetEntries()
        hdata.Add(h)
    #hdata.Draw()        
    #c2.SaveAs('test.png')
    Err=ROOT.Double(0.)
    Integr= hdata.IntegralAndError(-1,-1,Err)
    return (Integr,Err)

parser = OptionParser(usage="usage: %prog <analysis> <sample> [options]")
parser.add_option("-f", "--finstate", dest="FinState",
                  default='4l',
                  help="Chose finalstate 4l, 2e2mu, 4e, 4m")

(options, args) = parser.parse_args()

FinState=options.FinState

data = ['DoubleMuA','DoubleMuB','DoubleMuC','DoubleMuD','DoubleEleA','DoubleEleB','DoubleEleC','DoubleEleD','MuEGA','MuEGB','MuEGC','MuEGD']
bkgIrr = ['TTZJets','WWZJets','TTWWJets']
bkgRed= ['reducible_background_from_DoubleEle','reducible_background_from_DoubleMu','reducible_background_from_MuEG']

inputdir="results/ZZjAnalyzer_SR/"
inputdirRed="results/ZZjAnalyzer_CR/"

files={}

(IntegrData,ErrorData) = GetNevAndErr(data,inputdir,FinState)
(IntegrIrr,ErrorIrr)   = GetNevAndErr(bkgIrr,inputdir,FinState)
(IntegrRed,ErrorRed)   = GetNevAndErr(bkgRed,inputdirRed,FinState)

print 'N events data',IntegrData,' +- ',ErrorData,'\nN events bkg irr ',IntegrIrr,' +- ',ErrorIrr,'\nN events bkg Red',IntegrRed,' +- ',ErrorRed

ErrTriggEffRel = 0.015 #%
ErrIsoIdeRel =0.015

Lumi=19712 #pb-1
ErrLumiRel = 0.044 #%

ErrAccRel = 0.05  
ErrEffRel = 0.015 # 

BRele = 0.03363
BRmu = 0.03366

if FinState=='4e':
    Acc = 0.35323 #%
    BR=BRele*BRele
    ErrRed=0.71524/IntegrRed 
elif FinState=='4m':
    Acc = 0.4835 #%
    BR=BRmu*BRmu
    ErrRed=0.60674/IntegrRed 
elif FinState=='2e2m':
    Acc=0.4133
    BR=2*BRmu*BRele
    ErrRed = 1.094/IntegrRed 
elif FinState=='4l':
    Acc=0.4255
    BR = (BRele+BRmu)*(BRele+BRmu)
    ErrRed = 1.4442495974/IntegrRed 


#ErrorData/=IntegrData
#ErrorIrr/=IntegrIrr

#ErrNum = math.sqrt(ErrorData*ErrorData)

Cross = (IntegrData-IntegrIrr-IntegrRed)/(Lumi*Acc*BR)

#ErrStat = (ErrorData*ErrorData+ErrorIrr*ErrorIrr)
#ErrStat /= (Lumi*Acc*Lumi*Acc)*BR
#ErrStat = math.sqrt(ErrStat)
print 'ErrRid ',ErrRed,'\n'

ErrStat= (ErrorData)/(Lumi*Acc*BR)
ErrSist = math.sqrt((ErrLumiRel*ErrLumiRel+ErrAccRel*ErrAccRel+ErrTriggEffRel*ErrTriggEffRel+ErrIsoIdeRel*ErrIsoIdeRel)*Cross*Cross+(ErrorIrr*ErrorIrr+ErrRed*ErrRed)/(Lumi*Lumi*Acc*Acc*BR*BR))
#CloseVar=input("digit anything you want to end the script \n")        

print 'CrossSection = ',Cross,'+-',ErrStat,'(stat) +- ',ErrSist,'(sist)'
