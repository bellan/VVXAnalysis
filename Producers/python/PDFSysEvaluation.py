#! /usr/bin/env python

import ROOT

from ROOT import gSystem
gSystem.Load("libFWCoreFWLite")
from ROOT import AutoLibraryLoader
AutoLibraryLoader.enable()

import sys
import math
import operator

def evaluateMax(a,b,c):   
    max = 0.
    if( a > b ):
        if( a > c ):
            max = a
        else :
            max = c
    else:
        if( b > c ):
            max = b
        else :
            max = c
#    if not max == 0 : print "For the MAX. There are {0:.3f}, {1:.3f}, {2:.3f}. The MAX is {3:.3f}".format(a,b,c,max)
    return max


def evaluateMin(a,b,c):
    min = 0.
    if( a < b ):
        if( a < c ):
            min = a
        else :
            min = c
    else:
        if( b < c ):
            min = b
        else :
            min = c
#    if not min == 0 : print "For the MIN. There are {0:.3f}, {1:.3f}, {2:.3f}. The MIN is {3:.3f}".format(a,b,c,min)
    return min

def Prec(prec,a):
    Corr = int((a * prec) + 0.5) /prec
    return Corr

def GetPdfResult(fileIn,i):

    nnpdfFlag = False
    if i==1: h1 = file.Get("pdfSystematics/hPdf_CT10")
    elif i==2: h1 = file.Get("pdfSystematics/hPdf_MSTW2008nlo68cl")
    elif i==3: h1 = file.Get("pdfSystematics/hPdf_NNPDF20")
    
    Wh1=[]
    WhSel1=[]
    Wh2=[]
    WhSel2=[]
    
    Orig= h1.GetBinContent(1,8)
    TotWhSel= h1.GetBinContent(1,6)
    TotWh= h1.GetBinContent(1,7)

    npars=0
    acc_central = TotWhSel/TotWh
    for j in range(1,55):
        
        if h1.GetBinContent(j,1)==0: break
        npars+=1        
        Wh1.append(h1.GetBinContent(j,1))
        WhSel1.append(h1.GetBinContent(j,2))
        Wh2.append(h1.GetBinContent(j,3))
        WhSel2.append(h1.GetBinContent(j,4))

    if i==3:
        nnpdfFlag = True

    wplus=0        
    wminus=0
    nplus=0
    nminus=0
    
    for j in range(0,npars):
        wa = 0.
	if Wh1[j]>0 : wa = WhSel1[j]/Wh1[j]/acc_central-1.
        wb = 0.
        if WhSel2[j]>0 : wb =  WhSel2[j]/Wh2[j]/acc_central-1.
        if nnpdfFlag :
            if wa>0.:
                wplus += wa*wa 
                nplus+=1
            else:
                wminus += wa*wa
                nminus+=1
                    
            if wb>0.: 
                wplus += wb*wb 
                nplus+=1
            else:
                wminus += wb*wb
                nminus+=1
        else:
            if wa>wb :
                if wa<0.: wa = 0.
                if wb>0.: wb = 0.
                wplus += wa*wa
                wminus += wb*wb
            else:
                if wb<0.: wb = 0.
                if wa>0.: wa = 0.
                wplus += wb*wb
                wminus += wa*wa
                  
    if wplus>0:
        wplus = math.sqrt( wplus )
    if wminus>0:
        wminus = math.sqrt( wminus )
    if nnpdfFlag:
        if nplus>0: wplus /= math.sqrt( nplus )
        if nminus>0: wminus /= math.sqrt( nminus )
        
    HistoName=h1.GetName()
    PdfName=HistoName[5:]

    listAcc = []
    listAcc.append(PdfName)
    listAcc.append(Orig)
    listAcc.append(TotWh)
    listAcc.append(TotWhSel)
    listAcc.append(acc_central)
    listAcc.append(wplus)
    listAcc.append(wminus)
    return listAcc

file =  ROOT.TFile(sys.argv[1])

Acc1 = GetPdfResult(file,1)
Acc2 = GetPdfResult(file,2)
Acc3 = GetPdfResult(file,3)

print " opening ",sys.argv[1]
prec = math.pow(10,-1+abs(int(math.log10(abs(Acc1[5]/1.645)*Acc1[5]))))

print "\n","######################################## \n","######### PDF Systematic Errors ######## \n","######################################## \n \n"
print "######### Pdf ",Acc1[0],"######### \n Original Events  = ",Acc1[2]," \n Selected Events  = ",Acc1[3] ,"\n Acceptance      = ",Prec(prec,Acc1[4])," + ", Prec(prec,Acc1[5]/1.645*Acc1[4])," - ", Prec(prec,Acc1[6]/1.645*Acc1[4]),"  \n"
print "######### Pdf ",Acc2[0],"######### \n Original Events  = ",Acc2[2]," \n Selected Events  = ",Acc2[3] ,"\n Acceptance      = ", Prec(prec,Acc2[4])," + ", Prec(prec,Acc2[5]*Acc2[4])," - ", Prec(prec,Acc2[6]*Acc2[4]),"  \n"
print "######### Pdf ",Acc3[0],"######### \n Original Events  = ",Acc3[2]," \n Selected Events  = ",Acc2[3] ,"\n Acceptance      = ", Prec(prec,Acc3[4])," + ", Prec(prec,Acc3[5]*Acc3[4])," - ", Prec(prec,Acc3[6]*Acc3[4]),"  \n"

CentralValue = 0.5*( evaluateMax(Acc1[4]*(1 + Acc1[5]/1.645) ,Acc2[4]*( 1 + Acc2[5]) , Acc3[4]*(1+Acc3[5]) ) + evaluateMin(Acc1[4]*(1 - Acc1[5]/1.645),Acc2[4]*(1 - Acc2[5]) , Acc3[4]*( 1- Acc3[5])))

Err = 0.5*( evaluateMax(Acc1[4]*(1 + Acc1[5]/1.645) ,Acc2[4]*( 1 + Acc2[5]) , Acc3[4]*(1+Acc3[5]) ) - evaluateMin(Acc1[4]*(1 - Acc1[5]/1.645),Acc2[4]*(1 - Acc2[5]) , Acc3[4]*( 1- Acc3[5])))

print "Envelope Value = ", Prec(prec,CentralValue)," +/- ", Prec(prec,Err),"\n "

    
