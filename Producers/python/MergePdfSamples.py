#! /usr/bin/env python
import ROOT

from ROOT import gSystem
gSystem.Load("libFWCoreFWLite")
from ROOT import AutoLibraryLoader
AutoLibraryLoader.enable()

import sys
import math
import operator
from PDFSysEvaluationTools import *
from ZZAnalysis.AnalysisStep.readSampleInfo import *


folder = sys.argv[1]
samples=["ZZTo4mu","ZZTo4e","ZZTo2e2mu","ZZTo4muJJ_SMHContinInterf","ZZTo4eJJ_SMHContinInterf","ZZTo2e2muJJ_SMHContinInterf","ggTo2e2mu_SMHContinInterf-MCFM67","ggTo4mu_SMHContinInterf-MCFM67","ggTo4e_SMHContinInterf-MCFM67","ZZZJets","ttH126","WH126","ZH126","VBFH126","WZZJets","powheg15H126"]

FinStates=["2e2mu","4mu","4e","Tot"]
PDFs=["1_CT10","2_MSTW2008nlo68cl","3_NNPDF20"]

FileOut=ROOT.TFile("Tot.root","recreate")
FileOut.mkdir("pdfSystematics")
FileOut.Close()
for pdf in PDFs:
    print 'PDF',pdf,'\n'
    for fin in FinStates:
        hTot = ROOT.TH2F("htot","TOT",55, 0., 55.,8,1.,9.)
        print "final state",fin,"\n"
        for s in samples:
            print 'opening',sys.argv[1]+s+".root"
            file =  ROOT.TFile(sys.argv[1]+s+".root")

            name="pdfSystematics/hPdf"+fin+"Set"+pdf
            h1 = file.Get(name)
            hTot.SetName(h1.GetName())
            Nev= h1.GetBinContent(1,8)
            if "_" in s:
                s+="_H125.6"
            Cross=crossSection(s,'../python/samples_8TeV.csv') 
            wh=Cross/Nev
            print "cross",Cross,"N ev",Nev,"Weight",wh,'\n'
            hTot.Add(h1,wh)
            file.Close()
        FileOut=ROOT.TFile("Tot.root","update")
        FileOut.cd("pdfSystematics")
        print hTot.GetEntries()
        hTot.Write()
        hTot.Delete()
        FileOut.Close()    
