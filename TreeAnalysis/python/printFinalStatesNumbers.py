#! /usr/bin/env python

import ROOT, copy, sys

from ROOT import TH1F

analysis = "VVXAnalyzer"
cregion = "CR2P2F"

runs = ['A','B','C','D']
datasets = ['DoubleMu','DoubleEle','MuEG']
finalstates = ['4m', '4e', '2e2m']
plots = ['ZZTo'+fs+'_mZZTo'+fs for fs in finalstates]

tableline = []

for dataset in datasets:
    for run in runs:
        fname = 'results/{0:s}_{1:s}/{2:s}{3:s}.root'.format(analysis,cregion,dataset,run)
        f = ROOT.TFile(fname)
        for finalstate in finalstates:
            plot = 'ZZTo'+finalstate+'_mZZTo'+finalstate
            h = f.Get(plot)
            if h: 
                tableline.append((dataset,run,finalstate,h.Integral())) 
                #print dataset+run, finalstate, h.Integral()

for finalstate in finalstates:        
    print "\nFinal state:", finalstate
    counter = 0
    for line in tableline:
        if line[2] == finalstate:
            print "{0:s}{1:s}\t{2:.5f}".format(line[0],line[1],line[3])
            counter += line[3]
    print "TOTAL \t\t{0:.5f}\n".format(counter)

print "STOP here, uncomment to run over MC too"
sys.exit(0)

print "MC"

datasets = [#'ZZTo4eJJ_SMHContinInterf_H125.6','ZZTo2e2muJJ_SMHContinInterf_H125.6','ZZTo4muJJ_SMHContinInterf_H125.6',
            #'ggTo4e_SMHContinInterf-MCFM67_H125.6','ggTo4mu_SMHContinInterf-MCFM67_H125.6','ggTo2e2mu_SMHContinInterf-MCFM67_H125.6',
            #'ZZTo4mu','ZZTo4e','ZZTo2mu2tau','ZZTo2e2tau','ZZTo2e2mu','ZZTo4tau',
            'ZZJetsTo4L',
            'VBFH126','ttH126','WH126','ZH126','powheg15H126',
            'DYJetsToLLTuneZ2M50','DYJetsToLLTuneZ2M10',
            #'WZ', 'WWJets', 'WGToLNuG', 'WZZJets','WWWJets','WWZJets','ZZZJets','WWGJets',
            'TTTo2L2Nu2B'
            #'TTWJets', 'TTZJets', 'TTWWJets','TTGJets'
            ]

tableline = []

for dataset in datasets:
    fname = 'results/{0:s}_{1:s}/{2:s}.root'.format(analysis,cregion,dataset)
    f = ROOT.TFile(fname)
    for finalstate in finalstates:
        plot = 'ZZTo'+finalstate+'_mZZTo'+finalstate
        h = f.Get(plot)
        if h: 
            tableline.append((dataset,finalstate,h.GetEntries(),h.Integral()))


for finalstate in finalstates:        
    print "\nFinal state:", finalstate
    counter = 0
    counterw = 0.
    for line in tableline:
        if line[1] == finalstate:
            print "{0:s}\t{1:.0f}\t{2:.3f}".format(line[0],line[2],line[3])
            counter  += line[2]
            counterw += line[3]
    print "TOTAL \t\t{0:.0f}\t{1:.3f}\n".format(counter,counterw)

