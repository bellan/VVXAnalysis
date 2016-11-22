#! /usr/bin/env python
import os
import sys

""" Mv root files from AAAOK to samples folder """

from os import walk

Targetfile= sys.argv[1]

What= sys.argv[2]
#datalist = ["DoubleEG2016B","DoubleMu2016B","MuonEG2016B","SingleEle2016B","SingleMuon2016B","DoubleEG2016C","DoubleMu2016C","MuonEG2016C","SingleEle2016C","SingleMuon2016C","DoubleEG2016D","DoubleMu2016D","MuonEG2016D","SingleEle2016D","SingleMuon2016D","DoubleEG2016E","DoubleMu2016E","MuonEG2016E","SingleEle2016E","SingleMuon2016E","DoubleEG2016F","DoubleMu2016F","MuonEG2016F","SingleEle2016F","SingleMuon2016G","DoubleEG2016G","DoubleMu2016G","MuonEG2016G","SingleEle2016G","SingleMuon2016F","DoubleEG2016H","DoubleMu2016H","MuonEG2016H","SingleEle2016H","SingleMuon2016H"]

datalist = ["DoubleEG2016C","DoubleMu2016C","MuonEG2016C","SingleEle2016C","SingleMuon2016C","DoubleEG2016D","DoubleMu2016D","MuonEG2016D","SingleEle2016D","SingleMuon2016D"]

datalist_2015 = ["DoubleEG2015C","DoubleMu2015C","MuonEG2015C","SingleEle2015C","DoubleEG2015D","DoubleMu2015D","MuonEG2015D","SingleEle2015D"]

MCList = ["ZZTo4l","ZZTo4lamcatnlo","ggZZ4e","ggZZ4mu","ggZZ2e2mu","ggZZ4tau","ggZZ2e2tau","ggZZ2mu2tau","DYJetsToLL_M50","TTJets","TTTo2L2Nu","WZ","ZGTo2LG","ZH125","ttH125","VBFH125","ZZTo2e2muJJ","ZZTo4muJJ","ZZTo4eJJ","WWZ","TTGJets","TTWJets","WWW","ZZZ","WZZ"]


AddbashCommand="hadd -f "+Targetfile+"/Data.root "
bashCommand="mv "

if What=="All" or What=="Data":
    for Dir in datalist:
        AddbashCommand=AddbashCommand+" "+Targetfile+"/"+Dir+".root"
        print bashCommand+" "+Dir+"/ZZ4lAnalysis.root "+Targetfile+"/"+Dir+".root"
        os.system(bashCommand+" "+Dir+"/ZZ4lAnalysis.root "+Targetfile+"/"+Dir+".root")

        print AddbashCommand
        os.system(AddbashCommand)
 
if What=="All" or What=="MC":
    for s in MCList:
        print bashCommand+" "+s+"/ZZ4lAnalysis.root "+Targetfile+"/"+s+".root"
        os.system(bashCommand+" "+s+"/ZZ4lAnalysis.root "+Targetfile+"/"+s+".root")
        
