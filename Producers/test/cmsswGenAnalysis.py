#! /usr/bin/env python

import ROOT

from ROOT import TH1F, gSystem, gPad

gSystem.Load("libFWCoreFWLite")
from ROOT import AutoLibraryLoader
AutoLibraryLoader.enable()
gSystem.Load("libDataFormatsFWLite")

import math
import sys, os, copy
from DataFormats.FWLite import Events, Handle
from array import array
import FWCore.ParameterSet.Config as cms

genParticlesH = Handle("std::vector<reco::GenParticle>")
genLabel       = "genParticles"
files = []

files.append("ZZ_MGNLO_py_LHE_GEN_1.root")


events = Events(files)
maxevents = 100
counter = 0
for event in events:
    
    if counter == maxevents: break
    counter += 1

    event.getByLabel(genLabel,genParticlesH)
    genParticles = genParticlesH.product()
    print "-----------------------------------------------------------"
    for gen in genParticles:
        if gen.status() == 1 and (abs(gen.pdgId()) == 11 or abs(gen.pdgId()) == 13 or abs(gen.pdgId() == 15)) and gen.pt() > 10: print gen.pdgId(), gen.pt()
        
