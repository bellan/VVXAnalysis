#! /usr/bin/env python

import commands, ROOT, sys
from Colours import *
from ROOT import gROOT, TTree, gSystem

gROOT.ProcessLine(".L /afs/cern.ch/work/b/bellan/lxcms131/2/CMSSW_5_3_18/src/VVXAnalysis/DataFormats/src/loader.C++")

from ROOT import phys

def convert(input_file_name):
    
    input_file = file(input_file_name, 'r')

    LorentzVector  = ROOT.Math.LorentzVector('ROOT::Math::PxPyPzE4D<double>')
    ParticleVector = ROOT.std.vector('phys::Particle')
    
    genParticles   = ParticleVector()
    genParticlesIn = ParticleVector()
    
    output_tree = TTree("Cypress"      , "Cypress"     )
    output_tree.Branch("genParticles"  , genParticles  )
    output_tree.Branch("genParticlesIn", genParticlesIn)

    in_ev   = 0 # To know when to look for particles we must know when we are inside an event
    in_ev_1 = 0 # The first line after <event> is information so we must skip that as well
    header  = True

    for line in input_file:

        if line.startswith("#"): continue

        if in_ev_1 == 1:
            in_ev_1 = 0
            in_ev = 1
            continue
    
        if line.startswith("<event>"):
            header = False
            in_ev_1 = 1
            continue
  
        if header: continue
    
        if in_ev == 1 and line.startswith("</event>"):
            output_tree.Fill()
            genParticles.clear()
            genParticlesIn.clear()
            in_ev = 0
            continue
    
        if in_ev == 1:
            gp= phys.Particle(LorentzVector(float(line.split()[6]), 
                                            float(line.split()[7]), 
                                            float(line.split()[8]), 
                                            float(line.split()[9])), 
                              phys.Particle.computeCharge(int(line.split()[0])),
                              int(line.split()[0]) )
        
            # Check the status of this particle
            if int(line.split()[1]) is 1:
                # We have a final state particle on this line
                genParticles.push_back(gp)
            elif int(line.split()[1]) is 2:
                genParticles.push_back(gp)            
            elif int(line.split()[1]) is -1:
                genParticlesIn.push_back(gp)
                
                
    output_file_name = input_file_name.replace(".lhe",".root")   
    
    f = ROOT.TFile(output_file_name,"RECREATE")
    output_tree.Write()
    f.Close()
    print "Output file is: ", Green(output_file_name)

 

input_file_name = sys.argv[1]
convert(input_file_name)
