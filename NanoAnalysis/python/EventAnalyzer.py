from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection
from ZZAnalysis.NanoAnalysis.tools import getLeptons, get_genEventSumw

from VVXAnalysis.NanoAnalysis.Histogrammer import *

from VVXAnalysis.NanoAnalysis.Regions import Flags as Regions

import ROOT
import math

import subprocess

class EventAnalyzer:

    registry = {}

    def __init_subclass__(cls, analysis_name=None, **kwargs):
        super().__init_subclass__(**kwargs)

        if analysis_name is not None:

            if analysis_name in cls.registry:
                raise ValueError(
                    f"Duplicate analysis name: {analysis_name}"
                )
            cls.analysis_name = analysis_name
            cls.registry[analysis_name] = cls


    
    def __init__(self, regions):#, base_configuration):

        #self.regions = regions
        self.histogrammers = Histogrammers(regions)

        self.genEventSumw = 1.
        self.weight = 1.

        
    ## Init per event quantities
    def init(self,event,genEventSumw,luminosity = -1., isMC=False):

        self.genEventSumw = genEventSumw
        self.analyzeMC    = isMC
        self.luminosity = luminosity
        self.weight = 1.   ## weight will contain ALL corrections, including dataMC weights and FR weights
        self.dataMCWeight = 1.
        
        self.event = event
        #Turn off all branches
        self.event.SetBranchStatus("*", 0)
        
        #Read only the one that are usefull
        self.event.SetBranchStatus("run", 1)
        self.event.SetBranchStatus("luminosityBlock", 1)
        self.event.SetBranchStatus("*Lepton*", 1)
        self.event.SetBranchStatus("*Muon*", 1)
        self.event.SetBranchStatus("*Electron*", 1)
        self.event.SetBranchStatus("*Photon*", 1)
        self.event.SetBranchStatus("*Jet*", 1)
        self.event.SetBranchStatus("*ZZCand*", 1)
        self.event.SetBranchStatus("bestCandIdx", 1)
        self.event.SetBranchStatus("HLT_passZZ4l", 1)
        self.event.SetBranchStatus("regionWord", 1)
        self.event.SetBranchStatus("*GenPart*", 1)
        self.event.SetBranchStatus("*GenZZ*", 1)
        
        if self.analyzeMC:
            self.event.SetBranchStatus("overallEventWeight",1)


    def eventSetup(self):

            
        self.ZZs        = Collection(self.event, 'ZZCand')
        self.electrons  = Collection(self.event, 'Electron')
        self.muons      = Collection(self.event, 'Muon')
        self.leptons    = list(self.electrons) + list(self.muons)
        
        self.regionWord = self.event.regionWord
        self.hEvent = self.histogrammers.checkRegions(self.regionWord)


        if self.analyzeMC:
            # lumi is expressed in 1/fb, while xsections are in pb
            self.weight = 10e3*self.luminosity*self.event.overallEventWeight/self.genEventSumw
            # compute the dataMC correction. Store it in a separate weight, for checks
            self.dataMCWeight = 1.
            self.dataMCWeight *= math.prod(lep.dataMC for lep in self.leptons if lep.ZZFullSel)
            self.weight *= self.dataMCWeight
            
        ## Analyze the type of event and prepare the composite particles consquently
        if Regions.check(self.regionWord, Regions.L4P):
            if(self.event.bestCandIdx != -1 and self.event.HLT_passZZ4l):
                self.ZZ = self.ZZs[self.event.bestCandIdx]
                
                #self.regionWord |= Regions.ZTL2P_ZTL2P

        
        
        
        
    def end(self, sample):
        
        # for region in regions:
        #     odir = outputdir_format %region
        #     outputdirs[region] = odir
        #     subprocess.check_call(['mkdir', '-p', odir])  # Use os.makedirs(odir, exist_ok=True) after switching to py3
        
 #       odir = f"results/{sample.year}" 
 #       subprocess.check_call(['mkdir', '-p', odir])  # Use os.makedirs(odir, exist_ok=True) after switching to py3
        
 #       outFile = ROOT.TFile.Open(
#            f"{odir}/{self.analysis_name}.root",
#            "recreate")
        
        self.histogrammers.write(
            base_odir = f"results/{sample.year}",
            analysis_name = self.analysis_name,
            sample_name   = sample.name
        )

    def begin(self):
        pass
