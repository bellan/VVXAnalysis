from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection
from ZZAnalysis.NanoAnalysis.tools import getLeptons, get_genEventSumw

from VVXAnalysis.NanoAnalysis.Histogrammer import *

from VVXAnalysis.NanoAnalysis.Regions import Flags as Regions

import ROOT

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
    def init(self,event,genEventSumw,isMC=False):

        self.genEventSumw = genEventSumw
        self.analyzeMC    = isMC
        
        self.event = event
        #Turn off all branches
        self.event.SetBranchStatus("*", 0)
        
        #Read only the one that are usefull
        self.event.SetBranchStatus("run", 1)
        self.event.SetBranchStatus("luminosityBlock", 1)
        self.event.SetBranchStatus("*Muon*", 1)
        self.event.SetBranchStatus("*Electron*", 1)
        self.event.SetBranchStatus("*Photon*", 1)
        self.event.SetBranchStatus("*Jet*", 1)
        self.event.SetBranchStatus("*ZZCand*", 1)
        self.event.SetBranchStatus("bestCandIdx", 1)
        self.event.SetBranchStatus("HLT_passZZ4l", 1)
        self.event.SetBranchStatus("regionWord", 1)
        
        if self.analyzeMC:
            self.event.SetBranchStatus("overallEventWeight",1)


    def getCollections(self):
        self.ZZs = Collection(self.event, 'ZZCand') 
        self.regionWord = self.event.regionWord
        self.hEvent = self.histogrammers.checkRegions(self.regionWord)
        
        
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
