from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection
from ZZAnalysis.NanoAnalysis.tools import getLeptons, get_genEventSumw

from VVXAnalysis.NanoAnalysis.EventAnalyzer import EventAnalyzer
from VVXAnalysis.NanoAnalysis.Histogrammer import *
from VVXAnalysis.NanoAnalysis.Regions import Flags as Regions

import math

class ZZGammaAnalyzer(EventAnalyzer, analysis_name="ZZGammaAnalyzer"):

    def __init__(self, regions):
        super().__init__(regions)

    def analyze(self):
        bestCandIdx = self.event.bestCandIdx
    
        # Check that the event contains a selected candidate, and that
        # passes the required triggers (which is necessary for samples
        # processed with TRIGPASSTHROUGH=True)
        if(bestCandIdx != -1 and self.event.HLT_passZZ4l): 
            weight = 1.
            
            # Collections

            #GenZZs = Collection(self.event, 'GenZZ')
            GenParts = Collection(self.event, 'GenPart')
            GenLeptons = [p for p in GenParts if abs(p.pdgId) == 11 or abs(p.pdgId == 13)]
            GenPhotons = [p for p in GenParts if abs(p.pdgId) == 22]

            ZZs = Collection(self.event, 'ZZCand') ## move it in EventAnalyzer::init(event) ??
            theZZ = ZZs[bestCandIdx]
            Leptons = Collection(self.event, 'Lepton')
            Photons = Collection(self.event, 'Photon')
            FsrPhotons = Collection(self.event, 'FsrPhoton')

            if self.analyzeMC: self.weight = (self.event.overallEventWeight*theZZ.dataMCWeight/self.genEventSumw)

            
            
           