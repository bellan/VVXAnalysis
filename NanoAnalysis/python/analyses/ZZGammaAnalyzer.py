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
            
            # ZZmass pre and post FSR 
            ZZs = Collection(self.event, 'ZZCand') ## move it in EventAnalyzer::init(event) ??
            theZZ = ZZs[bestCandIdx]
            if self.analyzeMC: self.weight = (self.event.overallEventWeight*theZZ.dataMCWeight/self.genEventSumw)
            
            mZZ = theZZ.mass
            self.hEvent.fill1D("ZZMass_10GeV", "ZZMass_10GeV", 93, 70., 1000., mZZ, self.weight)

            mZZPreFSR = theZZ.massPreFSR
            self.hEvent.fill1D("ZZMassPreFSR_10GeV", "ZZMassPreFSR_10GeV", 93, 70., 1000., mZZPreFSR, self.weight)
            
            #GenZZ = Collection(self.event, 'GenZZ')
            #mGenZZ = GenZZ.mass
            #self.hevent.fill1D("GenZZMass_10GeV", "GenZZMass_10GeV", 93, 70., 1000., mGenZZ, self.weight)

            #definizione di segnale
            GenParts = Collection(self.event, 'GenPart')
            GenLeptons = [p for p in GenParts if abs(p.pdgId) == 11 or abs(p.pdgId == 13)]
            GenPhotons = [p for p in GenParts if abs(p.pdgId) == 22]

            def ZZGammaSignalDefinition():
                if not self.analyzeMC: return False
                
                # ZZ requirements
                if len(theZZ) != 1 or len(GenLeptons) != 4: return False
                if sum(l.pt > 5 for l in GenLeptons) < 4: return False
                if sum(l.pt > 10 for l in GenLeptons) < 2: return False
                if sum(l.pt > 20 for l in GenLeptons) < 1: return False
                if any(abs(l.eta) > 2.5 for l in GenLeptons): return False
                
                # photon requirements
                if len(GenPhotons) < 1: return False
                if sum(p.pt > 20 for p in GenPhotons) < 1: return False
                if sum(abs(p.eta) <=2.4 and not 1.444 < abs(p.eta) < 1.566 for p in GenPhotons) < 1: return False
                

         