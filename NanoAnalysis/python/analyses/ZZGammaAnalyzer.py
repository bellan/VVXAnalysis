from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection
from ZZAnalysis.NanoAnalysis.tools import getLeptons, get_genEventSumw

from VVXAnalysis.NanoAnalysis.EventAnalyzer import EventAnalyzer
from VVXAnalysis.NanoAnalysis.Histogrammer import *
from VVXAnalysis.NanoAnalysis.Regions import Flags as Regions

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
            
            #to understand how ZZ Collection works 
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

            

