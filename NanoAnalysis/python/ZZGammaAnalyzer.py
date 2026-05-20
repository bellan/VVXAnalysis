from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection
from ZZAnalysis.NanoAnalysis.tools import getLeptons, get_genEventSumw

from VVXAnalysis.NanoAnalysis.EventAnalyzer import EventAnalyzer
from VVXAnalysis.NanoAnalysis.Histogrammer import Histogrammer

class ZZGammaAnalyzer(EventAnalyzer):

    def __init__(self,event, sampleName, isMC=True, genEventSumw=1.):
        super().__init__(event, sampleName, isMC, genEventSumw)

    def analyze(self):
        bestCandIdx = self.event.bestCandIdx
    
        # Check that the event contains a selected candidate, and that
        # passes the required triggers (which is necessary for samples
        # processed with TRIGPASSTHROUGH=True)
        if(bestCandIdx != -1 and self.event.HLT_passZZ4l): 
            weight = 1.
            ZZs = Collection(self.event, 'ZZCand') ## move it in EventAnalyzer::init(event) ??
            theZZ = ZZs[bestCandIdx]
            if self.isMC: self.weight = (self.event.overallEventWeight*theZZ.dataMCWeight/self.genEventSumw)
            
            m4l = theZZ.mass
            
            self.histogrammer.fill1D("ZZMass_10GeV_"+self.sampleName, "ZZMass_10GeV_"+self.sampleName, 93, 70., 1000., m4l, self.weight)

