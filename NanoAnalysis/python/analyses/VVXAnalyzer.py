from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection
from ZZAnalysis.NanoAnalysis.tools import getLeptons, get_genEventSumw

from VVXAnalysis.NanoAnalysis.EventAnalyzer import EventAnalyzer
from VVXAnalysis.NanoAnalysis.Histogrammer import *
from VVXAnalysis.NanoAnalysis.Regions import Flags as Regions

class VVXAnalyzer(EventAnalyzer, analysis_name="VVXAnalyzer"):

    def __init__(self, regions):
        super().__init__(regions)

    def analyze(self):
        #print(Regions.check(self.regionWord,Regions.L4P))

    
        # Check that the event contains a selected candidate, and that
        # passes the required triggers (which is necessary for samples
        # processed with TRIGPASSTHROUGH=True)

            #weight = 1.
            #ZZs = Collection(self.event, 'ZZCand') ## move it in EventAnalyzer::init(event) ??
#            theZZ = 
            #if self.analyzeMC: self.weight = (self.event.overallEventWeight*theZZ.dataMCWeight/self.genEventSumw)
        if Regions.check(self.regionWord, Regions.ZTL2P_ZTL2P):
            m4l = self.ZZ.mass
            
            self.hEvent.fill1D("ZZMass_10GeV", "ZZMass_10GeV", 93, 70., 1000., m4l, self.weight)

        #print(Regions.check(self.regionWord, Regions.L2P_P1cutT))
