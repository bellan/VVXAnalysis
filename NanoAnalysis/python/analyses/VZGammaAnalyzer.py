from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection
from ZZAnalysis.NanoAnalysis.tools import getLeptons, get_genEventSumw

from VVXAnalysis.NanoAnalysis.EventAnalyzer import EventAnalyzer
from VVXAnalysis.NanoAnalysis.Histogrammer import *
from VVXAnalysis.NanoAnalysis.Regions import Flags as Regions

class VZGammaAnalyzer(EventAnalyzer, analysis_name="VZGammaAnalyzer"):

    def __init__(self, regions):
        super().__init__(regions)

    def analyze(self):
        if(self.event.HLT_passZZ4l): 
            weight = 1.

            if self.analyzeMC: self.weight = (self.event.overallEventWeight/self.genEventSumw)
            

            m4l = 125.
            self.hEvent.fill1D("ZZMass_10GeV", "ZZMass_10GeV", 93, 70., 1000., m4l, self.weight)

