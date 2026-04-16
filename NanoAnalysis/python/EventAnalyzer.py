from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection
from ZZAnalysis.NanoAnalysis.tools import getLeptons, get_genEventSumw

from VVXAnalysis.NanoAnalysis.Histogrammer import Histogrammer


class EventAnalyzer:
    def __init__(self, event, sampleName, isMC=True, genEventSumw=1.):

        self.event = event
        self.sampleName = sampleName
        self.isMC = isMC
        
        self.event.SetBranchStatus("*", 0)
        self.event.SetBranchStatus("run", 1)
        self.event.SetBranchStatus("luminosityBlock", 1)
        self.event.SetBranchStatus("*Muon*", 1)
        self.event.SetBranchStatus("*Electron*", 1)
        self.event.SetBranchStatus("*ZZCand*", 1)
        self.event.SetBranchStatus("bestCandIdx", 1)
        self.event.SetBranchStatus("HLT_passZZ4l", 1)

        if isMC:
            self.event.SetBranchStatus("overallEventWeight",1)
            self.genEventSumw = genEventSumw
        else:
            self.genEventSumw = 1.
        
        self.weight = 1.
            
        self.histogrammer = Histogrammer()

    ## Init per event quantities
    def init(self):
        pass
        
    def end(self, outFile):
         self.histogrammer.write(outFile)

    def begin(self):
        pass
