from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module
from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection

from VVXAnalysis.NanoAnalysis.Selector import Selector
from VVXAnalysis.NanoAnalysis.Regions import Regions


class VVXEventTaggerAndFilter(Module):

    def __init__(self, regions):
        self.regions = regions
    
    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        self.out = wrappedOutputTree

        self.out.branch("regionWord", "I", title="Word that contains the regions that passed the selection")
    
        
    def analyze(self, event):
        """process event, return True (go to next module) or False (fail, go to next event)"""

        # Collections
        #electrons = Collection(event, "Electron")
        #muons     = Collection(event, "Muon")
        photons   = Collection(event, "Photon")

        regionWord = 0
        for region in self.regions:
           
            if self.check(leptons, region.get("leptons")) and self.check(photons,region.get("photons")) and self.check(jets, region.get("jets")):
                name = region.get("name")
                if name in Regions.__members__:
                    regionWord |= Regions[name].value

        
        self.out.fillBranch("regionWord", regionWord)

        return not regionWord == 0 
