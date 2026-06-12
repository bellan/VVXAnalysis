from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module
from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection

from VVXAnalysis.NanoAnalysis.Selector import Selector
from VVXAnalysis.NanoAnalysis.Regions import Flags


class FSEventTaggerAndFilter(Module):

    def __init__(self, flags, regions):

        self.flags   = flags
        self.regions = regions
        
    def check(self, collection, selection):

        if selection == None:
            return True # because there is not requirement on this collection
        
        s = Selector(selection)
        return s.applySelection(collection)

    
    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        self.out = wrappedOutputTree

        self.out.branch("regionWord", "I", title="Word that contains the regions that passed the selection")
    
        
    def analyze(self, event):
        """process event, return True (go to next module) or False (fail, go to next event)"""

        # Collections
        #electrons = Collection(event, "Electron")
        #muons     = Collection(event, "Muon")
        photons   = Collection(event, "Photon")
        leptons   = Collection(event, "Lepton")
        jets      = Collection(event, "Jet") # FIXME: for the time being, AK4 only

        regionWord = 0
        for flag in self.flags:
            print(flag)
            if self.check(leptons, flag.get("leptons")) and self.check(photons,flag.get("photons")) and self.check(jets, flag.get("jets")):
                name = flag.get("name")
                if name in Flags.__members__:
                    regionWord |= Flags[name]

        
        self.out.fillBranch("regionWord", regionWord)

        return not regionWord == 0 
