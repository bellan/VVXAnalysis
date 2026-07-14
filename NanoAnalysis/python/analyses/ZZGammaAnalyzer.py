from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection
from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Object
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

            GenParts = Collection(self.event, 'GenPart')
            GenLeptons = [p for p in GenParts if abs(p.pdgId) == 11 or abs(p.pdgId == 13)]
            GenPhotons = [p for p in GenParts if abs(p.pdgId) == 22]

            ZZs = Collection(self.event, 'ZZCand') ## move it in EventAnalyzer::init(event) ??
            theZZ = ZZs[bestCandIdx]
            Leptons = Collection(self.event, 'Lepton')
            Photons = Collection(self.event, 'Photon')
            FsrPhotons = Collection(self.event, 'FsrPhoton')

            if self.analyzeMC: self.weight = (self.event.overallEventWeight*theZZ.dataMCWeight/self.genEventSumw)

            # ZZMass pre and post FSR            
            
            ZZMass = theZZ.mass
            self.hEvent.fill1D("ZZMass_10GeV", "ZZMass_10GeV", 93, 20., 1000., ZZMass, self.weight)

            ZZMassPreFsr = theZZ.massPreFSR
            self.hEvent.fill1D("ZZMassPreFSR_10GeV", "ZZMassPreFSR_10GeV", 93, 20., 1000., ZZMassPreFsr, self.weight)

            GenZZMass = self.event.GenZZ_mass
            self.hEvent.fill1D("GenZZMass_10GeV", "GenZZMass_10GeV", 93, 20., 1000., GenZZMass, self.weight)
          
            # GenZZ: se guardo dentro ZZRo4l.root GenZZ ha indici, ma se li richiamo così vengono tutti 0; stesso problema per massa
            
            idx1 = self.event.GenZZ_Z1l1Idx
            self.hEvent.fill1D("GenZZidx1_10GeV", "GenZZIdx1_10GeV", 93, -100., 100., idx1, self.weight)
          
            # signal definition

            def ZZGammaSignalDefinition():
                if not self.analyzeMC: return False

                # leptons kinematic requirements
                if len(GenLeptons) != 4: return False
                if sum(l.pt > 5 for l in GenLeptons) < 4: return False
                if sum(l.pt > 10 for l in GenLeptons) < 2: return False
                if sum(l.pt > 20 for l in GenLeptons) < 1: return False
                if any(abs(l.eta) > 2.5 for l in GenLeptons): return False

                # photons kinematic requirments
                if len(GenPhotons) < 1: return False
                if sum(p.pt > 20 for p in GenPhotons) < 1: return False
                if any(abs(p.eta) > 2.4 or 1.444 < abs(p.eta) < 1.566 for p in GenPhotons): return False

            # plots

            def GetGenZZMass():
                m = self.event.GenZZ_mass
                self.hEvent.fill1D("GenZZMass_10GeV", "GenZZMass_10GeV", 93, 20., 1000., m, self.weight)
                return m

            def GetllGammaMassMin(ph):
                mllGamma1 = (GenParts[self.event.GenZZ_Z1l1Idx].p4 + GenParts[self.event.GenZZ_Z1l2Idx].p4 + GenParts[ph.genFsrIdx].p4).mass
                mll1 = (GenParts[self.event.GenZZ_Z1l1Idx].p4 + GenParts[self.event.GenZZ_Z1l2Idx].p4).mass
                mllGamma2 = (GenParts[self.event.GenZZ_Z2l1Idx].p4 + GenParts[self.event.GenZZ_Z2l2Idx].p4 + GenParts[ph.genFsrIdx].p4).mass
                mll2 = (GenParts[self.event.GenZZ_Z2l1Idx].p4 + GenParts[self.event.GenZZ_Z2l2Idx].p4).mass
                mllGammaMin = min(mllGamma1, mllGamma2)
                self.hEvent.fill1D("llGammaMassMin_10GeV", "llGammaMassMin_10GeV", 93, 20., 1000., mllGammaMin, self.weight)
                if mllGammaMin == mllGamma1: mllMin = mll1
                else: mllMin = mll2
                self.hEvent.fill2D("llGammaMassMin2D_10GeV", "llGammaMassMin2D_10GeV", 93, 20., 1000., 93, 20., 1000., mllGammaMin, mllMin, self.weight)
                return mllGammaMin

            #GenZZMass = GetGenZZMass()  !! len(GenPart) = 0
           
            #for ph in FsrPhotons:  !! len(GenPart) = 0
             #   mllGammaMin = GetllGammaMassMin(ph)
