from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection
from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Object
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

            # associate genlepton and lepton to fsr photon 

            def GetFsrAssociatedGenLepton(ph):
                if not self.analyzeMC: return None
                if ph.genFsrIdx != -1:
                    GenFsrPhoton = GenParts[ph.genFsrIdx]
                    GenFsrLpt = GenParts[GenFsrPhoton.genPartIdxMother]
                    if abs(GenFsrLpt.pdgId) == 11 or abs(GenFsrLpt.pdgId) == 13:
                        return GenFsrLpt
                return None

            def GetFsrAssociatedLepton(ph):
                dR = 0.5
                AssociatedLepton = None
                for l in Leptons:
                    dr = ph.DeltaR(l)
                    if dr < dR:
                        dR = dr
                        AssociatedLepton = l
                return AssociatedLepton    

            # plots

            def GetGenZZMass():
                m = self.event.GenZZ_mass
                self.hEvent.fill1D("GenZZMass_10GeV", "GenZZMass_10GeV", 93, 20., 1000., m, self.weight)
                return m

            def GetllGammaMassMin(ph):
                mllGamma1 = (GenParts[self.event.GenZZ_Z1l1Idx].p4() + GenParts[self.event.GenZZ_Z1l2Idx].p4() + ph.p4()).M()
                mll1 = (GenParts[self.event.GenZZ_Z1l1Idx].p4() + GenParts[self.event.GenZZ_Z1l2Idx].p4()).M()
                mllGamma2 = (GenParts[self.event.GenZZ_Z2l1Idx].p4() + GenParts[self.event.GenZZ_Z2l2Idx].p4() + ph.p4()).M()
                mll2 = (GenParts[self.event.GenZZ_Z2l1Idx].p4() + GenParts[self.event.GenZZ_Z2l2Idx].p4()).M()
                mllGammaMin = min(mllGamma1, mllGamma2)
                self.hEvent.fill1D("llGammaMassMin_10GeV", "llGammaMassMin_10GeV", 93, 20., 300., mllGammaMin, self.weight)
                if mllGammaMin == mllGamma1: mllMin = mll1
                else: mllMin = mll2
                self.hEvent.fill2D("llGammaMassMin2D_10GeV", "llGammaMassMin2D_10GeV", 93, 20., 500., 93, 20., 500., mllGammaMin, mllMin, self.weight)
                return mllGammaMin
            
            # ZZMass pre and post FSR          
            
            ZZMass = theZZ.mass
            self.hEvent.fill1D("ZZMass_10GeV", "ZZMass_10GeV", 93, 20., 1000., ZZMass, self.weight)

            ZZMassPreFsr = theZZ.massPreFSR
            self.hEvent.fill1D("ZZMassPreFSR_10GeV", "ZZMassPreFSR_10GeV", 93, 20., 1000., ZZMassPreFsr, self.weight)

            GenZZMass = GetGenZZMass()

            GenZZ4lMass = (GenParts[self.event.GenZZ_Z1l1Idx].p4() + GenParts[self.event.GenZZ_Z1l2Idx].p4() + GenParts[self.event.GenZZ_Z2l1Idx].p4() + GenParts[self.event.GenZZ_Z2l2Idx].p4()).M()
            self.hEvent.fill1D("GenZZ4lMass_10GeV", "GenZZ4lMass_10GeV", 93, 20., 1000., GenZZ4lMass, self.weight)

            # signal definition

            def SignalDefinition(RequireResonantZ2 = None, RequireThreeBosonregion = None):
                if not self.analyzeMC: return False

                # leptons kinematic requirements
                if len(GenLeptons) != 4: return False
                if (sum(l.pdgId == 11 for l in GenLeptons) != sum(l.pdgId == -11 for l in GenLeptons) or sum(l.pdgId == 13 for l in GenLeptons) != sum(l.pdgId == -13 for l in GenLeptons)): return False
                if sum(l.pt > 5 for l in GenLeptons) < 4: return False
                if sum(l.pt > 10 for l in GenLeptons) < 2: return False
                if sum(l.pt > 20 for l in GenLeptons) < 1: return False
                if any(abs(l.eta) > 2.5 for l in GenLeptons): return False

                # photons kinematic requirments
                if len(GenPhotons) < 1: return False
                if not any(p.pt > 20 and abs(p.eta) < 2.4 and not 1.444 < abs(p.eta) < 1.566 for p in GenPhotons): return False

                # Z1 mass
                Z1Mass = (GenParts[self.event.GenZZ_Z1l1Idx].p4() + GenParts[self.event.GenZZ_Z1l2Idx].p4()).M()
                if not 60 < Z1Mass < 120: return False
                self.hEvent.fill1D("Z1Mass_10GeV", "Z1Mass_10GeV", 93, 60., 120., Z1Mass, self.weight)
                
                if RequireResonantZ2 is not None:
                    if ResonantZ2() != RequireResonantZ2:
                        return False
                
                if RequireThreeBosonregion is not None:
                    if ThreeBosonRegion() != RequireThreeBosonregion:
                        return False

                return True

            def ResonantZ2():
                Z2Mass = (GenParts[self.event.GenZZ_Z2l1Idx].p4() + GenParts[self.event.GenZZ_Z2l2Idx].p4()).M()
                self.hEvent.fill1D("Z2Mass_10GeV", "Z2Mass_10GeV", 93, 20., 120., Z2Mass, self.weight)
                if 60 < Z2Mass < 120: return True
                if 20 < Z2Mass < 120: return False
                return None
                
            def ThreeBosonRegion():
                return any(
                    GetllGammaMassMin(ph) > 100 and
                    ph.DeltaR(GenParts[self.event.GenZZ_Z1l1Idx]) > 0.5 and
                    ph.DeltaR(GenParts[self.event.GenZZ_Z1l2Idx]) > 0.5 and
                    ph.DeltaR(GenParts[self.event.GenZZ_Z2l1Idx]) > 0.5 and
                    ph.DeltaR(GenParts[self.event.GenZZ_Z2l2Idx]) > 0.5
                    for ph in GenPhotons
                )
                    
                
            ZZGamma = SignalDefinition(True, True)
            Fsr = SignalDefinition(True, False)
            Higgs = SignalDefinition(False, True) # da controllare quale processo è
            quarto = SignalDefinition(False, False) # da controllare quale processo è