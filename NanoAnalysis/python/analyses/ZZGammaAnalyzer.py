from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection
from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Object
from ZZAnalysis.NanoAnalysis.tools import getLeptons, get_genEventSumw

from VVXAnalysis.NanoAnalysis.EventAnalyzer import EventAnalyzer
from VVXAnalysis.NanoAnalysis.Histogrammer import *
from VVXAnalysis.NanoAnalysis.Regions import Flags as Regions

import ROOT
from ROOT import TDatabasePDG

class ZZGammaAnalyzer(EventAnalyzer, analysis_name="ZZGammaAnalyzer"):

    def __init__(self, regions):
        super().__init__(regions)

    def analyze(self):
        bestCandIdx = self.event.bestCandIdx
    
        if(bestCandIdx != -1 and self.event.HLT_passZZ4l): 
            weight = 1.
            
            # Collections

            GenParts = Collection(self.event, 'GenPart')
            GenLeptons = [p for p in GenParts if abs(p.pdgId) == 11 or abs(p.pdgId) == 13 and p.status == 1]
            GenPhotons = [p for p in GenParts if p.pdgId == 22 and p.status == 1]
            self.hEvent.fill1D("nGenLeptons", "nGenLeptons", 93, 0., 10., len(GenLeptons), self.weight)
                        
            ZZs = Collection(self.event, 'ZZCand')
            theZZ = ZZs[bestCandIdx]
            Leptons = Collection(self.event, 'Lepton')
            Photons = Collection(self.event, 'Photon')
            FsrPhotons = Collection(self.event, 'FsrPhoton')

            if self.analyzeMC: self.weight = (self.event.overallEventWeight*theZZ.dataMCWeight/self.genEventSumw)

            # Signal Definition

            class GenZ:
                def __init__(self, Leptons):
                    self.leptons = Leptons
                    self.p4 = Leptons[0].p4() + Leptons[1].p4()
                    self.mass = self.p4.M()
                    self.pt = self.p4.Pt()
                    self.eta = self.p4.Eta()
                    self.phi = self.p4.Phi()
            
            class GenZZ:
                def __init__(self, Z1, Z2):
                    self.Z1 = Z1
                    self.Z2 = Z2
                    self.p4 = Z1.p4 + Z2.p4
                    self.mass = self.p4.M()
                    self.pt = self.p4.Pt()
                    self.eta = self.p4.Eta()
                    self.phi = self.p4.Phi()

            def GetZZ(Leptons):
                ZMassPDG = TDatabasePDG.Instance().GetParticle(23).Mass()
                Z1Mass = float("inf")
                Z1Leps = None
                for i, l1 in enumerate(Leptons):
                    for j in range(i+1, len(Leptons)):
                        l2 = Leptons[j]
                        if l1.pdgId != -l2.pdgId: continue
                        LeptonsMass = (l1.p4() + l2.p4()).M()
                        if abs(LeptonsMass - ZMassPDG) < abs(Z1Mass - ZMassPDG):
                            Z1Mass = LeptonsMass
                            Z1Leptons = [l1, l2]
                if Z1Leptons is None: return None
                Z2Leptons = [l for l in Leptons if l not in Z1Leptons]
                
                GenZ1Cand = GenZ(Z1Leptons)
                GenZ2Cand = GenZ(Z2Leptons)
                GenZZCand = GenZZ(GenZ1Cand, GenZ2Cand)

                return GenZZCand

            def GetZGammaMassMin(ph, GenZ1, GenZ2):
                Z1GammaMass = (GenZ1.leptons[0].p4() + GenZ1.leptons[1].p4() + ph.p4()).M()
                Z1Mass = GenZ1.mass
                Z2GammaMass = (GenZ2.leptons[0].p4() + GenZ2.leptons[1].p4() + ph.p4()).M()
                Z2Mass = GenZ2.mass
                ZGammaMassMin = min(Z1GammaMass, Z2GammaMass)
                self.hEvent.fill1D("ZGammaMassMin", "ZGammaMassMin", 93, 60., 300., ZGammaMassMin, self.weight)
                if ZGammaMassMin == Z1GammaMass: ZMassMin = Z1Mass
                else: ZMassMin = Z2Mass
                self.hEvent.fill2D("ZGammaMassMin2D", "ZGammaMassMin2D", 93, 60., 120., 93, 60., 300., ZMassMin, ZGammaMassMin, self.weight)
                return ZGammaMassMin

            def GetBestGamma(GoodPhotons):
                BestGamma = max(GoodPhotons, key=lambda ph: ph.pt, default=None)
                if BestGamma is None: return None
                BestGammaPt = BestGamma.pt
                self.hEvent.fill1D("BestGammaPt", "BestGammaPt", 50, 20., 220., BestGammaPt, self.weight)
                Gamma2 = max((ph for ph in GoodPhotons if ph != BestGamma), key=lambda ph: ph.pt, default=None)
                if Gamma2 is not None:
                    Gamma2Pt = Gamma2.pt
                    self.hEvent.fill1D("Gamma2Pt", "Gamma2Pt", 50, 20., 120., Gamma2Pt, self.weight)
                return BestGamma

            def SignalDefinition(GenLeptons, GenPhotons):

                if not self.analyzeMC: return False

                if len(GenLeptons) != 4: return False
                if (sum(l.pdgId == 11 for l in GenLeptons) != sum(l.pdgId == -11 for l in GenLeptons) or sum(l.pdgId == 13 for l in GenLeptons) != sum(l.pdgId == -13 for l in GenLeptons)): return False
                if sum(l.pt > 5 for l in GenLeptons) < 4: return False
                if sum(l.pt > 10 for l in GenLeptons) < 2: return False
                if sum(l.pt > 20 for l in GenLeptons) < 1: return False
                if any(abs(l.eta) > 2.5 for l in GenLeptons): return False

                GenZZ = GetZZ(GenLeptons)
                if GenZZ is None: return False
                if not 60 < GenZZ.Z1.mass < 120: return False
                if not 60 < GenZZ.Z2.mass < 120: return False
                self.hEvent.fill1D("GenZ1Mass", "GenZ1Mass", 100, 60, 120, GenZZ.Z1.mass, self.weight)
                self.hEvent.fill1D("GenZ2Mass", "GenZ2Mass", 100, 60, 120, GenZZ.Z2.mass, self.weight)
                                
                if len(GenPhotons) < 1: return False
                GoodGenPhotons = [ph for ph in GenPhotons if ph.pt > 20 and abs(ph.eta) < 2.4 and not 1.444 < abs(ph.eta) < 1.566]
                self.hEvent.fill1D("nGoodGenPhotons", "nGoodGenPhotons",5,-0.5, 4.5, len(GoodGenPhotons), self.weight)
                if len(GoodGenPhotons) == 0: return False
                GenGammaCands = [ph for ph in GoodGenPhotons if (ph.DeltaR(l) > 0.5 for l in GenLeptons) and (GetZGammaMassMin(ph, GenZZ.Z1, GenZZ.Z2) > 100)]
                self.hEvent.fill1D("nGenGammaCands", "nGenGammaCands",5,-0.5, 4.5, len(GenGammaCands), self.weight)
                if len(GenGammaCands) < 1: return False
                GenGamma = GetBestGamma(GenGammaCands)

                return True
                
            ZZGammaSignal = SignalDefinition(GenLeptons, GenPhotons)