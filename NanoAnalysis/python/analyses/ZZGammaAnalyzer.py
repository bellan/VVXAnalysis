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
    
        # Check that the event contains a selected candidate, and that
        # passes the required triggers (which is necessary for samples
        # processed with TRIGPASSTHROUGH=True)
        if(bestCandIdx != -1 and self.event.HLT_passZZ4l): 
            weight = 1.
            
            # Collections

            GenParts = Collection(self.event, 'GenPart')
            GenLeptons = [p for p in GenParts if abs(p.pdgId) == 11 or abs(p.pdgId) == 13]
            GenPhotons = [p for p in GenParts if p.pdgId == 22]

            ZZs = Collection(self.event, 'ZZCand') ## move it in EventAnalyzer::init(event) ??
            theZZ = ZZs[bestCandIdx]
            Leptons = Collection(self.event, 'Lepton')
            Photons = Collection(self.event, 'Photon')
            FsrPhotons = Collection(self.event, 'FsrPhoton')

            if self.analyzeMC: self.weight = (self.event.overallEventWeight*theZZ.dataMCWeight/self.genEventSumw)

            # associate leptons + fsrphotons

            def AssociateLeptonFsr(Leps, Phs):
                Associations = []
                for ph in Phs:
                    AssociatedLepton = None
                    dR = 0.5
                    for l in Leps:
                        dr = ph.DeltaR(l)
                        if dr < dR:
                            dR = dr
                            AssociatedLepton = l
                    if AssociatedLepton is not None:
                        Associations.append((AssociatedLepton, ph))
                return Associations

            def GetFsrsFromLepton(Lep, Associations):
                FsrPhotons = []
                for l, ph in Associations:
                    if Lep == l:
                        FsrPhotons.append(ph)
                return FsrPhotons

            def GetLeptonFromFsr(Photon, Associations):
                for l, ph in Associations:
                    if Photon == ph:
                        return l
                return None

            # ZZGen

            class GnZ:
                def __init__(self, leptons, fsr=None):
                    self.leptons = leptons
                    self.fsr = fsr

                    self.p4PreFsr = leptons[0].p4() + leptons[1].p4()
                    self.massPreFsr = self.p4PreFsr.M()

                    self.p4 = leptons[0].p4() + leptons[1].p4()
                    if fsr is not None:
                        for photons in fsr:
                            for ph in photons:
                                self.p4 += ph.p4()
                    self.mass = self.p4.M()
                    self.pt = self.p4.Pt()
                    self.eta = self.p4.Eta()
                    self.phi = self.p4.Phi()
                
            class GnZZ:
                def __init__(self, Z1, Z2):
                    self.Z1 = Z1
                    self.Z2 = Z2

                    self.p4PreFsr = Z1.p4PreFsr +Z2.p4PreFsr
                    self.massPreFsr = self.p4PreFsr.M()
                    self.p4 = Z1.p4 +Z2.p4
                    self.mass = self.p4.M()
                    self.pt = self.p4.Pt()
                    self.eta = self.p4.Eta()
                    self.phi = self.p4.Phi()
                    
            def GetGnZZ(GenLeps, GenFsrs):
                if len(GenLeps) != 4: return None
                if (sum(l.pdgId == 11 for l in GenLeps) != sum(l.pdgId == -11 for l in GenLeps) or sum(l.pdgId == 13 for l in GenLeps) != sum(l.pdgId == -13 for l in GenLeps)): return None
                if sum(l.pt > 5 for l in GenLeps) < 4: return None
                if sum(l.pt > 10 for l in GenLeps) < 2: return None
                if sum(l.pt > 20 for l in GenLeps) < 1: return None
                if any(abs(l.eta) > 2.5 for l in GenLeps): return None

                ZMassPDG = TDatabasePDG.Instance().GetParticle(23).Mass()
                Z1MassPreFsr = float("inf")
                Z1Leps = None
                for i, l1 in enumerate(GenLeps):
                    for j in range(i+1, len(GenLeps)):
                        l2 = GenLeps[j]
                        if l1.pdgId != -l2.pdgId: continue
                        LepsMass = (l1.p4() + l2.p4()).M()
                        if abs(LepsMass - ZMassPDG) < abs(Z1MassPreFsr - ZMassPDG):
                            Z1MassPreFsr = LepsMass
                            Z1Leps = [l1, l2]
                if Z1Leps is None: return None
                
                Z2Leps = [l for l in GenLeps if l not in Z1Leps]
                
                AssociationsLptFsr = AssociateLeptonFsr(GenLeps, GenFsrs)
                Z1Fsrs = [GetFsrsFromLepton(Z1Leps[0], AssociationsLptFsr), GetFsrsFromLepton(Z1Leps[1], AssociationsLptFsr)]
                Z2Fsrs = [GetFsrsFromLepton(Z2Leps[0], AssociationsLptFsr), GetFsrsFromLepton(Z2Leps[1], AssociationsLptFsr)]
                         
                GnZ1 = GnZ(Z1Leps, Z1Fsrs)
                GnZ2 = GnZ(Z2Leps, Z2Fsrs)

                return GnZZ(GnZ1, GnZ2)

            GnZZcand = GetGnZZ(GenLeptons, GenPhotons)

            #signal definition

            def GetllGammaMassMin(ph, GnZ1, GnZ2):
                mllGamma1 = (GnZ1.leptons[0].p4() + GnZ1.leptons[1].p4() + ph.p4()).M()
                mll1 = (GnZ1.leptons[0].p4() + GnZ1.leptons[1].p4()).M()
                mllGamma2 = (GnZ2.leptons[0].p4() + GnZ2.leptons[1].p4() + ph.p4()).M()
                mll2 = (GnZ2.leptons[0].p4() + GnZ2.leptons[1].p4()).M()
                mllGammaMin = min(mllGamma1, mllGamma2)
                self.hEvent.fill1D("llGammaMassMin", "llGammaMassMin", 93, 60., 120., mllGammaMin, self.weight)
                if mllGammaMin == mllGamma1: mllMin = mll1
                else: mllMin = mll2
                self.hEvent.fill2D("llGammaMassMin2D", "llGammaMassMin2D", 93, 60., 120., 93, 20., 200., mllGammaMin, mllMin, self.weight)
                return mllGammaMin

            def GetBestGamma(GoodPhotons):
                BestGamma = max(GoodPhotons, key=lambda ph: ph.pt, default=None)
                if BestGamma is None: return None
                BestGammapt = BestGamma.pt
                self.hEvent.fill1D("BestGammapt", "BestGammapt", 93, 20., 1000., BestGammapt, self.weight)
                Gamma2 = max((ph for ph in GoodPhotons if ph != BestGamma), key=lambda ph: ph.pt, default=None)
                if Gamma2 is not None:
                    Gamma2pt = Gamma2.pt
                    self.hEvent.fill1D("Gamma2pt", "Gamma2pt", 93, 20., 1000., Gamma2pt, self.weight)
                return BestGamma
                
            def ResonantZ2(GnZ2):
                if 60 < GnZ2.mass < 120: return True
                return False
                
            def ThreeBosonRegion(Gamma, GnZ1, GnZ2):
                if any(Gamma.DeltaR(l) < 0.5 for l in GnZ1.leptons + GnZ2.leptons): return False
                if GetllGammaMassMin(Gamma, GnZ1, GnZ2) < 100: return False
                return True

            def SignalDefinition(RequireResonantZ2 = None, RequireThreeBosonregion = None):
                if not self.analyzeMC: return False

                if GnZZcand is None: return False
                GnZ1 = GnZZcand.Z1
                GnZ2 = GnZZcand.Z2

                if not 60 < GnZ1.mass < 120: return False
                if not 12 < GnZ2.mass < 120: return False

                if len(GenPhotons) < 1: return False
                GoodPhotons = [p for p in GenPhotons if p.pt > 20 and abs(p.eta) < 2.4 and not 1.444 < abs(p.eta) < 1.566]
                self.hEvent.fill1D("nGoodPhotons", "nGoodPhotons", 10, 0, 10, len(GoodPhotons), self.weight)
                if len(GoodPhotons) == 0: return False
                
                Gamma = GetBestGamma(GoodPhotons)

                if RequireResonantZ2 is not None:
                    if ResonantZ2(GnZ2) != RequireResonantZ2:
                        return False
                
                if RequireThreeBosonregion is not None:
                    if ThreeBosonRegion(Gamma, GnZ1, GnZ2) != RequireThreeBosonregion:
                        return False

                return True

            ZZGamma = SignalDefinition(True, True)
            ZZFsr = SignalDefinition(True, False)
            HiggsGamma = SignalDefinition(False, True)
            HiggsFsr = SignalDefinition(False, False)

            # ZZMass pre and post FSR          
            if GnZZcand is not None:
                ZZMass = theZZ.mass
                self.hEvent.fill1D("ZZMass_10GeV", "ZZMass_10GeV", 93, 20., 1000., ZZMass, self.weight)
    
                ZZMassPreFsr = theZZ.massPreFSR
                self.hEvent.fill1D("ZZMassPreFSR_10GeV", "ZZMassPreFSR_10GeV", 93, 20., 1000., ZZMassPreFsr, self.weight)

                GnZZMass = GnZZcand.mass    
                self.hEvent.fill1D("GnZZMass", "GnZZMass", 93, 20., 1000., GnZZMass, self.weight)

                GnZZMassPreFsr = GnZZcand.massPreFsr
                self.hEvent.fill1D("GnZZMassPreFSR", "GnZZMassPreFSR", 93, 20., 1000., GnZZMassPreFsr, self.weight)
