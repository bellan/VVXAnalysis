from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection
from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Object
from ZZAnalysis.NanoAnalysis.tools import getLeptons, get_genEventSumw

from VVXAnalysis.NanoAnalysis.EventAnalyzer import EventAnalyzer
from VVXAnalysis.NanoAnalysis.Histogrammer import *
from VVXAnalysis.NanoAnalysis.Regions import Flags as Regions

from enum import IntFlag

import ROOT
from ROOT import TDatabasePDG
from ROOT import TLorentzVector

class ZZGammaAnalyzer(EventAnalyzer, analysis_name="ZZGammaAnalyzer"):

    def __init__(self, regions):
        super().__init__(regions)
        
    def analyze(self):

        self.LoadCollections()
        
        self.CutPlot("Definition", 0)
        self.CutPlot("Region", 0)
        
        isSD = self.SignalDefinition(self.GenLeptons, self.GenPhotons)
        
        bestCandIdx = self.event.bestCandIdx
        if(bestCandIdx != -1 and self.event.HLT_passZZ4l): 
            theZZ = self.ZZs[bestCandIdx]
            self.CutPlot("Region", 1)
            isSR = self.SignalRegion(theZZ, self.Leptons, self.Photons)

        #if isSD:
        #    if isSR:
        #        ResultSDSR = 1
        #    else:
        #        ResultSDSR = 0
        #    self.hEvent.fill1D("CfrSignalDefSignalRegion", "cfrSignalDefSignalRegion", 2, -0.5, 1.5, ResultSDSR, self.weight)
                

    # =========================
    # SIGNAL DEF AND REGION
    # =========================

    def SignalDefinition(self, GenLeptons, GenPhotons):

        if not self.analyzeMC: return False

        if not self.LeptonsKinematicCuts(GenLeptons, "Definition"): return False
        self.CutPlot("Definition", 1)
        GenZZ = self.GetZZ(GenLeptons)
        if not self.ZZMassCuts(GenZZ): return False
        self.CutPlot("Definition", 2)
        
        KinPassed, GoodPhotons = self.PhotonsKinematicCuts(GenPhotons)
        if not KinPassed: return False
        self.CutPlot("Definition", 3)
        FsrPassed, GenGammaCands = self.PhotonsFsrCuts(GoodPhotons, GenLeptons, GenZZ)
        if not FsrPassed: return False
        self.CutPlot("Definition", 4)
        GenGamma = self.GetBestGamma(GenGammaCands, "Gen")
        if GenGamma is None: return False
        PtThetaGraph = self.GraphGammaPtThetaZZG(GenGamma, GenZZ, "Gen")

        return True

    SDCuts = ["EventiTotali", "CinematicaLeptoni", "MasseZ1Z2", "CinematicaFotoni", "NotFsrPhotons"]

    def SignalRegion(self, theZZ, Leptons, Photons):
                
        if not Regions.check(self.regionWord, Regions.L4P): return False
        self.CutPlot("Region", 2)
        if not self.LeptonsKinematicCuts(Leptons, "Region"): return False
        if not self.ZZMassCuts(theZZ): return False
        self.CutPlot("Region", 8)
        
        if not Regions.check(self.regionWord, Regions.P1cutL): return False
        self.CutPlot("Region", 9)
        KinPassed, GoodPhotons = self.PhotonsKinematicCuts(Photons)
        if not KinPassed: return False
        self.CutPlot("Region", 10)
        FsrPassed, RecoGammaCands = self.PhotonsFsrCuts(GoodPhotons, Leptons, theZZ)
        if not FsrPassed: return False
        self.CutPlot("Region", 11)
        RecoGamma = self.GetBestGamma(RecoGammaCands, "Reco")
        if RecoGamma is None: return False
        PtThetaGraph = self.GraphGammaPtThetaZZG(RecoGamma, theZZ, "Reco")

        return True

    SRCuts = ["EventiTotali", "passZZ4l+BestCand", "regione L4P", "4Leptoni", "pt5", "pt10", "pt20", "pteta", "Z1Z2Mass", "regioneP1cutL", "CinematicaFotoni", "NotFsrPhotons"]

    def CutPlot(self, Level, Cut):
        if Level == "Definition":
            Labels = self.SDCuts
        elif Level == "Region":
            Labels = self.SRCuts
        self.hEvent.fill1D_label(f"EventiPerTaglio_Signal{Level}", f"EventiPerTaglio_Signal{Level}", Labels, Labels[Cut], self.weight)
        

    # =========================
    # COLLECTIONS
    # =========================

    def LoadCollections(self):
        if self.analyzeMC:
            self.GenParts = Collection(self.event, 'GenPart')
            self.GenLeptons = [p for p in self.GenParts if abs(p.pdgId) in (11, 13)
                                                 and p.status == 1
                                                 and (p.statusFlags & GenStatusFlag.SEL) == GenStatusFlag.SEL]
            self.GenPhotons = [p for p in self.GenParts if p.pdgId == 22
                                                 and p.status == 1
                                                 and (p.statusFlags & GenStatusFlag.SEL) == GenStatusFlag.SEL]
        self.ZZs = Collection(self.event, 'ZZCand')
        self.Leptons = Collection(self.event, 'Lepton')
        self.Photons = Collection(self.event, 'Photon')
            

    # =========================
    # USEFUL FUNCTIONS: ZZ
    # =========================

    def LeptonsKinematicCuts(self, Leptons, Level):
        #if len(Leptons) != 4: return False
        #if Level == "Region": self.CutPlot("Region", 3)
        if Level == "Definition":
            if (sum(l.pdgId == 11 for l in Leptons) != sum(l.pdgId == -11 for l in Leptons)
                or sum(l.pdgId == 13 for l in Leptons) != sum(l.pdgId == -13 for l in Leptons)): return False
        if sum(l.pt > 5 for l in Leptons) < 4: return False
        if Level == "Region": self.CutPlot("Region", 4)
        if sum(l.pt > 10 for l in Leptons) < 2: return False
        if Level == "Region": self.CutPlot("Region", 5)
        if sum(l.pt > 20 for l in Leptons) < 1: return False
        if Level == "Region": self.CutPlot("Region", 6)
        if any(abs(l.eta) > 2.5 for l in Leptons): return False
        if Level == "Region": self.CutPlot("Region", 7)
        return True
    
    def GetZZ(self, Leptons):
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
        if (len(Z1Leptons) != 2
            or Z1Leptons[0].pdgId != - Z1Leptons[1].pdgId): return None
        Z2Leptons = [l for l in Leptons if l not in Z1Leptons]
        if (len(Z2Leptons) != 2
            or Z2Leptons[0].pdgId != - Z2Leptons[1].pdgId): return None
                
        GenZ1Cand = ZClass(Z1Leptons)
        GenZ2Cand = ZClass(Z2Leptons)
        GenZZCand = ZZClass(GenZ1Cand, GenZ2Cand)

        return GenZZCand

    def ZZMassCuts(self, ZZ):
        if ZZ is None: return False
        if hasattr(ZZ, "Z1"):
            Z1Mass = ZZ.Z1.mass
            Z2Mass = ZZ.Z2.mass
            prefix = "Gen"
        elif hasattr(ZZ, "Z1mass"):
            Z1Mass = ZZ.Z1mass
            Z2Mass = ZZ.Z2mass
            prefix = "Reco"
        if not 60 < Z1Mass < 120: return False
        if not 60 < Z2Mass < 120: return False
        self.hEvent.fill1D(f"{prefix}Level_Z1Mass", f"{prefix}Level_Z1Mass", 100, 60., 120., Z1Mass, self.weight)
        self.hEvent.fill1D(f"{prefix}Level_Z2Mass", f"{prefix}Level_Z2Mass", 75, 60., 120., Z2Mass, self.weight)
        return True

    
    # =========================
    # USEFUL FUNCTIONS: PHOTONS
    # =========================

    def GetP4(self, ph):
        php4 = TLorentzVector()
        php4.SetPtEtaPhiM(ph.pt, ph.eta, ph.phi, 0.)
        return php4

    def PhotonsKinematicCuts(self, Photons):
        if len(Photons) < 1: return False, None
        GoodPhotons = [ph for ph in Photons if ph.pt > 20
                                               and abs(ph.eta) < 2.4
                                               and not 1.444 < abs(ph.eta) < 1.566]
        if len(GoodPhotons) == 0: return False, None
        return True, GoodPhotons

    def PhotonsFsrCuts(self, Photons, Leptons, ZZ):
        if (hasattr(ZZ, "Z1") and hasattr(ZZ, "Z2")):
            GammaCands = [ph for ph in Photons if all(ph.DeltaR(l) > 0.5 for l in Leptons)
                                                  and self.GetllGammaMassMin(ph, ZZ) > 100]
        elif (hasattr(ZZ, "Z1l1Idx") and hasattr(ZZ, "Z1l2Idx") and hasattr(ZZ, "Z2l1Idx") and hasattr(ZZ, "Z2l2Idx")):
            ZZRecoLeptons = [Leptons[ZZ.Z1l1Idx], Leptons[ZZ.Z1l2Idx], Leptons[ZZ.Z2l1Idx], Leptons[ZZ.Z2l2Idx]]
            GammaCands = [ph for ph in Photons if all(ph.DeltaR(l) > 0.5 for l in ZZRecoLeptons)
                                                  and self.GetllGammaMassMin(ph, ZZ, Leptons) > 100]
        if len(GammaCands) < 1: return False, None
        return True, GammaCands

    def GetllGammaMassMin(self, ph, ZZ, Leptons = None):
        if (hasattr(ZZ, "Z1") and hasattr(ZZ, "Z2")):
            ll1GammaMass = (ZZ.Z1.p4 + ph.p4()).M()
            ll1Mass = ZZ.Z1.mass
            ll2GammaMass = (ZZ.Z2.p4 + ph.p4()).M()
            ll2Mass = ZZ.Z2.mass
            prefix = "Gen"
        elif (hasattr(ZZ, "Z1l1Idx") and hasattr(ZZ, "Z1l2Idx") and hasattr(ZZ, "Z2l1Idx") and hasattr(ZZ, "Z2l2Idx")):
            ll1GammaMass = (Leptons[ZZ.Z1l1Idx].p4() + Leptons[ZZ.Z1l2Idx].p4() + self.GetP4(ph)).M()
            ll1Mass = (Leptons[ZZ.Z1l1Idx].p4() + Leptons[ZZ.Z1l2Idx].p4()).M()
            ll2GammaMass = (Leptons[ZZ.Z2l1Idx].p4() + Leptons[ZZ.Z2l2Idx].p4() + self.GetP4(ph)).M()
            ll2Mass = (Leptons[ZZ.Z2l1Idx].p4() + Leptons[ZZ.Z2l2Idx].p4()).M()
            prefix = "Reco"
        llGammaMassMin = min(ll1GammaMass, ll2GammaMass)
        self.hEvent.fill1D(f"{prefix}Level_llGammaMinMass", f"{prefix}Level_llGammaMinMass", 52, 40., 300., llGammaMassMin, self.weight)
        if llGammaMassMin == ll1GammaMass: llMassMin = ll1Mass
        else: llMassMin = ll2Mass
        self.hEvent.fill2D(f"{prefix}Level_llGammaMassMin2D", f"{prefix}Level_llMinMass(llGammaMinMass)", 104, 40., 300., 120, 60., 120., llGammaMassMin, llMassMin, self.weight)
        return llGammaMassMin

    def GetBestGamma(self, GammaCands, prefix):
        BestGamma = max(GammaCands, key=lambda ph: ph.pt, default=None)
        if BestGamma is None: return None
        BestGammaPt = BestGamma.pt
        self.hEvent.fill1D(f"{prefix}Level_BestGammaPt", f"{prefix}Level_BestGammaPt", 50, 20., 220., BestGammaPt, self.weight)
        Gamma2 = max((ph for ph in GammaCands if ph != BestGamma), key=lambda ph: ph.pt, default=None)
        if Gamma2 is not None:
            Gamma2Pt = Gamma2.pt
            self.hEvent.fill1D(f"{prefix}Level_Gamma2Pt", f"{prefix}Level_Gamma2Pt", 20, 20., 220., Gamma2Pt, self.weight)
        return BestGamma

    def GraphGammaPtThetaZZG(self, Gamma, ZZ, prefix):
        if prefix == "Gen":
            thetaZZGamma = ZZ.p4.Vect().Angle(self.GetP4(Gamma).Vect())
        elif prefix == "Reco":
            thetaZZGamma = ZZ.p4().Vect().Angle(self.GetP4(Gamma).Vect())
        self.hEvent.fill2D(f"{prefix}Level_BestGammaPt-AngoloZZGamma", f"{prefix}Level_BestGammaPt-AngoloZZGamma", 50, 0., ROOT.TMath.Pi(), 5, 60., 100., thetaZZGamma , Gamma.pt, self.weight)
        return None

# =========================
# USEFUL CLASSES
# =========================

class GenStatusFlag(IntFlag):
    IS_PROMPT = 1 << 0
    FROM_HARD_PROCESS = 1 << 8

    SEL = IS_PROMPT | FROM_HARD_PROCESS
    #PH_SEL = IS_PROMPT

class ZClass:
    def __init__(self, Leptons):
        self.leptons = Leptons
        self.p4 = Leptons[0].p4() + Leptons[1].p4()
        self.mass = self.p4.M()
        self.pt = self.p4.Pt()
        self.eta = self.p4.Eta()
        self.phi = self.p4.Phi()
            
class ZZClass:
    def __init__(self, Z1, Z2):
        self.Z1 = Z1
        self.Z2 = Z2
        self.p4 = Z1.p4 + Z2.p4
        self.mass = self.p4.M()
        self.pt = self.p4.Pt()
        self.eta = self.p4.Eta()
        self.phi = self.p4.Phi()