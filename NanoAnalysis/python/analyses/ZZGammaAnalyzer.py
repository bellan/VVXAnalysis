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
        
        isSD = self.SignalDefinition(self.GenLeptons, self.GenPhotons)
        isSR = self.SignalRegion(self.Leptons, self.Photons)

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
        self.CutPlot("Definition", 0)

        if not self.LeptonsKinematicCuts("Gen", GenLeptons): return False
        self.CutPlot("Definition", 1)
        GenZZ = self.GetZZ(GenLeptons)
        if not self.ZZMassCuts("Gen", GenZZ): return False
        self.CutPlot("Definition", 2)
        
        KinPass, GoodPhotons = self.PhotonsKinematicCuts(GenPhotons)
        if not KinPass: return False
        self.CutPlot("Definition", 3)
        FsrPass, GenGammaCands = self.PhotonsFsrCuts("Gen", GoodPhotons, GenLeptons, GenZZ)
        if not FsrPass: return False
        self.CutPlot("Definition", 4)
        GenGamma = self.GetBestGamma("Gen", GenGammaCands)
        PtThetaGraph = self.GraphGammaPtThetaZZG("Gen", GenGamma, GenZZ)

        return True

    SDCuts = ["EventiTotali", "CinematicaLeptoni", "MasseZ1Z2", "CinematicaFotoni", "TribosonRegion"]

    def SignalRegion(self, Leptons, Photons):

        self.CutPlot("Region", 0)
        
        bestCandIdx = self.event.bestCandIdx
        if not (bestCandIdx != -1 and self.event.HLT_passZZ4l): return False
        self.CutPlot("Region", 1)
        theZZ = self.ZZs[bestCandIdx]
        LeptonsOI = [Leptons[theZZ.Z1l1Idx], Leptons[theZZ.Z1l2Idx], Leptons[theZZ.Z2l1Idx], Leptons[theZZ.Z2l2Idx]]    
        if not Regions.check(self.regionWord, Regions.L4P): return False
        self.CutPlot("Region", 2)
        if not self.LeptonsKinematicCuts("Reco", LeptonsOI): return False
        if not self.ZZMassCuts("Reco", theZZ): return False
        self.CutPlot("Region", 9)
        
        if not Regions.check(self.regionWord, Regions.P1cutL): return False
        self.CutPlot("Region", 10)
        KinPass, GoodPhotons = self.PhotonsKinematicCuts(Photons)
        if not KinPass: return False
        self.CutPlot("Region", 11)
        FsrPass, RecoGammaCands = self.PhotonsFsrCuts("Reco", GoodPhotons, LeptonsOI, theZZ) #!!!!!!!!!!!!!!! leptonsOI !!!!!!!!!!!!!!
        if not FsrPass: return False
        self.CutPlot("Region", 12)
        RecoGamma = self.GetBestGamma("Reco", RecoGammaCands)
        if RecoGamma is None: return False
        PtThetaGraph = self.GraphGammaPtThetaZZG("Reco", RecoGamma, theZZ)

        return True

    SRCuts = ["EventiTotali", "passZZ4l+BestCand", "regione L4P", "4Leps", "pdgId", "pt5", "pt10", "pt20", "eta", "Z1Z2Mass", "regioneP1cutL", "CinematicaFotoni", "TribosonRegion"]

    
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

    def LeptonsKinematicCuts(self, Level, Leptons):
        if len(Leptons) != 4: return False
        if Level == "Reco": self.CutPlot("Region", 3)
        if (sum(l.pdgId == 11 for l in Leptons) != sum(l.pdgId == -11 for l in Leptons)
            or sum(l.pdgId == 13 for l in Leptons) != sum(l.pdgId == -13 for l in Leptons)): return False
        if Level == "Reco": self.CutPlot("Region", 4)
        if sum(l.pt >= self.LepsCuts(Level, l)[0] for l in Leptons) < 4: return False
        if Level == "Reco": self.CutPlot("Region", 5)
        if sum(l.pt >= 10 for l in Leptons) < 2: return False
        if Level == "Reco": self.CutPlot("Region", 6)
        if sum(l.pt >= 20 for l in Leptons) < 1: return False
        if Level == "Reco": self.CutPlot("Region", 7)
        if any(abs(l.eta) > self.LepsCuts(Level, l)[1] for l in Leptons): return False
        if Level == "Reco": self.CutPlot("Region", 8)
        
        return True

    def LepsCuts(self, Level, l):
        if Level == "Gen":
            pt = 5
            eta = 2.5
        elif Level == "Reco":
            if abs(l.pdgId) == 11:
                pt = 7
                eta = 2.5
            elif abs(l.pdgId) == 13:
                pt = 5
                eta = 2.4
        return pt, eta
    
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

    def ZZMassCuts(self, Level, ZZ):
        if ZZ is None: return False
        if Level == "Gen":
            Z1Mass = ZZ.Z1.mass
            Z2Mass = ZZ.Z2.mass
        elif Level == "Reco":
            Z1Mass = ZZ.Z1mass
            Z2Mass = ZZ.Z2mass
        if not 60 < Z1Mass < 120: return False
        if not 60 < Z2Mass < 120: return False
        self.hEvent.fill1D(f"{Level}Level_Z1Mass", f"{Level}Level_Z1Mass", 100, 60., 120., Z1Mass, self.weight)
        self.hEvent.fill1D(f"{Level}Level_Z2Mass", f"{Level}Level_Z2Mass", 75, 60., 120., Z2Mass, self.weight)
        return True

    
    # =========================
    # USEFUL FUNCTIONS: PHOTONS
    # =========================

    def PhotonsKinematicCuts(self, Photons):
        if len(Photons) < 1: return False, None
        GoodPhotons = [ph for ph in Photons if ph.pt > 20
                                               and abs(ph.eta) < 2.4
                                               and not 1.444 < abs(ph.eta) < 1.566]
        if len(GoodPhotons) == 0: return False, None
        return True, GoodPhotons

    def PhotonsFsrCuts(self, Level, Photons, Leptons, ZZ):
        GammaCands = [ph for ph in Photons if all(ph.DeltaR(l) > 0.5 for l in Leptons)
                                                  and self.GetllGammaMassMin(Level, ph, ZZ, Leptons) > 100]
        if len(GammaCands) < 1: return False, None
        return True, GammaCands

    def GetllGammaMassMin(self, Level, ph, ZZ, Leptons = None):
        Z1Mass = ZZ.Z1mass
        Z2Mass = ZZ.Z2mass
        if Level == "Gen":
            Z1GammaMass = (ZZ.Z1.p4 + ph.p4()).M()
            Z2GammaMass = (ZZ.Z2.p4 + ph.p4()).M()
        #elif Level == "Reco":
        #    Z1GammaMass = (self.GetP4())

        elif Level == "Reco":
            Z1GammaMass = (Leptons[0].p4() + Leptons[1].p4() + self.GetP4(ph)).M()
            Z1Mass = (Leptons[0].p4() + Leptons[1].p4()).M()
            Z2GammaMass = (Leptons[2].p4() + Leptons[3].p4() + self.GetP4(ph)).M()
            Z2Mass = (Leptons[2].p4() + Leptons[3].p4()).M()
        ZGammaMassMin = min(Z1GammaMass, Z2GammaMass)
        self.hEvent.fill1D(f"{Level}Level_ZGammaMinMass", f"{Level}Level_ZGammaMinMass", 52, 40., 300., ZGammaMassMin, self.weight)
        if ZGammaMassMin == Z1GammaMass: ZMassMin = Z1Mass
        else: ZMassMin = Z2Mass
        self.hEvent.fill2D(f"{Level}Level_llGammaMassMin2D", f"{Level}Level_llMinMass(llGammaMinMass)", 104, 40., 300., 120, 60., 120., ZGammaMassMin, ZMassMin, self.weight)
        return ZGammaMassMin

    def GetBestGamma(self, Level, GammaCands):
        BestGamma = max(GammaCands, key=lambda ph: ph.pt, default=None)
        if BestGamma is None: return None
        BestGammaPt = BestGamma.pt
        self.hEvent.fill1D(f"{Level}Level_BestGammaPt", f"{Level}Level_BestGammaPt", 50, 20., 220., BestGammaPt, self.weight)
        Gamma2 = max((ph for ph in GammaCands if ph != BestGamma), key=lambda ph: ph.pt, default=None)
        if Gamma2 is not None:
            Gamma2Pt = Gamma2.pt
            self.hEvent.fill1D(f"{Level}Level_Gamma2Pt", f"{Level}Level_Gamma2Pt", 20, 20., 220., Gamma2Pt, self.weight)
        return BestGamma

    def GraphGammaPtThetaZZG(self, Level, Gamma, ZZ):
        if Level == "Gen":
            thetaZZGamma = ZZ.p4.Vect().Angle(self.GetP4(Gamma).Vect())
        elif Level == "Reco":
            thetaZZGamma = ZZ.p4().Vect().Angle(self.GetP4(Gamma).Vect())
        self.hEvent.fill2D(f"{Level}Level_BestGammaPt-AngoloZZGamma", f"{Level}Level_BestGammaPt-AngoloZZGamma", 50, 0., ROOT.TMath.Pi(), 5, 60., 100., thetaZZGamma , Gamma.pt, self.weight)
        return None

    # =========================
    # AUXILIARY FUNCTIONS
    # =========================
    
    def GetP4(self, p):
        p4 = TLorentzVector()
        if hasattr(p, "mass"):
            p4.SetPtEtaPhiM(p.pt, p.eta, p.phi, p.mass)
        else:
            p4.SetPtEtaPhiM(p.pt, p.eta, p.phi, 0.)
        return p4

    def CutPlot(self, Level, Cut):
        if Level == "Definition":
            Labels = self.SDCuts
        elif Level == "Region":
            Labels = self.SRCuts
        self.hEvent.fill1D_label(f"EventiPerTaglio_Signal{Level}", f"EventiPerTaglio_Signal{Level}", Labels, Labels[Cut], self.weight)
        

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
        self.Z1mass = Z1.mass
        self.Z2mass = Z2.mass
        self.pt = self.p4.Pt()
        self.eta = self.p4.Eta()
        self.phi = self.p4.Phi()