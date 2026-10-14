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

        self.isBgivenA(isSD, isSR, "SD", "SR")
        self.isBgivenA(isSR, isSD, "SR", "SD")

        
    # =========================
    # SIGNAL DEF AND REGION
    # =========================

    def SignalDefinition(self, GenLeptons, GenPhotons):

        if not self.analyzeMC: return False
        self.CutPlot("Definition", 0)

        if not self.LeptonsKinematicCuts("Gen", GenLeptons): return False
        self.CutPlot("Definition", 1)
        GenZZ = self.GetZZ(GenLeptons)
        if GenZZ is None: return False
        if not self.ZZMassCuts(GenZZ): return False
        self.CutPlot("Definition", 2)
        
        KinPass, GoodPhotons = self.PhotonsKinematicCuts(GenPhotons)
        if not KinPass: return False
        self.CutPlot("Definition", 3)
        IsolatedPass, IsolatedPhotons = self.GetIsolatedPhotons(GoodPhotons, GenLeptons)
        if not IsolatedPass: return False
        self.CutPlot("Definition", 4)
        TribosonPass, GenGamma = self.GetGamma("Gen", IsolatedPhotons, GenLeptons, GenZZ)
        if not TribosonPass: return False
        self.CutPlot("Definition", 5)
        PtThetaGraph = self.GraphGammaPtThetaZZG("Gen", GenGamma, GenZZ)

        self.PlotsOI("Gen", GenZZ, GenGamma)

        return True

    SDCuts = ["EventiTotali", "CinematicaLeptoni", "MasseZ1Z2", "CinematicaFotoni", "IsolatedPhotons", "TribosonRegion"]

    def SignalRegion(self, Leptons, Photons):

        self.CutPlot("Region", 0)
        
        bestCandIdx = self.event.bestCandIdx
        if not (bestCandIdx != -1 and self.event.HLT_passZZ4l): return False
        self.CutPlot("Region", 1)
        RecoZZ = self.ZZs[bestCandIdx]
        LeptonsOI = [Leptons[RecoZZ.Z1l1Idx], Leptons[RecoZZ.Z1l2Idx], Leptons[RecoZZ.Z2l1Idx], Leptons[RecoZZ.Z2l2Idx]]    
        if not Regions.check(self.regionWord, Regions.L4P): return False
        self.CutPlot("Region", 2)
        if not self.LeptonsKinematicCuts("Reco", LeptonsOI): return False
        self.CutPlot("Region", 3)
        if not self.ZZMassCuts(RecoZZ): return False
        self.CutPlot("Region", 4)
        
        RegionPass, TruePhotons = self.GetTruePhotons(Photons)
        if not RegionPass: return False
        self.CutPlot("Region", 5)
        KinPass, GoodPhotons = self.PhotonsKinematicCuts(TruePhotons)
        if not KinPass: return False
        self.CutPlot("Region", 6)
        IsolatedPass, IsolatedPhotons = self.GetIsolatedPhotons(GoodPhotons, LeptonsOI)
        if not IsolatedPass: return False
        self.CutPlot("Region", 7)
        TribosonPass, RecoGamma = self.GetGamma("Reco", IsolatedPhotons, LeptonsOI, RecoZZ)
        if not TribosonPass: return False
        self.CutPlot("Region", 8)
        PtThetaGraph = self.GraphGammaPtThetaZZG("Reco", RecoGamma, RecoZZ)

        self.PlotsOI("Reco", RecoZZ, RecoGamma)

        return True

    SRCuts = ["EventiTotali", "passZZ4l+BestCand", "regione L4P", "CinematicaLeptoni", "Z1Z2Mass",  "P1CutL+ElectronVeto", "CinematicaFotoni", "IsolatedPhotons","TribosonRegion"]

    def isBgivenA(self, A, B, nomeA, nomeB):

        if A:
            if B:
                result = 1
            else: result = 0
            self.hEvent.fill1D(f"is{nomeB}given{nomeA}", f"is{nomeB}given{nomeA}", 2, -0.5, 1.5, result, self.weight)


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
        if (sum(l.pdgId == 11 for l in Leptons) != sum(l.pdgId == -11 for l in Leptons)
            or sum(l.pdgId == 13 for l in Leptons) != sum(l.pdgId == -13 for l in Leptons)): return False

        if sum(l.pt > self.LepsCuts(Level, l)[0] for l in Leptons) < 4: return False
        if sum(l.pt > 10 for l in Leptons) < 2: return False
        if sum(l.pt > 20 for l in Leptons) < 1: return False

        if any(abs(l.eta) >= self.LepsCuts(Level, l)[1] for l in Leptons): return False
        
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
        Electrons = [l for l in Leptons if abs(l.pdgId) == 11]
        Muons = [l for l in Leptons if abs(l.pdgId) == 13]
        
        if len(Electrons) == len(Muons) == 2:
            ZElMass = (Electrons[0].p4() + Electrons[1].p4()).M()
            ZMuMass = (Muons[0].p4() + Muons[1].p4()).M()
            if abs(ZElMass - ZMassPDG) < abs(ZMuMass - ZMassPDG):
                Z1Mass, Z1Leptons = ZElMass, Electrons
                Z2Mass, Z2Leptons = ZMuMass, Muons
            else:
                Z1Mass, Z1Leptons = ZMuMass, Muons
                Z2Mass, Z2Leptons = ZElMass, Electrons
            
        else:
            Z1Mass = float("inf")
            Z1Leptons = None
            for i, l1 in enumerate(Leptons):
                for j in range(i+1, len(Leptons)):
                    l2 = Leptons[j]
                    if l1.pdgId != -l2.pdgId: continue
                    LeptonsMass = (l1.p4() + l2.p4()).M()
                    if abs(LeptonsMass - ZMassPDG) < abs(Z1Mass - ZMassPDG):
                        Z1Mass = LeptonsMass
                        Z1Leptons = [l1, l2]
            Z2Leptons = [l for l in Leptons if l not in Z1Leptons]
            
        if (Z1Mass is None or Z1Leptons is None or Z2Leptons is None): return None
        if (len(Z1Leptons) != 2 or Z1Leptons[0].pdgId != - Z1Leptons[1].pdgId): return None
        if (len(Z2Leptons) != 2 or Z2Leptons[0].pdgId != - Z2Leptons[1].pdgId): return None
                        
        GenZ1Cand = ZClass(Z1Leptons)
        GenZ2Cand = ZClass(Z2Leptons)
        GenZZCand = ZZClass(GenZ1Cand, GenZ2Cand)

        return GenZZCand

    def ZZMassCuts(self, ZZ):

        if ZZ is None: return False
        
        if not 60 < ZZ.Z1mass < 120: return False
        if not 60 < ZZ.Z2mass < 120: return False
        
        return True

    
    # =========================
    # USEFUL FUNCTIONS: PHOTONS
    # =========================

    def GetTruePhotons(self, Photons):
        
        if len(Photons) < 1: return False, None
        TruePhotons = [ph for ph in Photons if (ph.cutBased >= 1
                                                and ph.electronVeto)]
        if len(TruePhotons) < 1: return False, None        

        return True, TruePhotons

    def PhotonsKinematicCuts(self, Photons):

        GoodPhotons = [ph for ph in Photons if ph.pt > 20
                                               and abs(ph.eta) < 2.4
                                               and not 1.444 < abs(ph.eta) < 1.566]
        if len(GoodPhotons) < 1: return False, None
        
        return True, GoodPhotons

    def GetIsolatedPhotons(self, Photons, Leptons):
        
        IsolatedPhotons = [ph for ph in Photons if all(ph.DeltaR(l) > 0.5 for l in Leptons)]
        if len(IsolatedPhotons) < 1: return False, None
        
        return True, IsolatedPhotons
        

    def GetGamma(self, Level, Photons, Leptons, ZZ):

        GammaCands = []
        for ph in Photons:
            ZGammaMinMass, ZMinMass = self.GetZGammaMinMass(Level, ph, Leptons, ZZ)
            if ZGammaMinMass > 100:
                GammaCands.append(ph)
        if len(GammaCands) < 1: return False, None

        Gamma = self.GetBestGamma(Level, GammaCands)
        if Gamma is None: return False, None
        
        if len(GammaCands) > 1:
            ZGammaMinMass, ZMinMass = self.GetZGammaMinMass(Level, Gamma, Leptons, ZZ)
        self.hEvent.fill1D(f"{Level}Level_ZGammaMinMass_PostCut", f"{Level}Level_ZGammaMinMass_PostCut", 52, 40., 300., ZGammaMinMass, self.weight)
        self.hEvent.fill2D(f"{Level}Level_ZGammaMinMass2D_PostCut", f"{Level}Level_ZMinMass(ZGammaMinMass)_PostCut", 104, 40., 300., 120, 60., 120., ZGammaMinMass, ZMinMass, self.weight)

        return True, Gamma

    def GetZGammaMinMass(self, Level, ph, Leptons, ZZ):

        Z1Mass = ZZ.Z1mass
        Z2Mass = ZZ.Z2mass
        
        if Level == "Gen":
            Z1GammaMass = (ZZ.Z1.p4 + ph.p4()).M()
            Z2GammaMass = (ZZ.Z2.p4 + ph.p4()).M()
        
        elif Level == "Reco":
            Z1GammaMass = (Leptons[0].p4() + Leptons[1].p4() + self.GetP4(ph)).M()
            Z1Mass = (Leptons[0].p4() + Leptons[1].p4()).M()
            Z2GammaMass = (Leptons[2].p4() + Leptons[3].p4() + self.GetP4(ph)).M()
            Z2Mass = (Leptons[2].p4() + Leptons[3].p4()).M()

        ZGammaMinMass = min(Z1GammaMass, Z2GammaMass)
        ZMinMass = Z1Mass if ZGammaMinMass == Z1GammaMass else Z2Mass
        self.hEvent.fill1D(f"{Level}Level_ZGammaMinMass_PreCut", f"{Level}Level_ZGammaMinMass_PreCut", 52, 40., 300., ZGammaMinMass, self.weight)
        self.hEvent.fill2D(f"{Level}Level_ZGammaMinMass2D_PreCut", f"{Level}Level_ZMinMass(ZGammaMinMass)_PreCut", 104, 40., 300., 120, 60., 120., ZGammaMinMass, ZMinMass, self.weight)
        
        return ZGammaMinMass, ZMinMass

    def GetBestGamma(self, Level, GammaCands):

        if len(GammaCands) == 1:
            BestGamma = GammaCands[0]

        elif len(GammaCands) > 1:
            BestGamma = max(GammaCands, key=lambda ph: ph.pt)
            Gamma2 = max((ph for ph in GammaCands if ph != BestGamma), key=lambda ph: ph.pt)
            self.hEvent.fill1D(f"{Level}Level_Gamma2Pt", f"{Level}Level_Gamma2Pt", 20, 20., 220., Gamma2.pt, self.weight)

        self.hEvent.fill1D(f"{Level}Level_BestGammaPt", f"{Level}Level_BestGammaPt", 50, 20., 220., BestGamma.pt, self.weight)
            
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

    def PlotsOI(self, Level, ZZ, Gamma):

        self.hEvent.fill1D(f"{Level}Level_Z1Mass", f"{Level}Level_Z1Mass", 100, 60., 120., ZZ.Z1mass, self.weight)
        self.hEvent.fill1D(f"{Level}Level_Z2Mass", f"{Level}Level_Z2Mass", 75, 60., 120., ZZ.Z2mass, self.weight)
        self.hEvent.fill1D(f"{Level}Level_ZZMass", f"{Level}Level_ZZMass", 40, 120., 520., ZZ.mass, self.weight)

        if Level == "Gen":
            ZZGammaInvMass = (ZZ.p4 + Gamma.p4()).M()
        elif Level == "Reco":
            ZZGammaInvMass = (ZZ.p4() + self.GetP4(Gamma)).M()
        self.hEvent.fill1D(f"{Level}Level_ZZGammaMass", f"{Level}Level_ZZGammaMass", 50, 120., 620., ZZGammaInvMass, self.weight)

        self.hEvent.fill1D(f"{Level}Level_GammaPt", f"{Level}Level_GammaPt", 18, 0., 180., Gamma.pt, self.weight)
        self.hEvent.fill1D(f"{Level}Level_GammaEta", f"{Level}Level_GammaEta", 12, 0., 2.4, Gamma.eta, self.weight)
               
    
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