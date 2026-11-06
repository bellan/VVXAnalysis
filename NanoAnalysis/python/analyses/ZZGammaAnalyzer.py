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
        
        self.PhotonsID = "mva90"                     # choose between "cutL", "cutM", "cutT", "mva90" and "mva80"
        self.PhotonsOrder = "CutThenBestPt"          # choose between "CutThenBestPt" and "BestPtThenCut"
        self.AnalysisMode = "DividePromptPhotons"    # choose between "ZZGamma_Signal", "ZZGamma_NotSignal", "ZZGamma_EfficiencyPurity", "SR_only", "DividePromptPhotons"
        self.PromptPhotonsType = "other"             # choose between "noMother", "fromG", "fromQ", "fromLep", "fromPh", "other"

        self.SDCuts = ["TotalEvents", "LeptonsKinCuts", "Z1Z2Mass", "PhotonsKinCuts", "IsolatedPhotons", "TribosonRegion"]
        self.SRCuts = ["TotalEvents", "HLT_passZZ4l", "BestCand", "LeptonsKinCuts", "Z1Z2Mass", "ElectronVeto", "PhotonsKinCuts", f"PhotonsID_{self.PhotonsID}", "IsolatedPhotons", "TribosonRegion"]

    def analyze(self):
        
        self.LoadCollections()
        
        if self.AnalysisMode == "ZZGamma_Signal":
            isSD = self.SignalDefinition(self.GenLeptons, self.GenPhotons, self.PhotonsOrder)               # da modificare GenPhotons??
            if isSD:
                isSR = self.SignalRegion(self.Leptons, self.Photons, self.PhotonsOrder, self.PhotonsID)

        elif self.AnalysisMode == "ZZGamma_NotSignal":
            isSD = self.SignalDefinition(self.GenLeptons, self.GenPhotons, self.PhotonsOrder)               # da modificare GenPhotons??
            if not isSD:
                isSR = self.SignalRegion(self.Leptons, self.Photons, self.PhotonsOrder, self.PhotonsID)

        elif self.AnalysisMode == "ZZGamma_EfficiencyPurity":
            isSD = self.SignalDefinition(self.GenLeptons, self.GenPhotons, self.PhotonsOrder)               # da modificare GenPhotons??
            isSR = self.SignalRegion(self.Leptons, self.Photons, self.PhotonsOrder, self.PhotonsID)
            self.isBgivenA(isSD, isSR, "SD", "SR")
            self.isBgivenA(isSR, isSD, "SR", "SD")

        elif self.AnalysisMode == "SR_only":
            isSR = self.SignalRegion(self.Leptons, self.Photons, self.PhotonsOrder, self.PhotonsID)

        elif self.AnalysisMode == "DividePromptPhotons":

            phTypes = list(self.GenPromptPhotons.keys())
            idx = phTypes.index(self.PromptPhotonsType)
            GenPhotonsOI = self.GenPromptPhotons[self.PromptPhotonsType]

            if len(GenPhotonsOI) > 0 and not any(self.GenPromptPhotons[key] for key in phTypes[:idx]):
                
                isSD = self.SignalDefinition(self.GenLeptons, GenPhotonsOI, self.PhotonsOrder)
                isSR = self.SignalRegion(self.Leptons, self.Photons, self.PhotonsOrder, self.PhotonsID)

        
    # =========================
    # SIGNAL DEF and REGION
    # =========================

    def SignalDefinition(self, GenLeptons, GenPhotons, PhotonsOrder):

        if not self.analyzeMC: return False
        self.CutPlot("Definition", 0, PhotonsOrder)

        if not self.LeptonsKinematicCuts("Gen", GenLeptons): return False
        self.CutPlot("Definition", 1, PhotonsOrder)
        GenZZ = self.GetZZ(GenLeptons)
        if GenZZ is None: return False
        if not self.ZZMassCuts(GenZZ): return False
        self.CutPlot("Definition", 2, PhotonsOrder)
        
        KinPass, GoodPhotons = self.PhotonsKinematicCuts(GenPhotons)
        if not KinPass: return False
        self.CutPlot("Definition", 3, PhotonsOrder)
        IsolatedPass, IsolatedPhotons = self.GetIsolatedPhotons(GoodPhotons, GenLeptons)
        if not IsolatedPass: return False
        self.CutPlot("Definition", 4, PhotonsOrder)
        TribosonPass, GenGamma = self.TribosonRegion("Gen", IsolatedPhotons, GenLeptons, GenZZ, PhotonsOrder)
        if not TribosonPass: return False
        self.CutPlot("Definition", 5, PhotonsOrder)

        self.PlotsOI("Gen", GenZZ, GenGamma)

        return True

    def SignalRegion(self, Leptons, Photons, PhotonsOrder, PhotonsID):

        self.CutPlot("Region", 0, PhotonsOrder)

        if not self.event.HLT_passZZ4l: return False
        self.CutPlot("Region", 1, PhotonsOrder)
        bestCandIdx = self.event.bestCandIdx
        if not (bestCandIdx != -1): return False
        self.CutPlot("Region", 2, PhotonsOrder)
        RecoZZ = self.ZZs[bestCandIdx]
        LeptonsOI = [Leptons[RecoZZ.Z1l1Idx], Leptons[RecoZZ.Z1l2Idx], Leptons[RecoZZ.Z2l1Idx], Leptons[RecoZZ.Z2l2Idx]]
        if not self.LeptonsKinematicCuts("Reco", LeptonsOI): return False
        self.CutPlot("Region", 3, PhotonsOrder)
        if not self.ZZMassCuts(RecoZZ): return False
        self.CutPlot("Region", 4, PhotonsOrder)

        self.PlotsOI_ControlRegionL4P("Reco", RecoZZ)

        if not self.analyzeMC: return False

        ElectronVetoPass, TruePhotons = self.GetTruePhotons(Photons)
        if not ElectronVetoPass: return False
        self.CutPlot("Region", 5, PhotonsOrder)
        KinPass, GoodPhotons = self.PhotonsKinematicCuts(TruePhotons)
        if not KinPass: return False
        self.CutPlot("Region", 6, PhotonsOrder)
        
        self.hEvent.fill1D("nPhotonsPreTribosonregion", "nPhotonsPreTribosonregion", 11, -0.5, 10.5, len(Photons), self.weight)
        
        IDPass, GoodIDPhotons = self.GetGoodIDPhotons(GoodPhotons, PhotonsID)
        if not IDPass: return False
        self.CutPlot("Region", 7, PhotonsOrder)
        IsolatedPass, IsolatedPhotons = self.GetIsolatedPhotons(GoodIDPhotons, LeptonsOI)
        if not IsolatedPass: return False
        self.CutPlot("Region", 8, PhotonsOrder)
        TribosonPass, RecoGamma = self.TribosonRegion("Reco", IsolatedPhotons, LeptonsOI, RecoZZ, PhotonsOrder)
        if not TribosonPass: return False
        self.CutPlot("Region", 9, PhotonsOrder)

        self.PlotsOI_PhotonID(RecoGamma, PhotonsID)
        self.PlotsOI("Reco", RecoZZ, RecoGamma)

        return True
    
    
    # =========================
    # LEPTONS and ZZ
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
    # PHOTONS
    # =========================

    def GetTruePhotons(self, Photons):
        
        if len(Photons) < 1: return False, None
        TruePhotons = [ph for ph in Photons if ph.electronVeto]
        if len(TruePhotons) < 1: return False, None        

        return True, TruePhotons

    def PhotonsKinematicCuts(self, Photons):
        
        if len(Photons) < 1: return False, None
        GoodPhotons = [ph for ph in Photons if ph.pt > 20
                                               and abs(ph.eta) < 2.4
                                               and not 1.444 < abs(ph.eta) < 1.566]
        if len(GoodPhotons) < 1: return False, None
        
        return True, GoodPhotons

    def GetIsolatedPhotons(self, Photons, Leptons):
        
        if len(Photons) < 1: return False, None
        IsolatedPhotons = [ph for ph in Photons if all(ph.DeltaR(l) > 0.5 for l in Leptons)]
        if len(IsolatedPhotons) < 1: return False, None
        
        return True, IsolatedPhotons

    def GetGoodIDPhotons(self, Photons, PhotonsID):

        if len(Photons) < 1: return False, None

        if PhotonsID == "mva90":
            GoodIDPhotons = [ph for ph in Photons if ph.mvaID_WP90]
        elif PhotonsID == "mva80":
            GoodIDPhotons = [ph for ph in Photons if ph.mvaID_WP80]
        elif PhotonsID == "cutL":
            GoodIDPhotons = [ph for ph in Photons if ph.cutBased >= 1]
        elif PhotonsID == "cutM":
            GoodIDPhotons = [ph for ph in Photons if ph.cutBased >= 2]
        elif PhotonsID == "cutT":
            GoodIDPhotons = [ph for ph in Photons if ph.cutBased >= 3]
        else:
            raise ValueError(f"Unknown PhotonsID: {PhotonsID}")
            
        if len(GoodIDPhotons) < 1: return False, None

        return True, GoodIDPhotons
    
    
    # =========================
    # TRIBOSON REGION
    # =========================

    def TribosonRegion(self, Level, Photons, Leptons, ZZ, PhotonsOrder):

        if PhotonsOrder == "CutThenBestPt":
            TribosonPass, Gamma = self.GetGamma_CutThenBestPt(Level, Photons, Leptons, ZZ, PhotonsOrder)
            
        elif PhotonsOrder == "BestPtThenCut":
            TribosonPass, Gamma = self.GetGamma_BestPtThenCut(Level, Photons, Leptons, ZZ, PhotonsOrder)
            
        return TribosonPass, Gamma
        
    def GetGamma_CutThenBestPt(self, Level, Photons, Leptons, ZZ, PhotonsOrder):

        if len(Photons) < 1: return False, None

        ZGammaMinMass_MAX = 0
        ZMinMass_MAX = 0
        GammaCands = []

        for ph in Photons:
            ZGammaMinMass, ZMinMass = self.GetZGammaMinMass(Level, ph, Leptons, ZZ)
            
            if ZGammaMinMass > ZGammaMinMass_MAX:
                ZGammaMinMass_MAX = ZGammaMinMass
                ZMinMass_MAX = ZMinMass

            if ZGammaMinMass > 100:
                GammaCands.append(ph)
            
        self.hEvent.fill1D(f"{Level}Level_ZGammaMinMass_PreCut", f"{Level}Level_ZGammaMinMass_PreCut", 60, 0., 300., ZGammaMinMass_MAX, self.weight)
        self.hEvent.fill2D(f"{Level}Level_ZGammaMinMass2D_PreCut", f"{Level}Level_ZMinMass(ZGammaMinMass)_PreCut", 60, 0., 300., 24, 0., 120., ZGammaMinMass_MAX, ZMinMass_MAX, self.weight)
           
        if len(GammaCands) < 1: return False, None

        Gamma = self.GetBestGamma(Level, GammaCands, PhotonsOrder)
        ZGammaMinMass, ZMinMass = self.GetZGammaMinMass(Level, Gamma, Leptons, ZZ)
        self.hEvent.fill1D(f"{Level}Level_ZGammaMinMass_PostCut", f"{Level}Level_ZGammaMinMass_PostCut", 60, 0., 300., ZGammaMinMass, self.weight)
        self.hEvent.fill2D(f"{Level}Level_ZGammaMinMass2D_PostCut", f"{Level}Level_ZMinMass(ZGammaMinMass)_PostCut", 60, 0., 300., 24, 0., 120., ZGammaMinMass, ZMinMass, self.weight)

        return True, Gamma

    def GetGamma_BestPtThenCut(self, Level, Photons, Leptons, ZZ, PhotonsOrder):
        
        if len(Photons) < 1: return False, None

        BestGamma = self.GetBestGamma(Level, Photons, PhotonsOrder)
        ZGammaMinMass, ZMinMass = self.GetZGammaMinMass(Level, BestGamma, Leptons, ZZ)
        self.hEvent.fill1D(f"{Level}Level_ZGammaMinMass_PreCut_BestPtThenCut", f"{Level}Level_ZGammaMinMass_PreCut_BestPtThenCut", 60, 0., 300., ZGammaMinMass, self.weight)
        self.hEvent.fill2D(f"{Level}Level_ZGammaMinMass2D_PreCut_BestPtThenCut", f"{Level}Level_ZMinMass(ZGammaMinMass)_PreCut_BestPtThenCut", 60, 0., 300., 24, 0., 120., ZGammaMinMass, ZMinMass, self.weight)

        if ZGammaMinMass > 100:
            self.hEvent.fill1D(f"{Level}Level_ZGammaMinMass_PostCut_BestPtThenCut", f"{Level}Level_ZGammaMinMass_PostCut_BestPtThenCut", 60, 0., 300., ZGammaMinMass, self.weight)
            self.hEvent.fill2D(f"{Level}Level_ZGammaMinMass2D_PostCut_BestPtThenCut", f"{Level}Level_ZMinMass(ZGammaMinMass)_PostCut_BestPtThenCut", 60, 0., 300., 24, 0., 120., ZGammaMinMass, ZMinMass, self.weight)
            return True, BestGamma

        return False, None

    def GetZGammaMinMass(self, Level, ph, Leptons, ZZ):

        Z1Mass = ZZ.Z1mass
        Z2Mass = ZZ.Z2mass
        
        if Level == "Gen":
            Z1GammaMass = (ZZ.Z1.p4 + ph.p4()).M()
            Z2GammaMass = (ZZ.Z2.p4 + ph.p4()).M()
        
        elif Level == "Reco":
            Z1GammaMass = (Leptons[0].p4() + Leptons[1].p4() + self.GetP4(ph)).M()
            Z2GammaMass = (Leptons[2].p4() + Leptons[3].p4() + self.GetP4(ph)).M()
            
        ZGammaMinMass = min(Z1GammaMass, Z2GammaMass)
        ZMinMass = Z1Mass if ZGammaMinMass == Z1GammaMass else Z2Mass
        
        return ZGammaMinMass, ZMinMass

    def GetBestGamma(self, Level, GammaCands, PhotonsOrder):

        if len(GammaCands) < 1: return None

        elif len(GammaCands) == 1:
            BestGamma = GammaCands[0]

        elif len(GammaCands) > 1:
            BestGamma = max(GammaCands, key=lambda ph: ph.pt)
            Gamma2 = max((ph for ph in GammaCands if ph != BestGamma), key=lambda ph: ph.pt)
            if PhotonsOrder == "CutThenBestPt":
                self.hEvent.fill1D(f"{Level}Level_Gamma2Pt", f"{Level}Level_Gamma2Pt", 20, 20., 220., Gamma2.pt, self.weight)
            elif PhotonsOrder == "BestPtThenCut":
                self.hEvent.fill1D(f"{Level}Level_Gamma2Pt_{PhotonsOrder}", f"{Level}Level_Gamma2Pt_{PhotonsOrder}", 20, 20., 220., Gamma2.pt, self.weight)

        if PhotonsOrder == "CutThenBestPt":
            self.hEvent.fill1D(f"{Level}Level_BestGammaPt", f"{Level}Level_BestGammaPt", 50, 20., 220., BestGamma.pt, self.weight)
        elif PhotonsOrder == "BestPtThenCut":
            self.hEvent.fill1D(f"{Level}Level_BestGammaPt_{PhotonsOrder}", f"{Level}Level_BestGammaPt_{PhotonsOrder}", 50, 20., 220., BestGamma.pt, self.weight)
            
        return BestGamma

    
    # =========================
    # AUXILIARY FUNCTIONS
    # =========================

    def LoadCollections(self):

        if self.analyzeMC:

            self.GenParts = Collection(self.event, 'GenPart')
            self.GenLeptons = [p for p in self.GenParts if abs(p.pdgId) in (11, 13)
                                                           and p.status == 1
                                                           and (p.statusFlags & GenStatusFlag.SEL) == GenStatusFlag.SEL]

            self.AllGenPhotons = [p for p in self.GenParts if p.pdgId == 22]
                                                            
            self.GenPromptPhotons = {"noMother": [], "fromQ": [], "fromG": [], "fromLep": [], "fromPh": [], "other": []}

            for ph in self.AllGenPhotons:

                if (ph.statusFlags & GenStatusFlag.PH_SEL) != GenStatusFlag.PH_SEL:
                    continue

                motherIdx = ph.genPartIdxMother

                if motherIdx == -1:
                    phType = "noMother"

                else:
                    motherID = self.GenParts[motherIdx].pdgId

                    if abs(motherID) in (1, 2, 3, 4, 5, 6):
                        phType = "fromQ"
                    elif abs(motherID) == 21:
                        phType = "fromG"
                    elif abs(motherID) in (11, 13, 15):
                        phType = "fromLep"
                    elif motherID == 22:
                        phType = "fromPh"
                    else:
                        phType = "other"
                
                self.GenPromptPhotons[phType].append(ph)

            self.GenPhotons = [p for p in self.GenParts if p.pdgId == 22
                                                           and p.status == 1
                                                           and (p.statusFlags & GenStatusFlag.SEL) == GenStatusFlag.SEL]            
      
        self.ZZs = Collection(self.event, 'ZZCand')
        self.Leptons = Collection(self.event, 'Lepton')
        self.Photons = Collection(self.event, 'Photon')
            
    def GetP4(self, p):

        p4 = TLorentzVector()
        if hasattr(p, "mass"):
            p4.SetPtEtaPhiM(p.pt, p.eta, p.phi, p.mass)
        else:
            p4.SetPtEtaPhiM(p.pt, p.eta, p.phi, 0.)
        
        return p4

    def isBgivenA(self, A, B, nomeA, nomeB):

        if A:
            if B:
                result = 1
            else: result = 0
            self.hEvent.fill1D(f"is{nomeB}given{nomeA}", f"is{nomeB}given{nomeA}", 2, -0.5, 1.5, result, self.weight)

    def PhotonWPCutPlot(self, Gamma, PhotonID):

        if PhotonID == "mva90":
            if not Gamma.mvaID_WP90: return
            self.hEvent.fill1D_label("PhotonWP", "PhotonWP_CutFlow", self.PhotonWPLabels, self.PhotonWPLabels[0], self.weight)

        elif PhotonID == "mva80":
            if not Gamma.mvaID_WP80: return
            self.hEvent.fill1D_label("PhotonWP", "PhotonWP_CutFlow", self.PhotonWPLabels, self.PhotonWPLabels[1], self.weight)
        
        elif PhotonID == "cutL":
            if not Gamma.cutBased >= 1: return
            self.hEvent.fill1D_label("PhotonWP", "PhotonWP_CutFlow", self.PhotonWPLabels, self.PhotonWPLabels[3], self.weight)
        
        elif PhotonID == "cutM":
            if not Gamma.cutBased >= 2: return
            self.hEvent.fill1D_label("PhotonWP", "PhotonWP_CutFlow", self.PhotonWPLabels, self.PhotonWPLabels[4], self.weight)
        
        elif PhotonID == "cutT":
            if not Gamma.cutBased >= 3: return
            self.hEvent.fill1D_label("PhotonWP", "PhotonWP_CutFlow", self.PhotonWPLabels, self.PhotonWPLabels[5], self.weight)

    PhotonWPLabels = ["pass_mva90", "pass_mva80", "" , "pass_CutL", "pass_CutM", "pass_CutT"]

    def GraphGammaPtThetaZZG(self, Level, Gamma, ZZ):
        
        if Level == "Gen":
            thetaZZGamma = ZZ.p4.Vect().Angle(Gamma.p4().Vect())
        elif Level == "Reco":
            thetaZZGamma = ZZ.p4().Vect().Angle(self.GetP4(Gamma).Vect())
        
        self.hEvent.fill2D(f"{Level}Level_BestGammaPt-AngoloZZGamma", f"{Level}Level_BestGammaPt-AngoloZZGamma", 50, 0., ROOT.TMath.Pi(), 5, 60., 100., thetaZZGamma , Gamma.pt, self.weight)

    
    # =========================
    # PLOTS OF INTEREST
    # =========================

    def CutPlot(self, Level, Cut, Order = None):

        if Level == "Definition":
            Labels = self.SDCuts
        elif Level == "Region":
            Labels = self.SRCuts
        if Order is None:
            self.hEvent.fill1D_label(f"Signal{Level}_CutFlow", f"Signal{Level}_CutFlow", Labels, Labels[Cut], self.weight)
        else:
            self.hEvent.fill1D_label(f"Signal{Level}_CutFlow_{Order}", f"Signal{Level}_CutFlow_{Order}", Labels, Labels[Cut], self.weight)
        
    def PlotsOI_ControlRegionL4P(self, Level, ZZ):
        
        self.hEvent.fill1D(f"{Level}Level_Z1Mass_L4P", f"{Level}Level_Z1Mass_L4P", 60, 60., 120., ZZ.Z1mass, self.weight)
        self.hEvent.fill1D(f"{Level}Level_Z2Mass_L4P", f"{Level}Level_Z2Mass_L4P", 60, 60., 120., ZZ.Z2mass, self.weight)
        self.hEvent.fill1D(f"{Level}Level_ZZMass_L4P", f"{Level}Level_ZZMass_L4P", 10, 0., 1000., ZZ.mass, self.weight)

    def PlotsOI(self, Level, ZZ, Gamma):

        # self.GraphGammaPtThetaZZG(Level, Gamma, ZZ)
        
        self.hEvent.fill1D(f"{Level}Level_Z1Mass", f"{Level}Level_Z1Mass", 60, 60., 120., ZZ.Z1mass, self.weight)
        self.hEvent.fill1D(f"{Level}Level_Z2Mass", f"{Level}Level_Z2Mass", 60, 60., 120., ZZ.Z2mass, self.weight)
        self.hEvent.fill1D(f"{Level}Level_ZZMass", f"{Level}Level_ZZMass", 10, 0., 1000., ZZ.mass, self.weight)

        if Level == "Gen":
            ZZGammaInvMass = (ZZ.p4 + Gamma.p4()).M()
        elif Level == "Reco":
            ZZGammaInvMass = (ZZ.p4() + self.GetP4(Gamma)).M()
        self.hEvent.fill1D(f"{Level}Level_ZZGammaMass", f"{Level}Level_ZZGammaMass", 10, 0., 1000., ZZGammaInvMass, self.weight)

        self.hEvent.fill1D(f"{Level}Level_GammaPt", f"{Level}Level_GammaPt", 12, 0., 180., Gamma.pt, self.weight)
        self.hEvent.fill1D(f"{Level}Level_abs(GammaEta)", f"{Level}Level_abs(GammaEta)", 12, 0., 2.4, abs(Gamma.eta), self.weight)

    def PlotsOI_PhotonID(self, Gamma, PhotonsID):

        self.PhotonWPCutPlot(Gamma, PhotonsID)

        if PhotonsID == "mva90" or PhotonsID == "mva80":
            IDvar = Gamma.mvaID
            xmin = -1.1
            xmax = 1.1
            nbin = 22

        elif PhotonsID == "cutL" or PhotonsID == "cutM" or PhotonsID == "cutT":
            IDvar = Gamma.cutBased
            xmin = 0
            xmax = 4
            nbin = 4

        self.hEvent.fill1D(f"Photon_{PhotonsID}_all", f"Photon_{PhotonsID}_all", nbin, xmin, xmax, IDvar, self.weight)

        if Gamma.isScEtaEB:
            self.hEvent.fill1D(f"Photon_{PhotonsID}_barrel", f"Photon_{PhotonsID}_barrel", nbin, xmin, xmax, IDvar, self.weight)
        elif Gamma.isScEtaEE:
            self.hEvent.fill1D(f"Photon_{PhotonsID}_endcap", f"Photon_{PhotonsID}_endcap", nbin, xmin, xmax, IDvar, self.weight)
        
    
# =========================
# USEFUL CLASSES
# =========================

class GenStatusFlag(IntFlag):
    IS_PROMPT = 1 << 0
    FROM_HARD_PROCESS = 1 << 8

    SEL = IS_PROMPT | FROM_HARD_PROCESS
    PH_SEL = IS_PROMPT

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