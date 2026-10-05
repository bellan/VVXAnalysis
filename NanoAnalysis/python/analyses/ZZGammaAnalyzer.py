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
        bestCandIdx = self.event.bestCandIdx

        # =========================
        # GRAFICI CONTROLLO TAGLI
        # =========================

        LabelsTagliSignalDefinition = ["EventiTotali", "Trigger+BestCand", "CinematicaLeptoni","MasseZ1Z2", "CinematicaFotoniAndNotFsr"]
        self.hEvent.fill1D_label("EventiPerTaglio", "EventiPerTaglio", LabelsTagliSignalDefinition, LabelsTagliSignalDefinition[0], self.weight)
        
        LabelsTagliSignalRegion = ["No tagli", "Trigger+BestCand", "esiste theZZ", "ZZMass", "Cinematica leptoni"]
        self.hEvent.fill1D_label("EventiPerTaglioSignalRegion", "EventiPerTaglioSignalRegion", LabelsTagliSignalRegion, LabelsTagliSignalRegion[0], self.weight)

        if(bestCandIdx != -1 and self.event.HLT_passZZ4l): 
            
            # =========================
            # COLLECTIONS AND WEIGHTS
            # =========================

            GenParts = Collection(self.event, 'GenPart')
            GenLeptons = [p for p in GenParts if abs(p.pdgId) in (11, 13)
                                                 and p.status == 1
                                                 and (p.statusFlags & GenStatusFlag.SEL) == GenStatusFlag.SEL]
            GenPhotons = [p for p in GenParts if p.pdgId == 22
                                                 and p.status == 1
                                                 and (p.statusFlags & GenStatusFlag.SEL) == GenStatusFlag.SEL]
            # self.hEvent.fill1D("nGenLeptonsPostStatusFlags", "nGenLeptonsPostStatusFlags", 11, -0.5, 10.5, len(GenLeptons), self.weight)
            # self.hEvent.fill1D("nGenPhotonsPostStatusFlags", "nGenPhotonsPostStatusFlags", 6, -0.5, 5.5, len(GenPhotons), self.weight)
                       
            ZZs = Collection(self.event, 'ZZCand')
            theZZ = ZZs[bestCandIdx]
            Leptons = Collection(self.event, 'Lepton')
            Photons = Collection(self.event, 'Photon')
            FsrPhotons = Collection(self.event, 'FsrPhoton')
            # self.hEvent.fill1D("nLeptons", "nLeptons",11, -0.5, 10.5, len(Leptons), self.weight)
            # self.hEvent.fill1D("nPhotons", "nPhotons", 11, -0.5, 10.5, len(Photons), self.weight)
            # self.hEvent.fill1D("nFsrPhotons", "nFsrPhotons", 6, -0.5, 5.5, len(FsrPhotons), self.weight)

            self.hEvent.fill1D_label("EventiPerTaglio", "EventiPerTaglio", LabelsTagliSignalDefinition, LabelsTagliSignalDefinition[1], self.weight)
            self.hEvent.fill1D_label("EventiPerTaglioSignalRegion", "EventiPerTaglioSignalRegion", LabelsTagliSignalRegion, LabelsTagliSignalRegion[1], self.weight)

            # =========================
            # SIGNAL DEFINITION
            # =========================

            def SignalDefinition(GenLeptons, GenPhotons):

                if not self.analyzeMC: return False

                if not self.LeptonsKinematicCuts(GenLeptons): return False
                self.hEvent.fill1D_label("EventiPerTaglio", "EventiPerTaglio", LabelsTagliSignalDefinition, LabelsTagliSignalDefinition[2], self.weight)
                GenZZ = self.GetZZ(GenLeptons)
                if not self.ZZMassCuts(GenZZ): return False
                self.hEvent.fill1D_label("EventiPerTaglio", "EventiPerTaglio", LabelsTagliSignalDefinition, LabelsTagliSignalDefinition[3], self.weight)
        
                passed, GenGammaCands = self.PhotonsKinematicAndFsrCuts(GenPhotons, GenLeptons, GenZZ)
                if not passed: return False
                self.hEvent.fill1D_label("EventiPerTaglio", "EventiPerTaglio", LabelsTagliSignalDefinition, LabelsTagliSignalDefinition[4], self.weight)
                GenGamma = self.GetBestGamma(GenGammaCands)
                if GenGamma is None: return False
                PtThetaGraph = self.GraphGammaPtThetaZZG(GenGamma, GenZZ)
                
                return True

            ZZGammaSignalDefinition = SignalDefinition(GenLeptons, GenPhotons)

            
            # =========================
            # SIGNAL REGION
            # =========================

            # def MatchingControlZZ(theZZ)
            
            # def SignalRegion(Leptons, Photons):
                # if theZZ is None: return False
                # self.hEvent.fill1D_label("EventiPerTaglioSignalRegion", "EventiPerTaglioSignalRegion", LabelsTagliSignalRegion, LabelsTagliSignalRegion[2], self.weight)
                # if not 60 < theZZ.Z1mass < 120: return False
                # if not 60 < theZZ.Z2mass < 120: return False
                # self.hEvent.fill1D_label("EventiPerTaglioSignalRegion", "EventiPerTaglioSignalRegion", LabelsTagliSignalRegion, LabelsTagliSignalRegion[3], self.weight)
                
                # ll2Mass = (Leptons[theZZ.Z2l1Idx].p4() + Leptons[theZZ.Z2l2Idx].p4()).M()
                # self.hEvent.fill1D("ll2mass", "ll2mass", 100., -10, 120, ll2Mass, self.weight)
                
                # self.hEvent.fill1D("theZ1mass", "theZ1mass", 100., 60, 120, theZZ.Z1mass, self.weight)
                # self.hEvent.fill1D("theZ2mass", "theZ2mass", 75., 60, 120, theZZ.Z2mass, self.weight)               
                
                # return True

            # ZZGammaSignalRegion = SignalRegion(Leptons, Photons)
    

    # =========================
    # USEFUL FUNCTIONS: ZZ
    # =========================

    def LeptonsKinematicCuts(self, Leptons):
        if len(Leptons) != 4: return False
        if (sum(l.pdgId == 11 for l in Leptons) != sum(l.pdgId == -11 for l in Leptons)
            or sum(l.pdgId == 13 for l in Leptons) != sum(l.pdgId == -13 for l in Leptons)): return False
        if sum(l.pt > 5 for l in Leptons) < 4: return False
        if sum(l.pt > 10 for l in Leptons) < 2: return False
        if sum(l.pt > 20 for l in Leptons) < 1: return False
        if any(abs(l.eta) > 2.5 for l in Leptons): return False
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

    def ZZMassCuts(self, GenZZ):
        if GenZZ is None: return False
        if not 60 < GenZZ.Z1.mass < 120: return False
        if not 60 < GenZZ.Z2.mass < 120: return False
        self.hEvent.fill1D("GenZ1Mass", "GenZ1Mass", 100, 90., 92.4, GenZZ.Z1.mass, self.weight)
        self.hEvent.fill1D("GenZ2Mass", "GenZ2Mass", 75, 90., 92.4, GenZZ.Z2.mass, self.weight)
        return True

    
    # =========================
    # USEFUL FUNCTIONS: PHOTONS
    # =========================

    def PhotonsKinematicAndFsrCuts(self, Photons, Leptons, ZZ):
        if len(Photons) < 1: return False, None
        GoodPhotons = [ph for ph in Photons if ph.pt > 20
                                               and abs(ph.eta) < 2.4
                                               and not 1.444 < abs(ph.eta) < 1.566]
        if len(GoodPhotons) == 0: return False, None

        GammaCands = [ph for ph in GoodPhotons if all(ph.DeltaR(l) > 0.5 for l in Leptons)
                                                  and self.GetZGammaMassMin(ph, ZZ.Z1, ZZ.Z2) > 100]
        if len(GammaCands) < 1: return False, None
                
        return True, GammaCands
                
    def GetZGammaMassMin(self, ph, GenZ1, GenZ2):
        Z1GammaMass = (GenZ1.p4 + ph.p4()).M()
        Z1Mass = GenZ1.mass
        Z2GammaMass = (GenZ2.p4 + ph.p4()).M()
        Z2Mass = GenZ2.mass
        ZGammaMassMin = min(Z1GammaMass, Z2GammaMass)
        self.hEvent.fill1D("llGammaMassMin", "llGammaMassMin", 60, 60., 300., ZGammaMassMin, self.weight)
        if ZGammaMassMin == Z1GammaMass: ZMassMin = Z1Mass
        else: ZMassMin = Z2Mass
        self.hEvent.fill2D("llGammaMassMin2D", "llGammaMassMin2D", 50, 89., 94.4, 93, 60., 300., ZMassMin, ZGammaMassMin, self.weight)
        return ZGammaMassMin

    def GetBestGamma(self, GoodPhotons):
        BestGamma = max(GoodPhotons, key=lambda ph: ph.pt, default=None)
        if BestGamma is None: return None
        BestGammaPt = BestGamma.pt
        self.hEvent.fill1D("BestGammaPt", "BestGammaPt", 50, 20., 220., BestGammaPt, self.weight)
        Gamma2 = max((ph for ph in GoodPhotons if ph != BestGamma), key=lambda ph: ph.pt, default=None)
        if Gamma2 is not None:
            Gamma2Pt = Gamma2.pt
            self.hEvent.fill1D("Gamma2Pt", "Gamma2Pt", 50, 20., 120., Gamma2Pt, self.weight)
        return BestGamma

    def GraphGammaPtThetaZZG(self, Gamma, ZZ):
        Gamma_p4 = TLorentzVector()
        Gamma_p4.SetPtEtaPhiM(Gamma.pt, Gamma.eta, Gamma.phi, Gamma.mass)
        thetaZZGamma = ZZ.p4.Vect().Angle(Gamma_p4.Vect())
        self.hEvent.fill2D("BestGammaPt-AngoloZZGamma", "BestGammaPt-AngoloZZGamma", 50, 0., ROOT.TMath.Pi(), 5, 60., 100., thetaZZGamma , Gamma.pt, self.weight)
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