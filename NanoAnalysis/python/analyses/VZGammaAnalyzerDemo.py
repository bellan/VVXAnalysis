from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection
from ZZAnalysis.NanoAnalysis.tools import getLeptons, get_genEventSumw

from VVXAnalysis.NanoAnalysis.EventAnalyzer import EventAnalyzer
from VVXAnalysis.NanoAnalysis.Histogrammer import *
from VVXAnalysis.NanoAnalysis.Regions import Flags as Regions

import ROOT
import math
import itertools
from ROOT import TDatabasePDG 
from enum import IntFlag

class VZGammaAnalyzerDemo(EventAnalyzer, analysis_name="VZGammaAnalyzerDemo"):

    
    def __init__(self, regions):
        super().__init__(regions)


    def loadCollections(self):

        if self.analyzeMC:
            self.GenParts = Collection(self.event, 'GenPart')
            self.GenJets = Collection(self.event, 'GenJet')
            self.GenLeptons = []
            self.GenPhotons = []
            self.GenQuarks = []

            for p in self.GenParts:
                if (abs(p.pdgId) in (11,13) and p.status == 1 and (p.statusFlags & GenStatusFlag.SEL) == GenStatusFlag.SEL):
                    self.GenLeptons.append(p)
                elif (p.pdgId == 22 and p.status == 1 and (p.statusFlags & GenStatusFlag.SEL) == GenStatusFlag.SEL):
                    self.GenPhotons.append(p)
                elif (1 <= abs(p.pdgId) <= 6 and p.status == 23 and (p.statusFlags & GenStatusFlag.QUARK_SEL) == GenStatusFlag.QUARK_SEL):
                    self.GenQuarks.append(p)

        self.event.SetBranchStatus("*Lepton*", 1)
        self.event.SetBranchStatus("*ZCand*", 1)
        self.event.SetBranchStatus("bestZIdx", 1)
        self.event.SetBranchStatus("*FsrPhoton*", 1)

        #self.Leptons = Collection(self.event, 'Lepton')
        self.Leptons = self.leptons
        self.Photons = Collection(self.event, 'Photon')
        #self.Jets = Collection(self.event, 'Jet')
        NoIdJets = Collection(self.event, 'Jet')
        self.Jets = [j for j in NoIdJets if j.jetId >= 2]
        self.ZCands = Collection(self.event, 'ZCand')
        self.bestZIdx = self.event.bestZIdx
        self.FSRPhotons = Collection(self.event, 'FsrPhoton')
            
        
    def analyze(self):
        self.CutFlowLabelsSD=["TotalEvents","QuarkNum","LeptonNumKin","LeptonZCand","PhotonKin","JetKin","FromLepIso",
                              "PhotonJetIso","JetVCand"]
        self.CutFlowLabelsSR=["TotalEvents","Trigger","ZCandMass","LeptonNumKin","ElectronVeto","PhotonKin","FromLepIso",
                              "PhotonJetIso","PhotonId","JetKin","JetVCand"]#,"ZCandMass100"]

        self.fillCutFlow("SD","TotalEvents")
        self.fillCutFlow("SR","TotalEvents")

        self.loadCollections()
        
        is_SD = self.IsSignal()
        is_SR = self.IsInSignalRegion(studyPlot=True)

        if self.analyzeMC:
            self.hEvent.fill2D("Confusion_Matrix", "Confusion_Matrix;Passes Signal Definition (GEN);Passes Signal Region (RECO)", 
                               2, -0.5, 1.5, 2, -0.5, 1.5, int(is_SD), int(is_SR), self.weight) 

        
        
        











    #=================
    #SIGNAL DEFINITION
    #=================  

    def IsSignal(self):

        if not self.analyzeMC: return False

        # Number of Jets before lepton kinematic Cut
        self.hEvent.fill1D("nJet_BeforeLepCut_GEN", "nJet_BeforeLepCut_GEN; N_{j};", 15, -0.5, 14.5, len(self.GenJets), self.weight)
        # Number of Quarks before lepton kinematic Cut
        self.hEvent.fill1D("nQuarks_BeforeLepCut_GEN", "nQuarks_BeforeLepCut_GEN; N_{q};", 15, -0.5, 14.5, len(self.GenQuarks), self.weight)
        # Number of Leptons before lepton kinematic Cut
        self.hEvent.fill1D("nLepton_BeforeLepCut_GEN", "nLepton_BeforeLepCut_GEN; N_{l};", 15, -0.5, 14.5, len(self.GenLeptons), self.weight)

        # -- QUARKS -- #
        if len(self.GenQuarks) != 2: return False
        self.fillCutFlow("SD","QuarkNum")       

        # -- LEPTONS -- #
        self.GenKinCutLeptons = self.leptonKinCut(self.GenLeptons, "Definition")
        if len(self.GenKinCutLeptons) != 2: return False
        self.fillCutFlow("SD","LeptonNumKin")

        self.GenZCandLeptons = self.getZLeptons(self.GenKinCutLeptons)
        if len(self.GenZCandLeptons) < 2:  return False
        self.fillCutFlow("SD","LeptonZCand")

        # Invariant Mass of the two best Leptons after cuts Distribution
        m_ll = self.getInvMass(self.GenZCandLeptons[0], self.GenZCandLeptons[1])
        self.hEvent.fill1D("InvMass_LepZCand_PreCut_GEN", "InvMass_LepZCand_PreCut_GEN;m_{ll} [GeV];Events", 60, 0., 180., m_ll, self.weight)
        # Number of Leptons after Lepton Cuts
        self.hEvent.fill1D("nLepton_AfterLepCut_GEN", "nLepton_AfterLepCut_GEN; N_{l};", 15, -0.5, 14.5, len(self.GenZCandLeptons), self.weight)
        # Number of Jets after Lepton Cuts
        self.hEvent.fill1D("nJet_AfterLepCut_GEN", "nJet_AfterLepCut_GEN; N_{j};", 15, -0.5, 14.5, len(self.GenJets), self.weight)
        # Number of Quarks after Lepton Cuts
        self.hEvent.fill1D("nQuarks_AfterLepCut_GEN", "nQuarks_AfterLepCut_GEN; N_{q};", 15, -0.5, 14.5, len(self.GenQuarks), self.weight)

        # -- JETS & PHOTONS -- #
        self.GenKinCutPhotons = self.photonKinCut(self.GenPhotons, "Definition")
        if len(self.GenKinCutPhotons) == 0: return False
        self.fillCutFlow("SD","PhotonKin")

        self.GenKinCutJets=self.jetKinCut(self.GenJets)
        if len(self.GenKinCutJets) < 2: return False
        self.fillCutFlow("SD","JetKin")

        # Number of Jets after kinematic Cut
        self.hEvent.fill1D("nJet_Kin_GEN", "nJet_Kin_GEN; N_{j};", 15, -0.5, 14.5, len(self.GenKinCutJets), self.weight)
        # Number of Quarks after kinematic Cut
        self.hEvent.fill1D("nQuarks_Kin_GEN", "nQuarks_Kin_GEN; N_{q};", 15, -0.5, 14.5, len(self.GenQuarks), self.weight)

        # DeltaR between photons and 3-most-energetic-jets Distribution 
        DeltaRGammaJet = self.PhotonJetDeltaRCut(self.GenKinCutPhotons,self.GenKinCutJets)
        for dR in DeltaRGammaJet:
            self.hEvent.fill1D("DeltaR_GammaJet_GEN", "DeltaR_GammaJet_GEN; #DeltaR(#gamma, j);", 60, 0., 6., dR, self.weight)
        for i, j in enumerate(self.GenKinCutJets[:3]):
            for p in self.GenKinCutPhotons:
                dR = p.DeltaR(j)
                hist_name = f"DeltaR_GammaJet{i+1}_GEN"
                hist_title = f"DeltaR Gamma - Jet {i+1} GEN;#DeltaR(#gamma, j_{i+1});# pairs"
                self.hEvent.fill1D(hist_name, hist_title, 60, 0., 6., dR, self.weight)
        # Pt Distributions of the photons and jet that have DeltaR between them <0.5
        for p in self.GenKinCutPhotons:
            for j in self.GenKinCutJets:
                if p.DeltaR(j) < 0.4:    
                    self.hEvent.fill1D("Pt_Photon_Matched_GEN", "Pt_Photon_Matched_GEN;p_{T}^{#gamma} [GeV];Events", 50, 0., 200., p.pt, self.weight)
                    self.hEvent.fill1D("Pt_Jet_Matched_GEN", "Pt_Jet_Matched_GEN;p_{T}^{j} [GeV];Events", 50, 0., 200., j.pt, self.weight)
                    delta_pt_norm = (j.pt - p.pt)/j.pt
                    self.hEvent.fill1D("DeltaPt_JetPhoton_GEN", "DeltaPt_JetPhoton_GEN;p_{T}^{j} - p_{T}^{#gamma} [GeV];Events", 60, -15., 15., delta_pt_norm, self.weight)

        self.GenNoLepPhotons = [p for p in self.GenKinCutPhotons if all(p.DeltaR(l) > 0.5 for l in self.GenZCandLeptons)]
        self.GenNoLepJets = [j for j in self.GenKinCutJets if all(j.DeltaR(l) > 0.5 for l in self.GenZCandLeptons)]
        if len(self.GenNoLepPhotons) == 0 or len(self.GenNoLepJets) < 2: return False
        self.fillCutFlow("SD","FromLepIso")

        self.GenGoodPhotons, _, self.GenGoodJets, _ = self.disambiguatePhotonJet(
            self.GenNoLepPhotons, self.GenNoLepJets, deltaRThr=0.1)

        if len(self.GenGoodPhotons) == 0 or len(self.GenGoodJets) < 2: return False
        self.fillCutFlow("SD","PhotonJetIso")

        #self.GenIsoPhotons = self.photonIsolationCut(self.GenKinCutPhotons, self.GenZCandLeptons, self.GenKinCutJets)
        #if len(self.GenIsoPhotons) == 0: return False
        #self.fillCutFlow("SD","PhotonIso")

        #self.GenIsoJets = self.jetIsolationCut(self.GenKinCutJets,self.GenZCandLeptons,self.GenKinCutPhotons)
        #if len(self.GenIsoJets) < 2: return False
        #self.fillCutFlow("SD","JetIso")

        # Number of Jets after DeltaR Cut
        self.hEvent.fill1D("nJet_AfterDeltaRCut_GEN", "nJet_AfterDeltaRCut_GEN; N_{j};", 15, -0.5, 14.5, len(self.GenGoodJets), self.weight)
        # Number of Quarks after DeltaR Cut
        self.hEvent.fill1D("nQuarks_AfterDeltaRCut_GEN", "nQuarks_AfterDeltaRCut_GEN; N_{q};", 15, -0.5, 14.5, len(self.GenQuarks), self.weight)

        self.GenVCandJets = self.getVJets(self.GenGoodJets)
        if len(self.GenVCandJets) < 2: return False
        self.fillCutFlow("SD","JetVCand")

        

        # Invariant Mass of the two best Leptons after cuts Distribution
        m_ll = self.getInvMass(self.GenZCandLeptons[0], self.GenZCandLeptons[1])
        self.hEvent.fill1D("InvMass_LepZCand_GEN", "InvMass_LepZCand_GEN;m_{ll} [GeV];Events", 60, 0., 180., m_ll, self.weight)        
        # Invariant Mass of 2l Gamma Distribution
        leading_photon = max(self.GenGoodPhotons, key=lambda ph: ph.pt)
        m_llgamma = self.getInvMass(self.GenZCandLeptons[0], self.GenZCandLeptons[1], leading_photon)
        self.hEvent.fill1D("InvMass_llgamma_GEN", "InvMass_llgamma_GEN;m_{ll#gamma} [GeV];Events", 60, 60., 240., m_llgamma, self.weight)
        # Invariant Mass of 2l Gamma Distribution vs Invariant Mass of 2l Distribution
        self.hEvent.fill2D("InvMass_ll_llgamma_GEN", "InvMass_ll_llgamma_GEN;m_{ll} [GeV];m_{ll#gamma} [GeV];Events", 36, 0., 180., 50, 0., 250., m_ll, m_llgamma, self.weight)
        # Invariant Mass of the two best Jet candidate to form V-Boson Distribution
        m_jj = self.getInvMass(self.GenVCandJets[0], self.GenVCandJets[1])
        self.hEvent.fill1D("InvMass_VJets_GEN", "InvMass_VJets_GEN;m_{jj} [GeV];Events", 50, 40., 140., m_jj, self.weight)
        # Invariant Mass of 2j Gamma Distribution
        m_jjgamma = self.getInvMass(self.GenVCandJets[0], self.GenVCandJets[1], leading_photon)
        self.hEvent.fill1D("InvMass_jjgamma_GEN", "InvMass_jjgamma_GEN;m_{jj#gamma} [GeV];Events", 60, 60., 240., m_jjgamma, self.weight)
        # Invariant Mass of 2j Gamma Distribution vs Invariant Mass of 2j Distribution
        self.hEvent.fill2D("InvMass_jj_jjgamma_GEN", "InvMass_jj_jjgamma_GEN;m_{jj} [GeV];Events;m_{jj#gamma} [GeV]", 36, 0., 180., 50, 0., 250., m_jj, m_jjgamma, self.weight)
        # Invariant Mass of 2j 2l Gamma Distribution
        m_lljjgamma = self.getInvMass(self.GenZCandLeptons[0], self.GenZCandLeptons[1], leading_photon, self.GenVCandJets[0], self.GenVCandJets[1])
        self.hEvent.fill1D("InvMass_VZgamma_GEN", "InvMass_VZgamma_GEN;m_{lljj#gamma} [GeV];Events", 150, 100., 1500., m_lljjgamma, self.weight)


        return True

    #=================
    #SIGNAL REGION
    #=================
    def IsInSignalRegion(self, studyPlot):

        # Number of Jets before Lepton Cut
        self.hEvent.fill1D("nJet_BeforeLepCut_RECO", "nJet_BeforeLepCut_RECO; N_{j};", 15, -0.5, 14.5, len(self.Jets), self.weight)
        # Number of Leptons before Lepton Cut
        self.hEvent.fill1D("nLepton_BeforeLepCut_RECO", "nLepton_BeforeLepCut_RECO; N_{l};", 15, -0.5, 14.5, len(self.Leptons), self.weight)

        # -- TRIGGER -- #
        if not self.event.HLT_passZZ4l: return False
        self.fillCutFlow("SR","Trigger")
        
        # -- LEPTONS -- #
        if len(self.ZCands) == 0 or self.bestZIdx < 0: return False
        self.theZ = self.ZCands[self.bestZIdx]
        self.ZCandLeptons, self.IsZCand = self.ZCandMassCut(self.theZ)
        if not self.IsZCand: return False
        self.fillCutFlow("SR","ZCandMass")
        # Invariant Mass of the ZCand after cuts Distribution
        self.hEvent.fill1D("InvMass_ZCand_PreCut_RECO", "InvMass_ZCand_PreCut_RECO;m_{ll} [GeV];Events", 60, 0., 180., self.theZ.mass, self.weight)

        self.KinCutLeptons = self.leptonKinCut(self.ZCandLeptons, "Region")
        if len(self.KinCutLeptons) != 2: return False
        self.fillCutFlow("SR","LeptonNumKin")
        # Number of Jets after Lepton Cut
        self.hEvent.fill1D("nJet_AfterLepCut_RECO", "nJet_AfterLepCut_RECO; N_{j};", 15, -0.5, 14.5, len(self.Jets), self.weight)
        # Number of Leptons after Lepton Cut
        self.hEvent.fill1D("nLepton_AfterLepCut_RECO", "nLepton_AfterLepCut_RECO; N_{l};", 15, -0.5, 14.5, len(self.KinCutLeptons), self.weight)
        
        # -- PHOTONS -- #
        self.VetoPass, self.VetoPhotons = self.getVetoPhotons(self.Photons)
        if not self.VetoPass: return False
        self.fillCutFlow("SR","ElectronVeto")

        # Inspection on pt-ratio/DeltaR relation of Photons that pass Kin Cuts
        self.ptRatioDeltaRPhotonJetRelation(self.VetoPhotons, self.Jets, studyPlot)
        # Tuning of deltaR threshold to select the photon
        self.deltaRThresholdScan(self.VetoPhotons, self.KinCutLeptons, self.Jets, studyPlot)
        # Inspection on disambiguation photons-jets (before KinCut)
        self.studyPhotonJetDisambiguation(self.VetoPhotons, self.Jets, self.KinCutLeptons, studyPlot)

        self.KinCutPhotons = self.photonKinCut(self.VetoPhotons, "Region")
        if len(self.KinCutPhotons) == 0: return False
        self.fillCutFlow("SR","PhotonKin")

        # Tuning of deltaR threshold to select the photon
        self.NoLepPhotons = [p for p in self.KinCutPhotons if all(p.DeltaR(l) > 0.5 for l in self.KinCutLeptons)]
        self.NoLepJets = [j for j in self.Jets if all(j.DeltaR(l) > 0.5 for l in self.KinCutLeptons)]
        if len(self.NoLepPhotons) == 0 or len(self.NoLepJets) < 2: return False
        self.fillCutFlow("SR","FromLepIso")

        # Disambiguation Photons-Jets
        self.GoodPhotons, self.BadPhotons, self.GoodJets, self.BadJets = self.disambiguatePhotonJet(
            self.NoLepPhotons, self.NoLepJets, deltaRThr=0.1
        )
        if len(self.GoodPhotons) == 0 or len(self.GoodJets) < 2: return False
        self.fillCutFlow("SR","PhotonJetIso")

        self.IdPass, self.TruePhotons = self.getTruePhotons(self.GoodPhotons, IdType="MVA_wp90")
        if not self.IdPass: return False
        self.fillCutFlow("SR","PhotonId")

        for ph in self.TruePhotons:
            self.hEvent.fill1D("Pt_Photon_PostPhCut_RECO","Pt_Photon_PostPhCut_RECO;#pt_{T}^{#gamma};# photons", 50, 0., 150., ph.pt, self.weight)

        
        #self.IsoPhotons = self.photonIsolationCut(self.TruePhotons, self.KinCutLeptons, self.Jets)
        #if len(self.IsoPhotons) == 0: return False
        #self.fillCutFlow("SR","PhotonIso")

        # Number of Jets after Photon Cut
            self.hEvent.fill1D("nJet_AfterPhotonCut_RECO", "nJet_AfterPhotonCut_RECO; N_{j};", 15, -0.5, 14.5, len(self.Jets), self.weight)

        # -- JETS -- #
        for j in self.Jets:
            self.hEvent.fill1D("Eta_Jet_PostPhCut_RECO","Eta_Jet_PostPhCut_RECO;#eta_{j};# jets", 50, -5., 5., j.eta, self.weight)       
        for j in self.Jets:
            self.hEvent.fill1D("Pt_Jet_PostPhCut_RECO","Pt_Jet_PostPhCut_RECO;p_{T}^{j};# jets", 50, 0., 150., j.pt, self.weight)
                
        self.KinCutJets = self.jetKinCut(self.GoodJets)
        if len(self.KinCutJets) < 2: return False
        self.fillCutFlow("SR","JetKin")

        for j in self.KinCutJets:
            self.hEvent.fill1D("Eta_Jet_PostJetKinCut_RECO","Eta_Jet_PostJetKinCut_RECO;#eta_{j};# jets", 50, -5., 5., j.eta, self.weight)
        for j in self.KinCutJets:
            self.hEvent.fill1D("Pt_Jet_PostJetKinCut_RECO","Pt_Jet_PostJetKinCut_RECO;p_{T}^{j};# jets", 50, 0., 150., j.pt, self.weight)
                
        #self.IsoJets = self.jetIsolationCut(self.KinCutJets,self.KinCutLeptons,self.TruePhotons)
        #if len(self.IsoJets) < 2: return False
        #self.fillCutFlow("SR","JetIso")        

        self.VCandJets = self.getVJets(self.KinCutJets)
        if len(self.VCandJets) != 2: return False
        self.fillCutFlow("SR","JetVCand")

        # Number of Jets after Jet Cut
        self.hEvent.fill1D("nJet_AfterJetCut_RECO", "nJet_AfterJetCut_RECO; N_{j};", 15, -0.5, 14.5, len(self.VCandJets), self.weight)

        
        
        
        
        
        # Invariant Mass of the ZCand (2l) after cuts Distribution
        self.hEvent.fill1D("InvMass_ZCand_RECO", "InvMass_ZCand_RECO;m_{ll} [GeV];Events", 60, 0., 180., self.theZ.mass, self.weight)        
        # Invariant Mass of 2l Gamma Distribution
        leading_photon = max(self.TruePhotons, key=lambda ph: ph.pt)
        m_llgamma = self.getInvMass(self.theZ, leading_photon)
        self.hEvent.fill1D("InvMass_llgamma_RECO", "InvMass_llgamma_RECO;m_{ll#gamma} [GeV];Events", 60, 60., 240., m_llgamma, self.weight)
        # Invariant Mass of 2l Gamma Distribution vs Invariant Mass of 2l Distribution
        self.hEvent.fill2D("InvMass_ll_llgamma_RECO", "InvMass_ll_llgamma_RECO;m_{ll} [GeV];m_{ll#gamma} [GeV];Events", 36, 0., 180., 50, 0., 250., self.theZ.mass, m_llgamma, self.weight)
        # Invariant Mass Distribution of the two best Jet candidate to form V-Boson
        m_jj = self.getInvMass(self.VCandJets[0], self.VCandJets[1])
        self.hEvent.fill1D("InvMass_VCandJets_RECO", "InvMass_VCandJets_RECO;m_{jj} [GeV];Events", 50, 40., 140., m_jj, self.weight)        
        # Invariant Mass of 2j Gamma Distribution
        m_jjgamma = self.getInvMass(self.VCandJets[0], self.VCandJets[1], leading_photon)
        self.hEvent.fill1D("InvMass_jjgamma_RECO", "InvMass_jjgamma_RECO;m_{jj#gamma} [GeV];Events", 60, 60., 240., m_jjgamma, self.weight)
        # Invariant Mass of 2j Gamma Distribution vs Invariant Mass of 2j Distribution
        self.hEvent.fill2D("InvMass_jj_jjgamma_RECO", "InvMass_jj_jjgamma_RECO;m_{jj} [GeV];m_{jj#gamma} [GeV];Events", 36, 0., 180., 50, 0., 250., m_jj, m_jjgamma, self.weight)
        # Invariant Mass of VZ Distribution
        m_VZ = self.getInvMass(self.theZ, self.VCandJets[0], self.VCandJets[1])
        self.hEvent.fill1D("InvMass_VZCand_RECO", "InvMass_VZCand_RECO;m_{ll} [GeV];Events", 100, 100., 1000., m_VZ, self.weight)
        # Invariant Mass of 2j 2l Gamma Distribution
        m_lljjgamma = self.getInvMass(self.theZ, leading_photon, self.VCandJets[0], self.VCandJets[1])
        self.hEvent.fill1D("InvMass_VZgamma_RECO", "InvMass_VZgamma_RECO;m_{2l2j#gamma} [GeV];Events", 150, 100., 1500., m_lljjgamma, self.weight)
        

        ## Cuts for Control Region
        #if self.theZ.mass < 80: return False
        #self.fillCutFlow("SR","ZCandMass100")

        # Vedere pt-ratio e deltaR per quei fotoni che non passano Id o KinCut, disambiguazione da jet
        # Vedere dopo aver fatto il taglio JetIso, fare loop sui fotoni e vedere quali hanno passato entrambi i tagli,
        # solo un taglio e entrambi i tagli (questo ultimo mi aspetto che siano tutti i fotoni)
        








        return True
  
            
    
    # ========================
    # USEFUL FUNCTION: LEPTONS         
    # ========================  
        
    def leptonKinCut(self,Leptons, Level):
        if Level == "Definition":
            baseLeptons = [p for p in Leptons if p.pt>10 and abs(p.eta)<2.5]
            baseLeptons.sort(key=lambda p: p.pt, reverse=True)
            if len(baseLeptons) >= 2 and baseLeptons[0].pt > 20:
                return baseLeptons
            else:
                return []
        elif Level == "Region":
            baseLeptons = [p for p in Leptons if p.pt>10 and (abs(p.pdgId)==11 and abs(p.eta)<2.5 or abs(p.pdgId)==13 and abs(p.eta)<2.4)]
            baseLeptons.sort(key=lambda p: p.pt, reverse=True)
            if len(baseLeptons) >= 2 and baseLeptons[0].pt > 20:
                return baseLeptons
            else:
                return []

    # GEN
    def getZLeptons(self, Leptons):
        bestLeptonPair = []
        min_diff = 40.0
        ZMassPDG = TDatabasePDG.Instance().GetParticle(23).Mass() 

        for l1, l2 in itertools.combinations(Leptons, 2):
            if l1.pdgId == -l2.pdgId:
                m_ll = self.getInvMass(l1, l2)
                if 60.0 < m_ll < 120.0:
                    diff = abs(m_ll - ZMassPDG)
                    if diff < min_diff:
                        min_diff = diff
                        bestLeptonPair = [l1, l2]

        return bestLeptonPair

    # RECO
    def ZCandMassCut(self, ZCand):
            if not 60 < ZCand.mass < 120 : return [], False
            return [self.Leptons[ZCand.l1Idx], self.Leptons[ZCand.l2Idx]], True
    
    # ========================
    # USEFUL FUNCTION: PHOTONS         
    # ========================

    def photonKinCut(self, Photons, Level):
        if Level == "Definition":
            CutPhotons = [p for p in Photons if p.pt>20 and abs(p.eta)<2.5]   
            return CutPhotons
        elif Level == "Region":
            CutPhotons = [p for p in Photons if p.pt>20 and abs(p.eta)<2.5 and not 1.442 < abs(p.eta) < 1.566]
            #NotCutPhotons = [p for p in Photons if p not in CutPhotons]
            return CutPhotons

    # GEN
    def photonIsolationCut(self, Photons, Leptons, Jets, deltaRThreshold=0.1):
        isolatedPhotons = []
        for p in Photons:
            lepIso = all(p.DeltaR(l) > 0.5 for l in Leptons)
            jetIso = True
            for j in Jets:
                if p.DeltaR(j) < 0.5:
                    dR = p.DeltaR(j)
                    if dR > deltaRThreshold: 
                        jetIso = False
                        break
            if lepIso and jetIso:
                isolatedPhotons.append(p)
        return isolatedPhotons

    # RECO
    def getVetoPhotons(self, Photons):
        if len(Photons) < 1: return False, None
        VetoPhotons = [ph for ph in Photons if (ph.cutBased >= 1 and ph.electronVeto)]
        if len(VetoPhotons) < 1: return False, None
        return True, VetoPhotons

    # RECO
    def disambiguatePhotonJet(self, Photons, Jets, deltaRThr=0.1):
        survivingJets = []
        discardedJets = []
        for j in Jets:
            if any(j.DeltaR(p) < deltaRThr for p in Photons):
                discardedJets.append(j)
            else:
                survivingJets.append(j)
                
        survivingPhotons = []
        discardedPhotons = []
        for p in Photons:
            if any(p.DeltaR(j) < 0.5 for j in survivingJets):
                discardedPhotons.append(p)
            else:
                survivingPhotons.append(p)
                
        return survivingPhotons, discardedPhotons, survivingJets, discardedJets

    '''# RECO
    def disambiguatePhotonJet(self, Photons, Jets, deltaRThr=0.1):
        survivingPhotons = []
        discardedPhotons = []
        jetIsDiscarded = [False] * len(Jets)

        for p in Photons:
            isPhotonDiscarded = False
            min_dR = 999.0
            closestJIdx = -1
            for jIdx, j in enumerate(Jets):
                dR = p.DeltaR(j)
                if dR < min_dR:
                    min_dR = dR
                    closestJIdx = jIdx
            if closestJIdx != -1 and min_dR < 0.5:
                if min_dR < deltaRThr:
                    jetIsDiscarded[closestJIdx] = True
                else:
                    isPhotonDiscarded = True
            if isPhotonDiscarded:
                discardedPhotons.append(p)
            else:
                survivingPhotons.append(p)

        survivingJets = []
        discardedJets = []
        for jIdx, j in enumerate(Jets):
            if not jetIsDiscarded[jIdx]:
                survivingJets.append(j)
            else:
                discardedJets.append(j)

        return survivingPhotons, discardedPhotons, survivingJets, discardedJets
    '''
    # RECO
    def getTruePhotons(self, Photons, IdType="MVA_wp90"):
        if IdType == "MVA_wp90":
            TruePhotons = [ph for ph in Photons if (ph.mvaID_WP90)]
        elif IdType == "Loose":
            TruePhotons = [ph for ph in Photons if (ph.cutBased >= 1)]
        elif IdType == "Medium":
            TruePhotons = [ph for ph in Photons if (ph.cutBased >= 2)]
        elif IdType == "Tight":
            TruePhotons = [ph for ph in Photons if (ph.cutBased >= 3)]
        if len(TruePhotons) < 1: return False, None
        return True, TruePhotons

    # INSPECTION
    def ptRatioDeltaRPhotonJetRelation(self, Photons, Jets, flag):
        if flag:
            for p in Photons:
                cat = self.categorizePhoton(p)
                self.hEvent.fill1D("Photon_Categories_RECO", "Photon Categories;0=Both, 1=KinOnly, 2=IdOnly, 3=None;Photons", 4, -0.5, 3.5, cat, self.weight)
                #if cat < 2:
                for j in Jets:
                    dR = p.DeltaR(j)
                    self.hEvent.fill1D("DeltaR_Overlap_Photons", 
                                        "DeltaR_Overlap_Photons;#DeltaR(#gamma, j);Events", 
                                        50, 0.05, 5.05, dR, self.weight)
                    pt_ratio = p.pt / j.pt
                    self.hEvent.fill2D("DeltaR_vs_PtRatio_Photons", 
                                        "DeltaR_vs_PtRatio_Photons;#DeltaR(#gamma, j);p_{T}^{#gamma}/p_{T}^{j};Events", 
                                        50, 0., 5.0, 60, 0., 2.0, dR, pt_ratio, self.weight)
                    if dR < 1.5:
                        self.hEvent.fill2D("DeltaR_vs_PtRatio_Photons_Zoom", 
                                           "DeltaR_vs_PtRatio_Photons_Zoom;#DeltaR(#gamma, j);p_{T}^{#gamma}/p_{T}^{j};Events", 
                                           20, 0., 0.5, 60, 0., 2.0, dR, pt_ratio, self.weight)
                    if dR < 0.1:
                        self.hEvent.fill1D("PtRatio_Overlap_GoodPhotons", 
                                           "PtRatio_Overlap_GoodPhotons (dR < 0.1);p_{T}^{#gamma}/p_{T}^{j};Events", 
                                           50, 0., 2.0, pt_ratio, self.weight)
                    if 0.1 < dR < 0.5:
                        self.hEvent.fill1D("PtRatio_Overlap_BadPhotons", 
                                           "PtRatio_Overlap_BadPhotons (0.1 < dR < 0.5);p_{T}^{#gamma}/p_{T}^{j};Events", 
                                           50, 0., 2.0, pt_ratio, self.weight)

            return True
        else:
            return False
        
    # INSPECTION
    def deltaRThresholdScan(self, Photons, Leptons, Jets, flag):
        if flag:
            thresholdsTest = [0.02, 0.04, 0.06, 0.08, 0.1, 0.12, 0.14, 0.16, 0.18]
            for thr in thresholdsTest:
                testPhotons = self.photonIsolationCut(Photons, Leptons, Jets, deltaRThreshold=thr)
                if len(testPhotons) > 0:
                    self.hEvent.fill1D("Tuning_DeltaR_Threshold", "Tuning_DeltaR_Threshold;Thr DeltaR(j,#gamma);SR Events", 
                                       10, 0.01, 0.21, thr, self.weight)
            return True
        else:
            return False

    # INSPECTION
    def studyPhotonJetDisambiguation(self, vetoPhotons, jets, kinLeptons, flag):
        if flag:
            """
            Studies the disambiguation Photons-Jets on photons that passed Electron Veto
            & classifies the results depending on KinCut & ID
            """
            NoLepPhotons = [p for p in vetoPhotons if all(p.DeltaR(l) > 0.5 for l in kinLeptons)]
            NoLepJets = [j for j in jets if all(j.DeltaR(l) > 0.5 for l in kinLeptons)]

            for j in NoLepJets:
                if len(NoLepPhotons) > 0:
                    min_dR = min(j.DeltaR(p) for p in NoLepPhotons)
                    self.hEvent.fill1D("MinDeltaR_Jet_Photon", "Min #DeltaR(j, #gamma);min #DeltaR;Jets", 
                                       100, 0., 5., min_dR, self.weight)
                    
                    self.hEvent.fill1D("MinDeltaR_Jet_Photon_Zoom", "Min #DeltaR(j, #gamma) (Zoom);min #DeltaR;Jets", 
                                       50, 0., 0.5, min_dR, self.weight)

            GoodPhotons, BadPhotons, GoodJets, BadJets = self.disambiguatePhotonJet(NoLepPhotons, NoLepJets, deltaRThr=0.1)

            nGoodPhotonsPassKin = 0
            for p in GoodPhotons:
                cat = self.categorizePhoton(p)
                self.hEvent.fill1D("Disambig_GoodPhotons_Categs", "Good Photons;0=Both, 1=KinOnly, 2=IdOnly, 3=None;Photons", 4, -0.5, 3.5, cat, self.weight)
                if cat in (0, 1):
                    nGoodPhotonsPassKin += 1

            nBadPhotonsPassKin = 0
            for p in BadPhotons:
                cat = self.categorizePhoton(p)
                self.hEvent.fill1D("Disambig_BadPhotons_Categs", "Bad Photons;0=Both, 1=KinOnly, 2=IdOnly, 3=None;Photons", 4, -0.5, 3.5, cat, self.weight)
                if cat in (0, 1):
                    nBadPhotonsPassKin += 1


            self.hEvent.fill1D("nGoodPhotons_Disambig_perEvent", "Good Photons per Event;# Photons;Events", 5, -0.5, 4.5, len(GoodPhotons), self.weight)
            self.hEvent.fill1D("nBadPhotons_Disambig_perEvent", "Bad Photons per Event;# Photons;Events", 5, -0.5, 4.5, len(BadPhotons), self.weight)

            self.hEvent.fill1D("nGoodPhotons_PassKin_perEvent", "Good Photons Pass Kin per Event;# Photons;Events", 5, -0.5, 4.5, nGoodPhotonsPassKin, self.weight)
            self.hEvent.fill1D("nBadPhotons_PassKin_perEvent", "Bad Photons Pass Kin per Event;# Photons;Events", 5, -0.5, 4.5, nBadPhotonsPassKin, self.weight)

            for j in jets:
                if j in GoodJets:
                    k = 1
                elif j in BadJets:
                    k = 2
                else:
                    k = 0
                self.hEvent.fill1D("nJet_AfterDisambiguation", "Jet Status;0=LepOverlap, 1=PassDisambig, 2=FailDisambig;Jets", 3, -0.5, 2.5, k, self.weight)

            self.hEvent.fill1D("nGoodJets_Disambig_perEvent", "Passed Jets per Event;# Jets;Events", 15, -0.5, 14.5, len(GoodJets), self.weight)
            self.hEvent.fill1D("nBadJets_Disambig_perEvent", "Discarded Jets per Event;# Jets;Events", 10, -0.5, 9.5, len(BadJets), self.weight)

            #return GoodPhotons, BadPhotons, GoodJets, BadJets
            return True
        else:
            return False

    # ========================
    # USEFUL FUNCTION: JETS         
    # ========================

    def jetKinCut(self,Jets):
        kinJet = [j for j in Jets if j.pt>30 and abs(j.eta)<4.7]
        return kinJet

    def jetIsolationCut(self,Jets,Leptons,Photons):
        jetIso = [j for j in Jets if all(j.DeltaR(l) > 0.5 for l in Leptons)
                                    and all(j.DeltaR(p) > 0.5 for p in Photons)]
        return jetIso

    def getVJets(self, Jets):
        bestJetPair = []
        global_min_diff = 50.0
        WMassPDG = TDatabasePDG.Instance().GetParticle(24).Mass()
        ZMassPDG = TDatabasePDG.Instance().GetParticle(23).Mass()

        for j1, j2 in itertools.combinations(Jets, 2):
            m_jj = self.getInvMass(j1, j2)
            if 50.0 < m_jj < 120.0:
                diff_W = abs(m_jj - WMassPDG)
                diff_Z = abs(m_jj - ZMassPDG)
                current_min_diff = min(diff_W, diff_Z)
                if current_min_diff < global_min_diff:
                    global_min_diff = current_min_diff
                    bestJetPair = [j1, j2]

        return bestJetPair

    # ====================
    # OTHER FUNCTIONS
    # ====================

    def categorizePhoton(self, p):
        pass_kin = (p.pt > 20 and abs(p.eta) < 2.5 and not 1.442 < abs(p.eta) < 1.566)
        pass_id = getattr(p, "mvaID_WP90", False)
        if pass_kin and pass_id: return 0
        elif pass_kin and not pass_id: return 1
        elif not pass_kin and pass_id: return 2
        else: return 3

    def PhotonJetDeltaRCut(self,Photons,Jets):
        deltaR = [p.DeltaR(j) for p in Photons for j in Jets]
        return deltaR

    def PhotonTopJetsDeltaR(self, Photons, Jets):
        deltaR_list = []
        top_jets = Jets[:3] 
        
        for p in Photons:
            for j in top_jets:
                deltaR_list.append(p.DeltaR(j))
                
        return deltaR_list

    def getInvMass(self, *particles):
        # Given a variable number of "daughter particles" it calculates the invariant mass
        ptot = ROOT.TLorentzVector()
        for part in particles:
            p4_temp = ROOT.TLorentzVector()
            
            pt   = getattr(part, "pt", 0.0)
            eta  = getattr(part, "eta", 0.0)
            phi  = getattr(part, "phi", 0.0)
            mass = getattr(part, "mass", 0.0)
            
            p4_temp.SetPtEtaPhiM(pt, eta, phi, mass)
            ptot += p4_temp
            
        return ptot.M()
        #for i, part in enumerate(particles):
        #    if i == 2 and not hasattr(part, "pdgId") :  
        #        p4_temp = ROOT.TLorentzVector()
        #        p4_temp.SetPtEtaPhiM(part.pt, part.eta, part.phi, 0.0)
        #        ptot += p4_temp
        #    else:
        #        ptot += part.p4()              
        #return ptot.M()

    def fillCutFlow(self, region, cut_label):
        # param region: "SD" (Signal Definition) o "SR" (Signal Region)
        if region == "SD":
            self.hEvent.fill1D_label("CutFlowSD", "CutFlowSD", self.CutFlowLabelsSD, cut_label, self.weight)
        elif region == "SR":
            self.hEvent.fill1D_label("CutFlowSR", "CutFlowSR", self.CutFlowLabelsSR, cut_label, self.weight)

# OTHER CLASSES

#class DressedLepton :
#    def __init__(self, nakedLepton, Photon = None ) :
#        self.nakedLepton = nakedLepton
#        self.pdgId = nakedLepton.pdgId
#        self.charge = nakedLepton.charge
#        self.fsrPhoton = Photon
#
#        p4_total = nakedLepton.p4()
#        if Photon is not None:
#            p4_total += Photon.p4()
#
#        self._p4 = p4_total
#        #kinematic
#        self.pt = self._p4.Pt()     
#        self.eta = self._p4.Eta()
#        self.phi = self._p4.Phi()
#        self.mass = self._p4.M()
#
#    def p4(self):
#        return self._p4
                  
class GenStatusFlag(IntFlag):
    IS_PROMPT = 1 << 0
    FROM_HARD_PROCESS = 1 << 8
    FROM_HARD_PROCESS_BEFORE_FSR = 1 << 11
    IS_FIRST_COPY = 1 << 12
    IS_LAST_COPY = 1 << 13

    SEL = IS_PROMPT | FROM_HARD_PROCESS
    QUARK_SEL = IS_PROMPT | FROM_HARD_PROCESS | IS_FIRST_COPY

