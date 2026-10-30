from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection
from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Object

from VVXAnalysis.NanoAnalysis.EventAnalyzer import EventAnalyzer
from VVXAnalysis.NanoAnalysis.Histogrammer import *
from VVXAnalysis.NanoAnalysis.Regions import Flags as Regions

import ROOT
import math
from ROOT import TDatabasePDG  
from enum import IntFlag

"I commented out all pieces of code that are not useful when running the analysis on background samples"

class WZGammaAnalyzer(EventAnalyzer, analysis_name="WZGammaAnalyzer"):
     
    def __init__(self, regions):
        super().__init__(regions)
        
   
    # creates particle collections
    def loadCollections(self) :
        if self.analyzeMC:  # gen-level
            self.GenParts = Collection(self.event, 'GenPart')
            self.GenChargedLeptons = [p for p in self.GenParts if abs(p.pdgId) in (11,13) and p.status == 1
                                and (p.statusFlags & GenStatusFlag.SEL) == GenStatusFlag.SEL]
            self.GenPhotons = [p for p in self.GenParts if p.pdgId == 22 and p.status == 1
                              and (p.statusFlags & GenStatusFlag.SEL) == GenStatusFlag.SEL]
            #self.GenPhotons_Prompt = [p for p in self.GenParts if p.pdgId == 22 and p.status == 1 and (p.statusFlags & (1<<0)) != 0]
            self.GenNeutrinos = [p for p in self.GenParts if abs(p.pdgId) in (12,14)  and p.status == 1
                            and (p.statusFlags & GenStatusFlag.SEL) == GenStatusFlag.SEL]

        # reco-level
        self.event.SetBranchStatus("*PuppiMET*", 1)
        self.event.SetBranchStatus("*PFMET*", 1)
        self.event.SetBranchStatus("*ZCand*", 1)
        self.event.SetBranchStatus("*FsrPhoton*", 1)
        
        self.Electrons = Collection(self.event, 'Electron')
        self.Muons = Collection(self.event, 'Muon')
        self.ChargedLeptons = list(self.Electrons) + list(self.Muons)
        self.Photons = Collection(self.event, 'Photon')
        self.ZCands = Collection(self.event, 'ZCand')          
        self.bestZIdx = self.event.bestZIdx
        self.PFMET = Object(self.event, "PFMET")    
        self.PuppiMET = Object(self.event, "PuppiMET")
        self.FSRPhotons = Collection(self.event, 'FsrPhoton')


    def analyze(self):  
        self.loadCollections()
        #self.SDCuts = ["Weighted_Events","Lep_KinCuts","WZ_MassLimits","Ph_KinCuts", "Ph_AngSeparation", "FSR_onZ", "FSR_onW"]
        #self.SRCuts=["Weighted_Events","HLT_passZZ4l","Lep_KinCuts","MET_Cut","Z_MassLimits","ZZFullSel","Ph_kinCuts",
                     "Ph_electronVeto","Ph_mvaWP80","Ph_AngSeparation", "FSR_onZ", "FSR_onW"]

        self.isSig = False
        if self.analyzeMC :
            self.isSig = self.isSignal()  
            genW, genZ = self.genCouple

        self.isInSR = self.isInSignalRegion(2024)
        recoZ, recoW = self.recoCouple                

        '''# sample's statistic values
        self.getPercentagePlot("Background_inSR")
        self.getPercentagePlot("ReconstructonEfficency")'''

            
        '''if isSig and isInSR :             
            # GEN-RECO LEPTONS MATCHING
            res = self.isRecoMatching(genZ.leptons, genW.leptons, (self.ChargedLeptons[recoZ.l1Idx], self.ChargedLeptons[recoZ.l2Idx]), recoW.lepton, 0.5)
            self.hEvent.fill1D("GenAndReco_Matching", "GenAndReco_Matching", 2, -0.5, 1.5, int(res), self.weight)'''    

        '''if self.hasGenFSR and self.lGammaMET_mt != None :
            self.hEvent.fill1D("lGammaMET_Mt_Cutted", "lGammaMET_Mt_Cutted", 30, 0., 300., self.lGammaMET_mt, self.weight)'''
   

    #---------------------      
    #SIGNAL DEFINITION
    #---------------------
    def isSignal(self) :        
        self.genCouple = None, None
        self.hasGenFSR = False

        self.CutPlot("Definition", 0)  
        GenBestLeptons = self.pass_ChLepKinCut(self.GenChargedLeptons, "Definition")  
        if not GenBestLeptons : return False
        self.CutPlot("Definition", 1)
        
        '''# COMPARISON OF WZ RECONSTRUCTION ALGORITHMS           
        if self.isEventEasy(self.GenBestLeptons) :        
            self.methodOutcome(False, False,"MethodZ_Outcome")  #Z-first
            self.methodOutcome(False, True, "MethodW_Outcome")   #W-first
            self.methodOutcome(True, False, "MethodWZ_Outcome")   #residues'''
       
        # W and Z RECONSTRUCTION (with the Z-first method)
        theW, theZ = self.WZRecon(GenBestLeptons, self.GenNeutrinos, False, False, False)
        if theW is None or theZ is None : return False
        self.CutPlot("Definition", 2)
        self.genCouple = (theW, theZ)
        lepW = next((l for l in theW.leptons if abs(l.pdgId) in (11,13)), None)
        nu  = next((n for n in theW.leptons if abs(n.pdgId) in  (12,14)), None)

        self.hEvent.fill1D("Gen_nPhoton", "Gen_nPhoton",11, -0.5, 10.5, len(self.GenPhotons), self.weight)
        GoodGenPhoton = self.photonsCuts(self.GenPhotons, self.GenChargedLeptons, "Definition")
        if not GoodGenPhoton : return False 

        '''# for graphic purposes only
        photon_maxPt = max(GoodGenPhoton, key=lambda p: p.pt, default=None) 
        Gen_invMassLLGamma_preCuts = self.getInvMass(theZ.leptons[0], theZ.leptons[1], photon_maxPt)
        Gen_invMassLNuGamma_preCuts = self.getInvMass(theW.leptons[0], theW.leptons[1], photon_maxPt)
        self.hEvent.fill2D("llG_vs_ll_invMass", "llG_vs_ll_invMass", 100, 40., 240., 30, 60., 120.,
Gen_invMassLLGamma_preCuts, theZ.mass, self.weight)
        self.hEvent.fill2D("lNuG_vs_lNu_invMass", "lNuG_vs_lNu_invMass", 100, 40., 240., 30, 50., 110., Gen_invMassLNuGamma_preCuts, theW.mass, self.weight)
        self.hEvent.fill1D("Gen_invMassLNuGamma", "Gen_invMassLNuGamma", 30, 0., 300., Gen_invMassLNuGamma_preCuts, self.weight)
         #computation of different transverse mass for comparison with reco-level
        if (nu.pt >= 30.0 ) :  #same condition as at reco-level
            W_mt_Gen = self.computeMt([(lepW.pt, lepW.eta, lepW.phi)], nu.pt, nu.phi)
            
            lGamma_list = [(lepW.pt, lepW.eta, lepW.phi), (photon_maxPt.pt, photon_maxPt.eta, photon_maxPt.phi)]
            lGammaMET_mt_Gen = self.computeMt(lGamma_list, nu.pt, nu.phi)
            
            lZ_list = [(lepW.pt, lepW.eta, lepW.phi), (theZ.leptons[0].pt, theZ.leptons[0].eta, theZ.leptons[0].phi),(theZ.leptons[1].pt, theZ.leptons[1].eta, theZ.leptons[1].phi)]
            lZMET_mt_Gen = self.computeMt(lZ_list, nu.pt, nu.phi)
            
            lZGamma_list = [(lepW.pt, lepW.eta, lepW.phi), (theZ.leptons[0].pt, theZ.leptons[0].eta, theZ.leptons[0].phi),
(theZ.leptons[1].pt, theZ.leptons[1].eta, theZ.leptons[1].phi) , (photon_maxPt.pt, photon_maxPt.eta, photon_maxPt.phi)]
            lZGammaMET_mt_Gen = self.computeMt(lZGamma_list, nu.pt, nu.phi)

            self.hEvent.fill1D("Gen_W_Mt", "Gen_W_Mt", 30, 0., 300., W_mt_Gen, self.weight)
            self.hEvent.fill1D("Gen_lGammaMET_Mt", "Gen_lGammaMET_Mt", 30, 0., 300., lGammaMET_mt_Gen, self.weight)
            self.hEvent.fill1D("Gen_lZGammaMET_Mt", "Gen_lZGammaMET_Mt", 50, 100., 600., lZGammaMET_mt_Gen, self.weight)
                        
        if not self.photonIsNotFSR(theW.leptons, photon_maxPt, 90) and nu.pt >= 30.0 :
            self.hEvent.fill1D("lGammaMET_Mt_CuttedGen", "lGammaMET_Mt_CuttedGen", 30, 0., 300., lGammaMET_mt_Gen, self.weight)
            self.hasGenFSR = True'''
        
        SurvivedGenPhotons = [gp for gp in GoodGenPhoton if self.photonIsNotFSR(theZ.leptons, gp, 100)]
        if not SurvivedGenPhotons : return False
        self.CutPlot("Definition", 5)  #FSR_onZ

        UltimateGenPhotons = [sp for sp in SurvivedGenPhotons if self.photonIsNotFSR(theW.leptons, sp, 90)]
        if not UltimateGenPhotons : return False
        else : genPhoton = max(UltimateGenPhotons, key=lambda p: p.pt, default=None)
        self.CutPlot("Definition", 6)  #FSR_onW

        invMassLLGamma_Gen = self.getInvMass(theZ.leptons[0], theZ.leptons[1], genPhoton)
        invMassLNuGamma_Gen = self.getInvMass(theW.leptons[0], theW.leptons[1], genPhoton)
        
        
        '''self.hEvent.fill2D("W_vs_lGamma_Mt_Gen", "W_vs_lGamma_Mt_Gen", 30, 0., 300., 12, 0., 120., lGammaMET_mt_Gen,
        W_mt_Gen, self.weight)
        self.hEvent.fill1D("Gen_invMassLLGamma", "Gen_invMassLLGamma", 30, 0., 300., invMassLLGamma_Gen, self.weight)            
        self.hEvent.fill1D("Gen_ZMass_91Gev", "ZGenMass_91Gev", 120, 89., 93, theZ.mass, self.weight)
        self.hEvent.fill1D("Gen_WMass_80Gev", "WGenMass_80Gev", 160, 50., 110., theW.mass, self.weight)
        '''

        return True
                                 
   
    #-----------------------
    # SIGNAL REGION
    #-----------------------
    def isInSignalRegion(self, year) :
        self.recoCouple = (None, None)      
        self.MET = None
        self.lGammaMET_mt = None

        self.CutPlot("Region", 0)  #Weighted Events
        if not self.event.HLT_passZZ4l : return False
        self.CutPlot("Region", 1)  #HLT_passZZ4l
        
        BestLeptons = self.pass_ChLepKinCut(self.ChargedLeptons, "Region")
        if not BestLeptons : return False
        self.CutPlot("Region", 2)  #Lep_KinCuts

        puppyYears = [2024]
        self.MET = self.PuppiMET if year in puppyYears else self.PFMET
        if self.MET.pt < 30 : return False
        self.CutPlot("Region", 3)  #MET_Ptcut
                                          
        theZ = self.ZCands[self.bestZIdx]
        if 60. < theZ.mass < 120. :
            lep1, lep2 = self.ChargedLeptons[theZ.l1Idx], self.ChargedLeptons[theZ.l2Idx]
            lepW_nd = next((l for l in BestLeptons if l not in (lep1, lep2)), None)
            #dressed leptons
            if theZ.fsr1Idx > -1 : lepZ1 = DressedLepton(lep1, self.FSRPhotons[theZ.fsr1Idx])
            else : lepZ1 = lep1
            if theZ.fsr2Idx > -1 : lepZ2 = DressedLepton(lep2, self.FSRPhotons[theZ.fsr2Idx])
            else : lepZ2 = lep2
        else : return False
        self.CutPlot("Region", 4)

        if not lepW_nd.ZZFullSel : return False
        self.CutPlot("Region", 5)

        W_mt = self.computeMt([(lepW_nd.pt, lepW_nd.eta, lepW_nd.phi)], self.MET.pt, self.MET.phi)
        if (lepW_nd.fsrPhotonIdx != -1 and lepW_nd.fsrPhotonIdx not in (theZ.fsr1Idx, theZ.fsr2Idx)) :
            fsrPhoton = self.FSRPhotons[lepW_nd.fsrPhotonIdx]   
            lepW = DressedLepton(lepW_nd, fsrPhoton)
            W_mt = self.computeMt([(lepW.pt, lepW.eta, lepW.phi)], self.MET.pt, self.MET.phi)
            theW = WReco(lepW, W_mt, fsrPhoton)
        else : 
            lepW = lepW_nd
            theW = WReco(lepW, W_mt)            
        self.recoCouple = (theZ, theW)
        
        '''BestPhotons = self.photonsCuts(self.Photons, self.ChargedLeptons, "Region")
        if not BestPhotons : return False

        photon_maxPt = max(BestPhotons, key=lambda p: p.pt, default=None)  # for graphic purposes only
        invMassLLGamma_preCuts = self.getInvMass(lepZ1, lepZ2, photon_maxPt)
        lGamma_list_pC = [(lepW.pt, lepW.eta, lepW.phi), (photon_maxPt.pt, photon_maxPt.eta, photon_maxPt.phi)]
        lGammaMET_mt_preCuts = self.computeMt(lGamma_list_pC, self.MET.pt, self.MET.phi)
        self.lGammaMET_mt = lGammaMET_mt_preCuts
        self.hEvent.fill1D("llGamma_invMass_preCuts", "llGamma_invMass_preCuts", 30, 0., 300., invMassLLGamma_preCuts, self.weight)
        self.hEvent.fill1D("lGammaMET_Mt_preCuts", "lGammaMET_Mt_preCuts", 30, 0., 300., lGammaMET_mt_preCuts, self.weight)
        self.hEvent.fill2D("W_vs_lGammaMET_Mt", "W_vs_lGammaMET_Mt_Gen", 30, 0., 300., 12, 0., 120., lGammaMET_mt_preCuts,
        W_mt, self.weight)
        SurvivedPhotons = [bp for bp in BestPhotons if self.photonIsNotFSR([lepZ1,lepZ2], bp, 100)]
        if not SurvivedPhotons : return False
        self.CutPlot("Region", 10)  #FSR_onZ
        
       
        UltimatePhotons = [sp for sp in SurvivedPhotons if self.photonIsNotWFSR(lepW, sp, 90)]
        if not UltimatePhotons : return False
        else : thePhoton = max(UltimatePhotons, key=lambda p: p.pt, default=None)
        self.CutPlot("Region", 11)  #FSR_onW

        invMassLLGamma = self.getInvMass(lepZ1, lepZ2, thePhoton)
        lGamma_list = [(lepW.pt, lepW.eta, lepW.phi), (thePhoton.pt, thePhoton.eta, thePhoton.phi)]
        lGammaMET_mt = self.computeMt(lGamma_list, self.MET.pt, self.MET.phi)
        
        lZ_list = [(lepW.pt, lepW.eta, lepW.phi), (lepZ1.pt, lepZ1.eta, lepZ1.phi), (lepZ2.pt, lepZ2.eta, lepZ2.phi)]
        lZMET_mt = self.computeMt(lZ_list, self.MET.pt, self.MET.phi)
        lZGamma_list = [(lepW.pt, lepW.eta, lepW.phi), (lepZ1.pt, lepZ1.eta, lepZ1.phi), (lepZ2.pt, lepZ2.eta, lepZ2.phi), (thePhoton.pt, thePhoton.eta, thePhoton.phi)]
        lZGammaMET_mt = self.computeMt(lZGamma_list, self.MET.pt, self.MET.phi)'''
                        
        self.hEvent.fill1D("LepZ_1_pt", "LepZ_1_pt", 20, 0., 200., lepZ1.pt, self.weight)
        self.hEvent.fill1D("LepZ_2_pt", "LepZ_2_pt", 20, 0., 200., lepZ2.pt, self.weight)
        self.hEvent.fill1D("LepW_pt", "LepW_pt", 20, 0., 200., lepW.pt, self.weight)
        self.hEvent.fill1D("LepW_eta", "LepW_eta", 12, -2.5, 2.5, lepW.eta, self.weight)
        self.hEvent.fill1D("MET_pt", "MET_pt", 17, 30., 200., self.MET.pt, self.weight)
        self.hEvent.fill1D("MET_phi", "MET_phi", 18, -3.2, 3.2, self.MET.phi, self.weight)
        #self.hEvent.fill1D("Photon_pt", "Photon_pt", 40, 0., 200., thePhoton.pt, self.weight)
        #self.hEvent.fill1D("Photon_eta", "Photon_eta", 12, -2.4, 2.4, thePhoton.eta, self.weight)
        self.hEvent.fill1D("ZCandMass_91Gev", "ZCandMass_91Gev", 60, 60., 120., theZ.mass, self.weight)
        self.hEvent.fill1D("W_Mt", "W_Mt", 10, 0., 200., W_mt, self.weight)
        #self.hEvent.fill1D("lGammaMET_Mt", "lGammaMET_Mt", 15, 0., 300., lGammaMET_mt, self.weight)
        #self.hEvent.fill1D("lZMET_Mt", "lZMET_Mt", 14, 0., 700., lZMET_mt, self.weight)
        #self.hEvent.fill1D("lZGammaMET_Mt", "lZGammaMET_Mt", 10, 0., 1000., lZGammaMET_mt, self.weight)

        return True
   
 
    #--------------------------------------
    # Functions for WZ RECOSTRUCTION
    #--------------------------------------
    # given a collection of leptons (charged and neutrino), this function reconstructs the W or the Z by the pdgId given
    def BosonRecon(self, GoodLeptons, partPdgId,  massThreshold, isForComparison) :
        CandMass = {}
        MassPdg = TDatabasePDG.Instance().GetParticle(partPdgId).Mass()

        if not isForComparison :
            for i in range(len(GoodLeptons)):
                for j in range(i + 1, len(GoodLeptons)):
                    isW = (self.getFamily(GoodLeptons[i].pdgId) == self.getFamily(GoodLeptons[j].pdgId)                                                  and self.getFamily(GoodLeptons[i].pdgId) != 0
                            and abs(GoodLeptons[i].pdgId + GoodLeptons[j].pdgId) == 1)
                    if (partPdgId == 24 and isW) :  
                        CandMass[(i,j)] = self.getInvMass(GoodLeptons[i], GoodLeptons[j])  
                    elif (partPdgId == 23 and GoodLeptons[i].pdgId + GoodLeptons[j].pdgId == 0) :
                        CandMass[(i,j)] = self.getInvMass(GoodLeptons[i], GoodLeptons[j])  
                    else : CandMass[(i,j)] = None

        #just for algorithm comparison, based only on charge conditions
        if isForComparison :
            for i in range(len(GoodLeptons)):
                for j in range(i + 1, len(GoodLeptons)):
                    if self.isChargeCompatible(GoodLeptons[i].pdgId, GoodLeptons[j].pdgId, partPdgId) :
                        CandMass[(i,j)] = self.getInvMass(GoodLeptons[i], GoodLeptons[j])
                    else : CandMass[(i,j)] = None
                           
        bestCouple = self.getCloserCand(MassPdg, CandMass)
        if(bestCouple != None and massThreshold[0] < CandMass[bestCouple] < massThreshold[1] or (isForComparison and bestCouple is not None)) :  
            if(partPdgId == 24) : boson = GenW([GoodLeptons[bestCouple[0]], GoodLeptons[bestCouple[1]]])
            elif(partPdgId == 23) : boson = GenZ([GoodLeptons[bestCouple[0]], GoodLeptons[bestCouple[1]]])
            return boson, bestCouple

        return None, None


    # return a W and a Z boson. Using 3 different algorithms based on the useResidues and recbyW options:
    # useResidues = True, recbyW = False for the residues method
    # useResidues = False, recbyW = True for the W-first method
    # useResidues = False, recbyW = False for the Z-first method
    def WZRecon(self, ChLeptons, Neutrinos, useResidues, recbyW, isForComparison) :
        if not any(l.pdgId < 0 for l in ChLeptons) or not any(l.pdgId > 0 for l in ChLeptons): return None, None         
        
        # for Z FIRST and RESIDUES METHOD
        Z1, coupleZ1 = (None, None) if recbyW else self.BosonRecon(ChLeptons, 23, (60,120), isForComparison)
        W1, coupleW1 = None, None
        if coupleZ1 is not None :
            RemainingChLeptons1 = [p for i, p in enumerate(ChLeptons) if i not in coupleZ1]
            Leptons1 = RemainingChLeptons1 + Neutrinos
            W1, coupleW1 = self.BosonRecon(Leptons1, 24, (50,110), isForComparison)     
        bosonCouple = (W1, Z1) if (coupleW1 is not None and coupleZ1 is not None) else (None, None)
        
        # for W FIRST and RESIDUES METHOD
        if useResidues or recbyW :
            Leptons2 = ChLeptons + Neutrinos
            W2, coupleW2 = self.BosonRecon(Leptons2, 24, (50,110),isForComparison)
            Z2, coupleZ2 = None, None
            if coupleW2 is not None :
                RemainingChLeptons2 = [p for i, p in enumerate(ChLeptons) if i not in coupleW2]
                Z2, coupleZ2 = self.BosonRecon(RemainingChLeptons2, 23, (60,120), isForComparison)
            couple2Valid = coupleW2 is not None and coupleZ2 is not None
            couple1Valid = bosonCouple[0] is not None and bosonCouple[1] is not None
            if couple2Valid and (not couple1Valid or self.getLowerRes(W1, Z1, W2, Z2)):
                bosonCouple = (W2, Z2)  
          
        return bosonCouple


    # Method that studies the algoritm results in terms of "error" (looks at the lepton coupling)
    def methodOutcome(self, useResidues, recbyW, histName) :
        bCouple, fail = self.WZRecon(self.GenBestLeptons, self.GenNeutrinos, useResidues, recbyW, True)
        theW, theZ = bCouple
        if theW is not None and theZ is not None : 
            if not self.pairIsTruthBoson(theZ.leptons, 23) or not self.pairIsTruthBoson(theW.leptons,                        24): outcome = 0   # error
            else: outcome = 1   # success
            self.hEvent.fill1D(histName, histName, 2, -0.5, 1.5, outcome, self.weight)

    
    #----------------------------------
    # Functions for CUTS
    #----------------------------------
    # if the cut is passed (and also the trigger 20-10-5) it returns the three charged leptons with the biggest pt
    def pass_ChLepKinCut(self, ChLeptons, level) :
        if level == "Definition" : GoodChLeptons = [p for p in ChLeptons if p.pt > 5 and abs(p.eta) < 2.5]
        elif level == "Region" : GoodChLeptons = [p for p in ChLeptons if (abs(p.pdgId) == 13 and p.pt > 5 or abs(p.pdgId) == 11 and p.pt > 7)  and  (abs(p.pdgId) == 13 and abs(p.eta) < 2.4 or abs(p.pdgId) == 11 and abs(p.eta) < 2.5)]
        else : raise ValueError(f"Not valid option: '{option}'.")
        if len(GoodChLeptons) > 2 :
            l1 = max(GoodChLeptons, key=lambda p: p.pt, default=None)
            if l1.pt > 20 :
                l2 = max((p for p in GoodChLeptons if p != l1), key=lambda p: p.pt, default=None)
                if l2.pt  > 10 :
                    l3 = max((p for p in GoodChLeptons if p not in (l1, l2)), key=lambda p: p.pt, default=None)
                    return [l1,l2,l3]                  
        return None

   
    # for kinematic cut on photons: returns the photons that passed all the cuts
    def photonsCuts(self, Photons, ChLeptons, Level) :
        GoodPhotons = [p for p in Photons if p.pt > 20 and abs(p.eta) < 2.4 and not 1.444 < abs(p.eta) < 1.566]
        if Level == "Region" :
            self.CutPlot(Level, 6)  #Ph_kinCuts
            BestPhotons = [gp for gp in GoodPhotons if gp.electronVeto]
            if not BestPhotons : return None
            self.CutPlot(Level, 7)  #Ph_electronVeto
            SelectedPhotons = [bp for bp in BestPhotons if bp.mvaID_WP80]
            if not SelectedPhotons : return None
            self.CutPlot(Level, 8)  #Ph_cutBased>=1
            FinalistPhotons = [sp for sp in SelectedPhotons if all(sp.DeltaR(l) > 0.5 for l in ChLeptons)]
            if not FinalistPhotons : return None
            self.CutPlot(Level, 9)  #Ph_AngularSep
            return FinalistPhotons
        elif Level == "Definition" :
            self.CutPlot(Level, 3)  #Ph_kinCuts
            FinalistPhotons = [gp for gp in GoodPhotons if all(gp.DeltaR(l) > 0.5 for l in ChLeptons)]
            if not FinalistPhotons : return None
            self.CutPlot(Level, 4)  #Ph_AngularSep
            return FinalistPhotons
        else : raise ValueError(f"Not valid option: '{Level}'.")
            

    # return True if the selected photon is NOT a FSR attached to the Z (by looking at the invariant mass)
    def photonIsNotFSR(self, Leptons, photon, massLimit) :
        if len(Leptons) == 2 :
            llGamma_invMass = self.getInvMass(Leptons[0], Leptons[1], photon)
            if llGamma_invMass > massLimit : return True
        return False
        
    #return True if the selected photon is NOT a FSR attached to the W (by looking at the invariant mass)
    def photonIsNotWFSR(self, lepW, photon, massLimit) :
        lGamma_list = [(lepW.pt, lepW.eta, lepW.phi), (photon.pt, photon.eta, photon.phi)]
        lGammaMET_mt = self.computeMt(lGamma_list, self.MET.pt, self.MET.phi)
        if lGammaMET_mt > massLimit : return True
        return False

       
           
   
    #--------------------------
    # UTILITY FUNCTIONS
    #--------------------------
    # Return True if the difference between the gen boson masses and the nominal one is lower for the second pair
    def getLowerRes(self, Wboson1, Zboson1, Wboson2, Zboson2) :
        WMassPdg = TDatabasePDG.Instance().GetParticle(24).Mass()
        ZMassPdg = TDatabasePDG.Instance().GetParticle(23).Mass()
       
        res = lambda b, mass: abs(b.mass - mass) if b is not None else float('inf')
        res1 = res(Wboson1, WMassPdg) + res(Zboson1, ZMassPdg)
        res2 = res(Wboson2, WMassPdg) + res(Zboson2, ZMassPdg)
        if res1 > res2 : return True
        return False
       
   
    # finds the candidate with similar invariant mass to "partMass"
    def getCloserCand(self, partMass, candMass) :
        best = None
        minDiff = 9999.
        for couple, mass in candMass.items() :
            if(mass != None) :  
                invM_diff = abs(partMass - mass)
                if (invM_diff < minDiff) :
                    minDiff = invM_diff
                    best = couple
        return best

   
    # given a variable number of "daughter particles" it calculates the invariant mass
    # it is important that if a photon is passed as attribute it has index = 2 (necessary cause the photon has no pdgId)
    def getInvMass(self, *particles) :
        ptot = ROOT.TLorentzVector()
        for i, part in enumerate(particles):
            if i == 2 and not hasattr(part, "pdgId") :  
                p4_temp = ROOT.TLorentzVector()
                p4_temp.SetPtEtaPhiM(part.pt, part.eta, part.phi, 0.0)
                ptot += p4_temp
            else:
                ptot += part.p4()              
        return ptot.M()


    #return the generation (family) of the lepton
    def getFamily(self, pdgId) :
        a = abs(pdgId)
        if a in (11, 12): return 1  # e
        if a in (13, 14): return 2  # mu
        return 0


    # returns true if the leptons are two from a generation and the third from the other
    def isEventEasy(self, Leptons) :
        if len(Leptons) != 3 : return False
        num_muons = 0
        num_electrons = 0
        for l in Leptons :
            if self.getFamily(l.pdgId) == 1 : num_electrons += 1
            elif self.getFamily(l.pdgId) == 2 : num_muons += 1
        if (num_electrons == 2 and num_muons == 1) or (num_electrons == 1 and num_muons == 2) : return True
        return False

   
    # just for a CHARGE COMPATIBILITY (used in the comparison between WZRecon algorithms)
    def isChargeCompatible(self, pdg1, pdg2, boson):
        CHARGED = {11, 13}
        NEUTRAL = {12, 14}
        a1, a2 = abs(pdg1), abs(pdg2)
        if boson == 24:  
            return (a1 in CHARGED and a2 in NEUTRAL) or (a1 in NEUTRAL and a2 in CHARGED)
        elif boson == 23:
            if a1 not in CHARGED or a2 not in CHARGED: return False
            return pdg1 * pdg2 < 0
        return None

   
    # it says if the two leptons have the correct pdgId characteristic to be "boson daughters"
    def pairIsTruthBoson(self, Leptons, boson_pdgId):
        if len(Leptons) != 2 : return None
        if boson_pdgId == 23 :
            return (Leptons[0].pdgId + Leptons[1].pdgId) == 0      
        if boson_pdgId == 24 :
            isWCouple = (self.getFamily(Leptons[0].pdgId) == self.getFamily(Leptons[1].pdgId)                                                        and self.getFamily(Leptons[0].pdgId) != 0
                        and abs(Leptons[0].pdgId + Leptons[1].pdgId) == 1)
            return isWCouple
        return None


    # computes the transverse mass of any visible sistem of particles + MET
    # where visParticles is a list of tuple (pt, eta, phi)
    def computeMt(self, visParticles, met_pt, met_phi):
        px_tot, py_tot, pz_tot, e_tot = 0.0, 0.0, 0.0, 0.0
        for pt, eta, phi in visParticles:
            px = pt * math.cos(phi)
            py = pt * math.sin(phi)
            pz = pt * math.sinh(eta)
            e = math.hypot(pt, pz)  #for all particles mass considered 0
            px_tot += px
            py_tot += py
            pz_tot += pz
            e_tot += e
    
        pt_sq = px_tot * px_tot + py_tot * py_tot
        p3_sq = pt_sq + pz_tot * pz_tot
        m_sq = max(0.0, e_tot * e_tot - p3_sq)
    
        et_vis = math.sqrt(pt_sq + m_sq)
        met_x = met_pt * math.cos(met_phi)
        met_y = met_pt * math.sin(met_phi)
    
        pt_dot_met = px_tot * met_x + py_tot * met_y  #scalar product pt_vis*met_pt
        mt2 = m_sq + 2.0 * (et_vis * met_pt - pt_dot_met)
        return math.sqrt(max(0.0, mt2))

   
    def areMatched(self, lep_gen, lep_reco, dRLimit) :
        dr = lep_gen.DeltaR(lep_reco)
        if dr < dRLimit and (lep_gen.pdgId - lep_reco.pdgId) == 0 : return True
        return False

    def isRecoMatching(self, lepGenZ, lepGenW, lepRecoZ, lepW_reco, dRLimit) :
        if lepW_reco is None or lepGenW is None or len(lepGenZ) != len(lepRecoZ) != 2 : return None
           
        lepW_gen = next((l for l in lepGenW if abs(l.pdgId) in (11, 13)), None)
        if not self.areMatched(lepW_gen, lepW_reco, dRLimit) : return False

        idx_GenLep = self.couplePart_toGen(lepRecoZ[0], lepGenZ, dRLimit)
        if idx_GenLep != None :
            remaining_lepGenZ = next((p for i, p in enumerate(lepGenZ) if i != idx_GenLep), None)
            if not self.areMatched(lepRecoZ[1], remaining_lepGenZ, dRLimit) :
                if self.areMatched(lepRecoZ[1], lepGenZ[idx_GenLep], dRLimit) :
                   if not self.areMatched(lepRecoZ[0], remaining_lepGenZ, dRLimit) : return False
                else : return False
        else : return False

        return True
       
   
    # given a single particle extract another from a list by looking at the smaller DR
    def couplePart_toGen(self, part, genParticles, dRLimit) :
        idx_GenPart = None
        for i, gen_part in enumerate(genParticles) :
            dr = part.DeltaR(gen_part)
            if dr < dRLimit and gen_part.pdgId == part.pdgId :
                dRLimit = dr
                idx_GenPart = i
        return idx_GenPart

    
    #--------------------------
    # Graphics functions
    #--------------------------
    #function for producing cutFlowStage plot
    def CutPlot(self, Level, Cut):
        if Level == "Definition": Labels = self.SDCuts
        elif Level == "Region": Labels = self.SRCuts
        else : raise ValueError(f"Not valid option: '{Level}'.")
        self.hEvent.fill1D_label(f"EventsPerCuts_Signal{Level}", f"EventsPerCuts_Signal{Level}", Labels, Labels[Cut], self.weight)

    # returns useful plot for computing percentages that characterise the sample : rec. efficiency, background percentage in SR, events entering SR
    def getPercentagePlot(self, option) :
        label = []
        if option == "ReconstructonEfficency" :
            label = ["isSignal", "isSig_and_inSigReg"]
            if self.isSig : 
                self.hEvent.fill1D_label(option, option, label, label[0], self.weight)
                if self.isInSR : self.hEvent.fill1D_label(option, option, label, label[1], self.weight)
        elif option == "Background_inSR" :
            label = ["totalN_events","is_inSigReg", "is_background"]
            self.hEvent.fill1D_label(option, option, label, label[0], self.weight)
            if self.isInSR : 
                self.hEvent.fill1D_label(option, option, label, label[1], self.weight)
                if not self.isSig : self.hEvent.fill1D_label(option, option, label, label[2], self.weight)
        else : raise ValueError(f"Not valid option: '{option}'.")

    
    '''# useful for plotting some common hist between Gen and Reco Level
    def commonPlots(self, Level, , Gamma):
        self.hEvent.fill1D(f"{Level}Level_Z1Mass", f"{Level}Level_Z1Mass", 100, 60., 120., ZZ.Z1mass, self.weight)
        self.hEvent.fill1D(f"{Level}Level_Z2Mass", f"{Level}Level_Z2Mass", 75, 60., 120., ZZ.Z2mass, self.weight)
        self.hEvent.fill1D(f"{Level}Level_ZZMass", f"{Level}Level_ZZMass", 40, 120., 520., ZZ.mass, self.weight)

        if Level == "Gen":
            ZZGammaInvMass = (ZZ.p4 + Gamma.p4()).M()
        elif Level == "Reco":
            ZZGammaInvMass = (ZZ.p4() + self.GetP4(Gamma)).M()
        self.hEvent.fill1D(f"{Level}Level_ZZGammaMass", f"{Level}Level_ZZGammaMass", 50, 120., 620., ZZGammaInvMass, self.weight)

        self.hEvent.fill1D(f"{Level}Level_GammaPt", f"{Level}Level_GammaPt", 18, 0., 180., Gamma.pt, self.weight)
        self.hEvent.fill1D(f"{Level}Level_GammaEta", f"{Level}Level_GammaEta", 12, 0., 2.4, Gamma.eta, self.weight)'''
               
   
#---------------------------
# other CLASSES
#---------------------------
class GenZ:  
    def __init__(self, Leptons) :
        self.leptons = Leptons
        self.p4 = Leptons[0].p4() + Leptons[1].p4()
        self.mass = self.p4.M()
        self.pt = self.p4.Pt()
        self.eta = self.p4.Eta()
        self.phi = self.p4.Phi()  

class GenW:
    def __init__(self, Leptons) :
        self.leptons = Leptons
        self.p4 = Leptons[0].p4() + Leptons[1].p4()
        self.mass = self.p4.M()
        self.pt = self.p4.Pt()
        self.eta = self.p4.Eta()
        self.phi = self.p4.Phi()

class WReco :
    def __init__(self, Lepton, transvers_mass, Fsr = None, LeptonDressed = None) :
        self.lepton = Lepton
        self.transvMass = transvers_mass
        self.fsrPhoton = Fsr
        self.leptonDressed = LeptonDressed


class DressedLepton :
    def __init__(self, nakedLepton, Photon = None ) :
        self.nakedLepton = nakedLepton
        self.pdgId = nakedLepton.pdgId
        self.charge = nakedLepton.charge
        self.fsrPhoton = Photon

        p4_total = nakedLepton.p4()
        if Photon is not None:
            p4_total += Photon.p4()

        self._p4 = p4_total
        #kinematic
        self.pt = self._p4.Pt()     
        self.eta = self._p4.Eta()
        self.phi = self._p4.Phi()
        self.mass = self._p4.M()

    def p4(self):
        return self._p4


#useful for a first check on gen particles
class GenStatusFlag(IntFlag) :
    IS_PROMPT = 1 << 0
    FROM_HARD_PROCESS = 1 << 8

    SEL = IS_PROMPT | FROM_HARD_PROCESS        
    #PH_SEL = IS_PROMPT   #removed from_hard_Process for taking also the FSR photon

    