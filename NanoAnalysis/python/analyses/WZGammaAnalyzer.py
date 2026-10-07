from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection
from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Object

from VVXAnalysis.NanoAnalysis.EventAnalyzer import EventAnalyzer
from VVXAnalysis.NanoAnalysis.Histogrammer import *
from VVXAnalysis.NanoAnalysis.Regions import Flags as Regions

import ROOT
import math
from ROOT import TDatabasePDG  
from enum import IntFlag


class WZGammaAnalyzer(EventAnalyzer, analysis_name="WZGammaAnalyzer"):
     
    def __init__(self, regions):
        super().__init__(regions)

   
    # creates particle collections
    def loadCollections(self) :
        if self.analyzeMC:
            self.GenParts = Collection(self.event, 'GenPart')
            self.GenChargedLeptons = [p for p in self.GenParts if abs(p.pdgId) in (11,13) and p.status == 1
                                and (p.statusFlags & GenStatusFlag.SEL) == GenStatusFlag.SEL]
            self.GenPhotons = [p for p in self.GenParts if p.pdgId == 22 and p.status == 1
                              and (p.statusFlags & GenStatusFlag.SEL) == GenStatusFlag.SEL]
            self.GenNeutrinos = [p for p in self.GenParts if abs(p.pdgId) in (12,14)  and p.status == 1
                            and (p.statusFlags & GenStatusFlag.SEL) == GenStatusFlag.SEL]
                               

        self.event.SetBranchStatus("*ZCand*", 1)
        self.event.SetBranchStatus("*PFMET*", 1)
        self.event.SetBranchStatus("*PuppiMET*", 1)
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

        self.cutFlowStage = 0.0
        self.cutFlowStageSR = 0.0  
        isSig = False
        if self.analyzeMC :
            isSig = self.isSignal()  
            genZ, genW = self.genCouple
               
        isInSR = self.isInSignalRegion(2024)
        recoZ, recoW = self.recoCouple

        if self.hasGenFSR and self.lGammaMET_mt != None :
            self.hEvent.fill1D("lGammaMET_Mt_Cutted", "lGammaMET_Mt_Cutted", 30, 0., 300., self.lGammaMET_mt, self.weight)

        if isSig :
            eventType_Labels = ["isSignal", "isSig_and_inSigReg"]
            event_type = 0
           
            if isSig and isInSR :
                event_type = 1
               
                # GEN-RECO LEPTONS MATCHING
                res = self.isRecoMatching(genZ.leptons, genW.leptons, (self.ChargedLeptons[recoZ.l1Idx], self.ChargedLeptons[recoZ.l2Idx]), recoW.lepton, 0.5)
                self.hEvent.fill1D("GenAndReco_Matching", "GenAndReco_Matching", 2, -0.5, 1.5, int(res), self.weight)
        

            self.hEvent.fill1D_label("Event_type", "Event_type", eventType_Labels, eventType_Labels[event_type], self.weight)

        cutLabels_sigDef = ["Events_passed", "Lep_Cuts", "Ph_Cuts", "WZ_MassLimits", "FSR_onWZ"]
        cutLabels_sigReg = ["Events_passed", "Lep_Cuts", "Ph_Cuts", "MET_Cut","Z_MassLimits", "FSR_onZ"]

        self.hEvent.fill1D_label("CutFlowStage_Gen", "CutFlowStage_Gen", cutLabels_sigDef, cutLabels_sigDef[int(round(self.cutFlowStage))], self.weight)
        self.hEvent.fill1D_label("CutFlowStage_Reco", "CutFlowStage_Reco", cutLabels_sigReg, cutLabels_sigReg[int(round(self.cutFlowStageSR))], self.weight)

        entersSigReg = 0
        if isInSR and self.cutFlowStageSR == 0.0 : entersSigReg = 1
        sigRegTag_Labels = ["isNOT_inSigReg", "is_inSigReg"]
        self.hEvent.fill1D_label("Signal_Region", "Signal_Region", sigRegTag_Labels, sigRegTag_Labels[entersSigReg], self.weight)
        self.hEvent.fill1D("Signal_Region_calc", "Signal_Region_calc", 2, -0.5, 1.5, entersSigReg, self.weight)
        
    
       

    #---------------------      
    #SIGNAL DEFINITION
    #---------------------
    def isSignal(self) :        
        self.cutFlowStage = 0.0
        self.genCouple = None, None
        self.W_mt_Gen = None
        self.lGammaMET_mt_Gen = None
        self.hasGenFSR = False
        self.GenBestLepton = None
           
        # KINEMATIC CUT ON GENERATED PARTICLES      
       
        self.hEvent.fill1D("nGenChargedLeptons", "nGenChargedLeptons",11, -0.5, 10.5, len(self.GenChargedLeptons), self.weight)
        self.GenBestLeptons = self.pass_ChLepKinCut(self.GenChargedLeptons)  
        if self.GenBestLeptons is None :
            self.cutFlowStage = 1.0
            return False
        self.hEvent.fill1D("nGenChargedLeptons_postCuts", "nGenChargedLeptons_postCuts",11, -0.5, 10.5, len(self.GenBestLeptons), self.weight)

        self.hEvent.fill1D("nGenPhoton", "nGenPhoton",11, -0.5, 10.5, len(self.GenPhotons), self.weight)
        genPhoton = self.pass_PhKinCut(self.GenPhotons, self.GenChargedLeptons)
        if genPhoton is None :
            self.cutFlowStage = 2.0
            return False
        self.hEvent.fill1D("Photon_pt_Gen", "Photon_pt_Gen", 50, 20., 200., genPhoton.pt, self.weight)


        # COMPARISON OF WZ RECONSTRUCTION ALGORITHMS           
        if self.isEventEasy(self.GenBestLeptons) :        
            self.methodOutcome(False, False,"MethodZ_Outcome")  #Z-first
            self.methodOutcome(False, True, "MethodW_Outcome")   #W-first
            self.methodOutcome(True, False, "MethodWZ_Outcome")   #residues

       
        # W and Z RECONSTRUCTION (with the Z-first method)
        bCouple, fail = self.WZRecon(self.GenBestLeptons, self.GenNeutrinos, False , False, False)
        theW, theZ = bCouple
        if theW is None or theZ is None :
            self.cutFlowStage = 3.0
            return False
           
        self.hEvent.fill1D("ZGenMass_91Gev", "ZGenMass_91Gev", 180, 60., 120., theZ.mass, self.weight)
        self.hEvent.fill1D("WGenMass_80Gev", "WGenMass_80Gev", 160, 50., 110., theW.mass, self.weight)
        self.genCouple = (theZ, theW)
        lepW = next((l for l in theW.leptons if abs(l.pdgId) in (11,13)), None)
        nu  = next((n for n in theW.leptons if abs(n.pdgId) in  (12,14)), None)

        #computation of different transverse mass
        if (nu.pt >= 30.0 ) :  #same condition as at reco-level
            self.W_mt_Gen = self.computeMt(lepW.pt, lepW.phi, nu.pt, nu.phi)
            self.lGammaMET_mt_Gen = self.computeMt(lepW.pt, lepW.phi, nu.pt, nu.phi, lepW.eta, genPhoton.pt, genPhoton.eta, genPhoton.phi)
            self.hEvent.fill1D("Gen_W_Mt", "Gen_W_Mt", 30, 0., 300., self.W_mt_Gen, self.weight)
            self.hEvent.fill1D("lGammaMET_Mt_Gen", "lGammaMET_Mt_Gen", 30, 0., 300., self.lGammaMET_mt_Gen, self.weight)
            #2d plot
            self.hEvent.fill2D("W_vs_lGamma_Mt_Gen", "W_vs_lGamma_Mt_Gen", 30, 0., 300., 12, 0., 120., self.lGammaMET_mt_Gen, self.W_mt_Gen, self.weight)
       
        # invariant mass "llGamma" analysis          
        notpassedZ, invMassLLGamma = self.isPhotonFSR(theZ.leptons, genPhoton, 100)
        notpassedW, invMassLNuGamma = self.isPhotonFSR(theW.leptons, genPhoton, 90)    # only cause the neutrino is gen
        self.hEvent.fill1D("LNuGammaGen_invMass", "LNuGammaGen_invMass", 50, 40., 240., invMassLNuGamma, self.weight)               # 2d plot
        invMassLL = theZ.mass  
        invMassLNu = theW.mass  
        self.hEvent.fill2D("llG_vs_ll_invMass", "llG_vs_ll_invMass", 100, 40., 240., 30, 60., 120., invMassLLGamma, invMassLL, self.weight)
        self.hEvent.fill2D("lNuG_vs_lNu_invMass", "lNuG_vs_lNu_invMass", 100, 40., 240., 30, 50., 110., invMassLNuGamma, invMassLNu, self.weight)
        if notpassedW and nu.pt >= 30.0 :
            self.hEvent.fill1D("lGammaMET_Mt_CuttedGen", "lGammaMET_Mt_CuttedGen", 30, 0., 300., self.lGammaMET_mt_Gen, self.weight)
            self.hasGenFSR = True
        if notpassedW or notpassedZ :
            self.cutFlowStage = 4.0
            return False

        #control on neutrinos number
        self.hEvent.fill1D("nGenNeutrino_postCuts", "nGenNeutrino_postCuts",7, -1.5, 5.5, len(self.GenNeutrinos), self.weight)

        return True

                                 
   
    #-----------------------
    # SIGNAL REGION
    #-----------------------
    def isInSignalRegion(self, year) :
        self.cutFlowStageSR = 0.0
        self.recoCouple = (None, None)
        self.W_mt = None
        self.lGammaMET_mt = None
        self.MET = None
       
        # KINEMATIC CUT ON PARTICLES
        BestLeptons = self.pass_ChLepKinCut(self.ChargedLeptons)
        if BestLeptons is None :
            self.cutFlowStageSR = 1.0
            return False
   
                                   
        bestPhoton = self.pass_PhKinCut(self.Photons, self.ChargedLeptons)
        if bestPhoton is None :
            self.cutFlowStageSR = 2.0
            return False

        self.hEvent.fill1D("Photon_pt_beforeCuts", "Photon_pt_beforeCuts", 30, 0., 300., bestPhoton.pt, self.weight)

        # CUT on MISSING ENERGY
        puppyYears = [2024]
        self.MET = self.PuppiMET if year in puppyYears else self.PFMET
       
        self.hEvent.fill1D("MET_pt_beforeCuts", "MET_pt_beforeCuts", 30, 0., 300., self.MET.pt, self.weight)
        if self.MET.pt <= 30 :
            self.cutFlowStageSR = 3.0
            return False

        # Z RECONSTRUCTION
        theZ = self.ZCands[self.bestZIdx]
        if 60. < theZ.mass < 120. :
            self.hEvent.fill1D("ZCandMass_91Gev_SR", "ZCandMass_91Gev_SR", 180, 60., 120., theZ.mass, self.weight)
            lep1, lep2 = self.ChargedLeptons[theZ.l1Idx], self.ChargedLeptons[theZ.l2Idx]
            lepW = next((l for l in BestLeptons if l not in (lep1, lep2)), None)
        else :
            self.cutFlowStageSR = 4.0
            return False
           
        notpassedZ, invMassLLGamma = self.isPhotonFSR([lep1,lep2], bestPhoton, 100)
        self.hEvent.fill1D("llGamma_invMass", "llGamma_invMass", 40, 40., 240., invMassLLGamma, self.weight)
        if notpassedZ :
            self.cutFlowStageSR = 5.0
            return False
           
        # W RECONSTRUCTION (looking at different kind of FSR now (ones inside 0.5 of DR))
        self.W_mt = self.computeMt(lepW.pt, lepW.phi, self.MET.pt, self.MET.phi)
        self.hEvent.fill1D("W_Mt_NotDressed", "W_Mt_NotDressed", 30, 0., 300., self.W_mt, self.weight)

        isPaired, fsrPhoton = self.isFsrPaired(lepW, theZ)
        if  isPaired :
            self.W_mt = self.computeMt(lepW.pt, lepW.phi, self.MET.pt, self.MET.phi, lepW.eta, fsrPhoton.pt, fsrPhoton.eta, fsrPhoton.phi)
            theW = WReco(lepW, self.W_mt, fsrPhoton)
            self.hEvent.fill1D("W_Mt_OnlyDressed", "W_Mt_OnlyDressed", 30, 0., 300., self.W_mt, self.weight)
        else : theW = WReco(lepW, self.W_mt)
        self.recoCouple = (theZ, theW)

        # Transverse mass lGamma con MET
        self.lGammaMET_mt = self.computeMt(lepW.pt, lepW.phi, self.MET.pt, self.MET.phi, lepW.eta, bestPhoton.pt, bestPhoton.eta, bestPhoton.phi)
        self.hEvent.fill1D("W_Mt", "W_Mt", 30, 0., 300., self.W_mt, self.weight)
        self.hEvent.fill1D("lGammaMET_Mt", "lGammaMET_Mt", 30, 0., 300., self.lGammaMET_mt, self.weight)
        self.hEvent.fill2D("W_vs_lGamma_Mt", "W_vs_lGamma_Mt", 30, 0., 300., 12, 0., 120., self.lGammaMET_mt, self.W_mt, self.weight)

        #plot post cuts
        self.hEvent.fill1D("Photon_pt_postCuts", "Photon_pt_postCuts", 30, 0., 300., bestPhoton.pt, self.weight)
        self.hEvent.fill1D("MET_pt_postCuts", "MET_pt_postCuts", 30, 0., 300., self.MET.pt, self.weight)

        return True

       
         
       
    #----------------------------------
    # Functions for KINEMATIC CUT
    #----------------------------------
    # if the cut is passed (and also the trigger 20-10-5) it returns the three charged leptons with the biggest pt
    def pass_ChLepKinCut(self, ChLeptons) :
        GoodChLeptons = [p for p in ChLeptons if p.pt > 5 and abs(p.eta) < 2.5]
        ret = {}
        if len(GoodChLeptons) > 2 :
            l1 = max(GoodChLeptons, key=lambda p: p.pt, default=None)
            if l1.pt > 20 :
                l2 = max((p for p in GoodChLeptons if p != l1), key=lambda p: p.pt, default=None)
                if l2.pt  > 10 :
                    l3 = max((p for p in GoodChLeptons if p not in (l1, l2)), key=lambda p: p.pt, default=None)
                    return [l1,l2,l3]                  
        return None

   
    # for kinematic cut on photons: returns the photons that passed the cuts with the biggest pt
    def pass_PhKinCut(self, Photons, ChLeptons) :
        GoodPhotons = [p for p in Photons if p.pt > 20 and abs(p.eta) < 2.4 and not 1.444 < abs(p.eta) < 1.566]    
        BestPhotons = [gp for gp in GoodPhotons if all(gp.DeltaR(l) > 0.5 for l in ChLeptons)]
        if len(BestPhotons) > 0 :
            return max(BestPhotons, key=lambda p: p.pt, default=None)
        else : return None

   
 
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
   
        fail = int(bosonCouple[0] is None or bosonCouple[1] is None)        
        return bosonCouple, fail


   
    # return True if the selected photon is a FSR by looking at the invariant mass
    def isPhotonFSR(self, Leptons, photon, massLimit) :
        res = None
        if len(Leptons) == 2 :
            res = True
            llGamma_mass = self.getInvMass(Leptons[0], Leptons[1], photon)
            if llGamma_mass > massLimit : res = False
        return res, llGamma_mass


    # Method that studies the algoritm results in terms of "error" (looks at the lepton coupling)
    def methodOutcome(self, useResidues, recbyW, histName) :
        bCouple, fail = self.WZRecon(self.GenBestLeptons, self.GenNeutrinos, useResidues, recbyW, True)
        theW, theZ = bCouple
        if theW is not None and theZ is not None : 
            if not self.pairIsTruthBoson(theZ.leptons, 23) or not self.pairIsTruthBoson(theW.leptons,                        24): outcome = 0   # error
            else: outcome = 1   # success
            self.hEvent.fill1D(histName, histName, 2, -0.5, 1.5, outcome, self.weight)
       
           
   
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

   
    # it says if the two leptons have the correct pdgId characteristich to be "boson daughters"
    def pairIsTruthBoson(self, Leptons, boson_pdgId):
        if len(Leptons) != 2 : return None
        if boson_pdgId == 23 :
            return (Leptons[0].pdgId + Leptons[1].pdgId) == 0      
        if boson_pdgId == 24 :
            isWCouple = (self.getFamily(Leptons[0].pdgId) == self.getFamily(Leptons[1].pdgId)                                                        and self.getFamily(Leptons[0].pdgId) != 0
                        and abs(Leptons[0].pdgId + Leptons[1].pdgId) == 1)
            return isWCouple
        return None

   
    # computes the transverse mass of the W boson and also in the case WGamma
    def computeMt(self, lep_pt, lep_phi, met_pt, met_phi, lep_eta=None, pho_pt=None, pho_eta=None, pho_phi=None):
        if pho_pt is None:
            return math.sqrt(2.0 * lep_pt * met_pt * (1.0 - math.cos(lep_phi - met_phi)))
   
        lx, ly, lz = lep_pt * math.cos(lep_phi), lep_pt * math.sin(lep_phi), lep_pt * math.sinh(lep_eta)
        gx, gy, gz = pho_pt * math.cos(pho_phi), pho_pt * math.sin(pho_phi), pho_pt * math.sinh(pho_eta)
   
        px, py = lx + gx, ly + gy
        pt2 = px * px + py * py
        m2 = max(0.0, (math.hypot(lx, ly, lz) + math.hypot(gx, gy, gz)) ** 2 - pt2 - (lz + gz) ** 2)   #E_tot^2 - p_t^2 - p_z^2
        Et = math.sqrt(pt2 + m2)
   
        return math.sqrt(max(0.0, m2 + 2.0 * (Et * met_pt - (px * met_pt * math.cos(met_phi) + py * met_pt * math.sin(met_phi)))))

   
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


    # returns True and the FSRPhoton if the FSRPhoton selected as candidate matches the index with the W lepton
    def isFsrPaired(self, lepW, theZ) :
        FsrPhotonCands = [p for i,p in enumerate(self.FSRPhotons) if i not in (theZ.fsr1Idx, theZ.fsr2Idx)]
        fsrPhotonCand = min(FsrPhotonCands, key=lambda p: p.dROverEt2, default=None)

        idx = lepW._index
        if fsrPhotonCand is None or (fsrPhotonCand.muonIdx == fsrPhotonCand.electronIdx == -1): return False, None      
        elif (fsrPhotonCand.electronIdx != -1 and abs(lepW.pdgId) == 11) :  
            if fsrPhotonCand.electronIdx == idx : return True, fsrPhotonCand
        elif (fsrPhotonCand.muonIdx != -1 and abs(lepW.pdgId) == 13) : 
            if fsrPhotonCand.muonIdx == idx : return True, fsrPhotonCand
        return False, None  

   
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
    def __init__(self, Lepton, transvers_mass, Fsr = None) :
        self.lepton = Lepton
        self.transvMass = transvers_mass
        self.fsrPhoton = Fsr



#useful for a first check on gen particles
class GenStatusFlag(IntFlag) :
    IS_PROMPT = 1 << 0
    FROM_HARD_PROCESS = 1 << 8

    SEL = IS_PROMPT | FROM_HARD_PROCESS        
    #PH_SEL = IS_PROMPT   #removed from_hard_Process for taking also the FSR photon

    
