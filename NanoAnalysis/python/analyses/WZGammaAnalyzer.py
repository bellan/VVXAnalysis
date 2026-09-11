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

    def analyze(self):   
        cutFlowStage = 0.
        cutFlowStageSR = 0.
        if self.event.HLT_passZZ4l :
            weight = 1.   #da implementare.....  

                
            GenParts = Collection(self.event, 'GenPart')
            GenChargedLeptons = [p for p in GenParts if abs(p.pdgId) in (11,13) and p.status == 1
                                 and (p.statusFlags & GenStatusFlag.SEL) == GenStatusFlag.SEL]
            GenPhotons = [p for p in GenParts if p.pdgId == 22 and p.status == 1 
                          and (p.statusFlags & GenStatusFlag.SEL) == GenStatusFlag.SEL]
            GenNeutrinos = [p for p in GenParts if abs(p.pdgId) in (12,14)  and p.status == 1 
                            and (p.statusFlags & GenStatusFlag.SEL) == GenStatusFlag.SEL] 
                            
            
            Electrons = Collection(self.event, 'Electron')
            Muons = Collection(self.event, 'Muon')
            ChargedLeptons = list(Electrons) + list(Muons)
            Photons = Collection(self.event, 'Photon')
            ZCands = Collection(self.event, 'ZCand')           
            bestZIdx = self.event.bestZIdx 
            PFMET = Object(self.event, "PFMET")    
            PuppiMET = Object(self.event, "PuppiMET")

            
            #---------------------      
            #SIGNAL DEFINITION
            #---------------------
            def isSignal() :                
                # KINEMATIC CUT ON GENERATED PARTICLES
                self.hEvent.fill1D("nGenChargedLeptons", "nGenChargedLeptons",11, -0.5, 10.5, len(GenChargedLeptons), self.weight)
                if self.pass_ChLepKinCut(GenChargedLeptons) is None : return False, 1.
                else : GenBestLeptons = self.pass_ChLepKinCut(GenChargedLeptons)  
                self.hEvent.fill1D("nGenChargedLeptons_PostCut", "nGenChargedLeptons_PostCut",11, -0.5, 10.5, len(GenBestLeptons), self.weight)

                self.hEvent.fill1D("nGenPhoton", "nGenPhoton",11, -0.5, 10.5, len(GenPhotons), self.weight)
                genPhoton = self.pass_PhKinCut(GenPhotons, GenChargedLeptons)
                if genPhoton is None : return False, 2.


                # COMPARISON OF WZ RECONSTRUCTION ALGORITHMS
                if self.isEventEasy(GenBestLeptons) :
                    error_ZRec = 0
                    error_WRec = 0
                    error_WZRec = 0
                    
                    bCouple_ZRec, fail_ZRec = self.WZRecon(GenBestLeptons, GenNeutrinos, False , False, True)
                    theW_ZRec, theZ_ZRec = bCouple_ZRec   
                    if fail_ZRec == 1 : outcome_ZRec = 2   # fail
                    elif not self.pairIsTruthBoson(theZ_ZRec.leptons, 23) or not self.pairIsTruthBoson(theW_ZRec.leptons,                        24): outcome_ZRec = 1   # error
                    else: outcome_ZRec = 0   # success
                    self.hEvent.fill1D("MethodZ_Outcome", "MethodZ_Outcome", 3, -0.5, 2.5, outcome_ZRec, self.weight)
                    
                    bCouple_WRec, fail_WRec = self.WZRecon(GenBestLeptons, GenNeutrinos, False , True, True)
                    theW_WRec, theZ_WRec = bCouple_WRec
                    if fail_WRec == 1 : outcome_WRec = 2   # fail
                    elif not self.pairIsTruthBoson(theZ_WRec.leptons, 23) or not self.pairIsTruthBoson(theW_WRec.leptons,                        24): outcome_WRec = 1   # error
                    else: outcome_WRec = 0   # success
                    self.hEvent.fill1D("MethodW_Outcome", "MethodW_Outcome", 3, -0.5, 2.5, outcome_WRec, self.weight)

                    bCouple_WZRec, fail_WZRec = self.WZRecon(GenBestLeptons, GenNeutrinos, True, False, True)
                    theW_WZRec, theZ_WZRec = bCouple_WZRec
                    if fail_WZRec == 1 : outcome_WZRec = 2   # fail
                    elif not self.pairIsTruthBoson(theZ_WZRec.leptons, 23) or not self.pairIsTruthBoson(theW_WZRec.leptons,                        24): outcome_WZRec = 1   # error
                    else: outcome_WZRec = 0   # success
                    self.hEvent.fill1D("MethodWZ_Outcome", "MethodWZ_Outcome", 3, -0.5, 2.5, outcome_WZRec, self.weight)

                
                # W and Z RECONSTRUCTION (with the winner method)
                bCouple, fail = self.WZRecon(GenBestLeptons, GenNeutrinos, True , False, False)
                theW, theZ = bCouple
                if theW is None or theZ is None : return False, 3.
                else : 
                    self.hEvent.fill1D("ZGenMass_91Gev", "ZGenMass_91Gev", 180, 60., 120., theZ.mass, self.weight)
                    self.hEvent.fill1D("WGenMass_80Gev", "WGenMass_80Gev", 160, 50., 110., theW.mass, self.weight)
                  
                            
                notpassedZ, invMassLLGamma = self.isPhotonFSR(theZ.leptons, genPhoton, 100)
                #self.hEvent.fill1D("llGammaGen_invMass", "llGammaGen_invMass", 100, 40., 250., invMassLLGamma, self.weight)
                # only cause the neutrino is gen
                notpassedW, invMassLNuGamma = self.isPhotonFSR(theW.leptons, genPhoton, 90)
                #self.hEvent.fill1D("lNuGammaGen_invMass", "lNuGammaGen_invMass", 100, 40., 250., invMassLNuGamma, self.weight)
                self.hEvent.fill1D("MinLLGammaGen_invMass", "MinLLGammaGen_invMass", 100, 40., 250., min(invMassLLGamma, invMassLNuGamma), self.weight)
                if notpassedW or notpassedZ : 
                    #printout evento
                    '''print(f"\n>>> DEBUG EVENT: Event={self.event}")
                    for i, p in enumerate(Collection(self.event, "GenPart")):
                        if (abs(p.pdgId) in (11, 12, 13, 14, 15, 16, 22, 23, 24)) :
                            print(f"Idx: {i:2d} | PDG: {p.pdgId:4d} | Status: {p.status:2d} | "
                                  f"pT: {p.pt:6.1f} | eta: {p.eta:5.2f} | MotherIdx: {p.genPartIdxMother:2d} | "
                                  f"Flags: {p.statusFlags:015b}")'''
                    return False, 4.

                #control on neutrinos number
                self.hEvent.fill1D("nGenNeutrino_postCuts", "nGenNeutrino_postCuts",7, -1.5, 5.5, len(GenNeutrinos), self.weight)

                return True, -1.
                                
                
            
            #-----------------------
            # SIGNAL REGION
            #-----------------------
            def isInSignalRegion() :
                # KINEMATIC CUT ON PARTICLES
                if self.pass_ChLepKinCut(ChargedLeptons) is None : return False, 1.
                else : BestLeptons = self.pass_ChLepKinCut(ChargedLeptons)  
                                           
                bestPhoton = self.pass_PhKinCut(Photons, ChargedLeptons)
                if bestPhoton is None : return False, 2.

                # CUT on MISSING ENERGY
                self.hEvent.fill1D("PFMET_pt", "PFMET_pt", 121, -0.5, 120.5, PFMET.pt, self.weight)
                self.hEvent.fill1D("PuppiMET_pt", "PuppiMET_pt", 121, -0.5, 120.5, PuppiMET.pt, self.weight)
                if PFMET.pt <= 30 or PuppiMET.pt <= 30 : return False, 3.

                # Z RECONSTRUCTION AND COMPARE
                ZRec, recCouple = self.BosonRecon(BestLeptons, 23, (60,120), False)
                if recCouple != None :
                    self.hEvent.fill1D("ZRecMass_91Gev_SR", "ZRecMass_91Gev_SR", 180, 60., 120., ZRec.mass, self.weight)
                
                theZ = ZCands[bestZIdx]
                if 60. < theZ.mass < 120. : 
                    self.hEvent.fill1D("ZCandMass_91Gev_SR", "ZCandMass_91Gev_SR", 180, 60., 120., theZ.mass, self.weight)
                    lep1, lep2 = ChargedLeptons[theZ.l1Idx], ChargedLeptons[theZ.l2Idx]
                    lepW = next((l for l in BestLeptons if l not in (lep1, lep2)), None)
                else : return False, 4.
                    

                notpassedZ, invMassLLGamma = self.isPhotonFSR([lep1,lep2], bestPhoton, 100)
                self.hEvent.fill1D("llGamma_invMass", "llGamma_invMass", 100, 40., 250., invMassLLGamma, self.weight)
                if notpassedZ : return False, 5.


                W_mt_pf = self.computeMt(lepW.pt, lepW.phi, PFMET.pt, PFMET.phi)
                self.hEvent.fill1D("W_Mt_PF", "W_Mt_PF", 71, -1.5, 140.5, W_mt_pf, self.weight)
                W_mt_puppi = self.computeMt(lepW.pt, lepW.phi, PuppiMET.pt, PuppiMET.phi)
                self.hEvent.fill1D("W_Mt_Puppi", "W_Mt_Puppi", 71, -1.5, 140.5, W_mt_puppi, self.weight)

                return True, -1.

            

            if self.analyzeMC : isSig, cutFlowStage = isSignal()   
            
            isInSR, cutFlowStageSR = isInSignalRegion()

        self.hEvent.fill1D("CutFlowStage", "CutFlowStage",8, -1.5, 6.5, cutFlowStage, self.weight)
        self.hEvent.fill1D("CutFlowStage_SignalRegion", "CutFlowStage_SignalRegion",8, -1.5, 6.5, cutFlowStageSR, self.weight)
        
    
          
        
    #----------------------------------
    # Functions for KINEMATIC CUT 
    #----------------------------------
    # if the cut is passed (and also the trigger 20-10-5) it returns the three charged leptons with the biggest pt
    def pass_ChLepKinCut(self, ChLeptons) :
        GoodChLeptons = [p for p in ChLeptons if p.pt > 5 and abs(p.eta) < 2.5]
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
        if(bestCouple != None and massThreshold[0] < CandMass[bestCouple] < massThreshold[1]) :  
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
        minDiff = 999. 
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

    
    # computes the transverse mass of the W boson
    def computeMt(self, lep_pt, lep_phi, met_pt, met_phi):
        return math.sqrt(2.0 * lep_pt * met_pt * (1.0 - math.cos(lep_phi - met_phi)))
        

    
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


#useful for a first check on gen particles
class GenStatusFlag(IntFlag) :
    IS_PROMPT = 1 << 0
    FROM_HARD_PROCESS = 1 << 8 

    SEL = IS_PROMPT | FROM_HARD_PROCESS         
    #PH_SEL = IS_PROMPT   #removed from_hard_Process for taking also the FSR photon
    