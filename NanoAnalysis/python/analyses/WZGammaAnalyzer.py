from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection
from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Object

from VVXAnalysis.NanoAnalysis.EventAnalyzer import EventAnalyzer
from VVXAnalysis.NanoAnalysis.Histogrammer import *
from VVXAnalysis.NanoAnalysis.Regions import Flags as Regions

import ROOT
import math
from ROOT import TDatabasePDG   


class WZGammaAnalyzer(EventAnalyzer, analysis_name="WZGammaAnalyzer"):
     
    def __init__(self, regions):
        super().__init__(regions)

    def analyze(self):   
        cutFlowStage = 0.
        cutFlowStageSR = 0.
        if self.event.HLT_passZZ4l :
            weight = 1.   #da implementare.....  
    
            GenParts = Collection(self.event, 'GenPart')
            GenChargedLeptons = [p for p in GenParts if abs(p.pdgId) in (11,13)]
            GenPhotons = [p for p in GenParts if p.pdgId == 22]
            GenNeutrinos = [p for p in GenParts if abs(p.pdgId) in (12,14)]
            
            Electrons = Collection(self.event, 'Electron')
            Muons = Collection(self.event, 'Muon')
            #ChargedLeptons = Collection(self.event, 'Lepton')
            ChargedLeptons = list(Electrons) + list(Muons)
            Photons = Collection(self.event, 'Photon')
            ZCands = Collection(self.event, 'ZCand')
            MET_pt = self.event.PFMET_pt 
            bestZIdx = self.event.bestZIdx 

            
            #---------------------      
            #SIGNAL DEFINITION
            #---------------------
            def isSignal() : 
                # KINEMATIC CUT ON GENERATED PARTICLES
                if self.pass_ChLepKinCut(GenChargedLeptons) is None : return False, 1.
                else : GenBestLeptons = self.pass_ChLepKinCut(GenChargedLeptons)  
                                           
                genPhoton = self.pass_PhKinCut(GenPhotons, GenChargedLeptons)
                if genPhoton is None : return False, 2.
               
                # W and Z RECONSTRUCTION
                theW, theZ = self.WZRecon(GenBestLeptons, GenNeutrinos, True, False)
                if theW is None or theZ is None : return False, 3.
                else : 
                    self.hEvent.fill1D("ZGenMass_91Gev", "ZGenMass_91Gev", 180, 60., 120., theZ.mass, self.weight)
                    self.hEvent.fill1D("WGenMass_80Gev", "WGenMass_80Gev", 160, 50., 110., theW.mass, self.weight)
                    #self.hEvent.fill1D("WPt", "WPt", 90, 0., 180., theW.pt, self.weight)
       
                notpassedZ, invMassLLGamma = self.isPhotonFSR(theZ.leptons, genPhoton, 110)
                #self.hEvent.fill1D("llGammaGen_invMass", "llGammaGen_invMass", 100, 40., 250., invMassLLGamma, self.weight)
                if notpassedZ : return False, 4.

                # only cause the neutrino is gen
                notpassedW, invMassLNuGamma = self.isPhotonFSR(theW.leptons, genPhoton, 100)
                #self.hEvent.fill1D("lNuGammaGen_invMass", "lNuGammaGen_invMass", 100, 40., 250., invMassLNuGamma, self.weight)
                if notpassedW : return False, 5.
                  
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
                else : self.hEvent.fill1D("Gamma_pt", "Gamma_pt", 122, -0.0, 121., bestPhoton.pt, self.weight)

                # CUT on MISSING ENERGY
                self.hEvent.fill1D("MET_pt", "MET_pt", 121, -0.5, 120.5, MET_pt, self.weight)
                if MET_pt <= 30: return False, 3.

                # Z RECONSTRUCTION AND COMPARISON (ricostruire anche massa trasversa W)
                '''theZRec, coupleZRec = self.BosonRecon(ChargedLeptons, 23, (60,120), False)
                if theZRec is None : return False, 3.
                else : self.hEvent.fill1D("ZRecMass_91Gev_SR", "ZRecMass_91Gev_SR", 180, 60., 120., theZRec.mass, self.weight)'''

                theZ = ZCands[bestZIdx]
                if 60. < theZ.mass < 120. : 
                    self.hEvent.fill1D("ZCandMass_91Gev_SR", "ZCandMass_91Gev_SR", 180, 60., 120., theZ.mass, self.weight)
                    lep1, lep2 = ChargedLeptons[theZ.l1Idx], ChargedLeptons[theZ.l2Idx]
                    lepW = next((l for l in BestLeptons if l not in (lep1, lep2)), None)
                    self.hEvent.fill1D("Wlep_pt", "Wlep_pt", 120, 0., 120., lepW.pt, self.weight)
                else : return False, 4.

                
                # ricostruire la massa trasversa di W, mi sa che posso solo vedere il limite superiore per info su MET che ho
                Mt_WMax = 2*math.sqrt(lepW.pt*MET_pt)
                self.hEvent.fill1D("Mt_W_Max", "Mt_W_Max", 120, 0., 150.,  Mt_WMax, self.weight)
                

                notpassedZ, invMassLLGamma = self.isPhotonFSR([lep1,lep2], bestPhoton, 110)
                self.hEvent.fill1D("llGamma_invMass", "llGamma_invMass", 100, 40., 250., invMassLLGamma, self.weight)
                if notpassedZ : return False, 5.
                    

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
        
        #come scelgo poi il fotone migliore? Per ora prendo quello con pt maggiore (dovrebbe venire dal coupling)
        if len(BestPhotons) > 0 :
            return max(BestPhotons, key=lambda p: p.pt, default=None)
        else : return None

        
        
  
    #--------------------------------------
    # Functions for WZ RECOSTRUCTION
    #--------------------------------------
    # given a collection of leptons (charged and neutrino), this function reconstructs the W or the Z by the pdgId given
    def BosonRecon(self, Leptons, partPdgId,  massThreshold, chooseNeutrinoByPt) :
        CandMass = {}
        MassPdg = TDatabasePDG.Instance().GetParticle(partPdgId).Mass()

        GoodLeptons = Leptons
        if chooseNeutrinoByPt and partPdgId == 24  :
            GoodLeptons = self.getLeptonList(Leptons)
               
        for i in range(len(GoodLeptons)):
            for j in range(i + 1, len(GoodLeptons)):
                if (partPdgId == 24 and abs(GoodLeptons[i].pdgId+ GoodLeptons[j].pdgId) == 1) :   
                    CandMass[(i,j)] = self.getInvMass(GoodLeptons[i], GoodLeptons[j])  
                elif (partPdgId == 23 and GoodLeptons[i].pdgId+ GoodLeptons[j].pdgId == 0) :
                     CandMass[(i,j)] = self.getInvMass(GoodLeptons[i], GoodLeptons[j])  
                else : CandMass[(i,j)] = None
                    
        bestCouple = self.getCloserCand(MassPdg, CandMass)
        if(bestCouple != None and massThreshold[0] < CandMass[bestCouple] < massThreshold[1]) :  
            if(partPdgId == 24) : boson = GenW([GoodLeptons[bestCouple[0]], GoodLeptons[bestCouple[1]]])
            elif(partPdgId == 23) : boson = GenZ([GoodLeptons[bestCouple[0]], GoodLeptons[bestCouple[1]]])
            return boson, bestCouple

        return None, None


    # return a W and a Z boson. It reconstruct them in order to minimize the difference with the theoretical mass
    def WZRecon(self, ChLeptons, Neutrinos, useResidues, chooseNeutrinoByPt) :
        if any(l.pdgId < 0 for l in ChLeptons) and any(l.pdgId > 0 for l in ChLeptons)  :  
            Z1, coupleZ1 = self.BosonRecon(ChLeptons, 23, (60,120), False)
            W1, coupleW1 = None, None
            if coupleZ1 != None :
                RemainingChLeptons1 = [p for i, p in enumerate(ChLeptons) if i not in coupleZ1]
                Leptons1 = RemainingChLeptons1 + Neutrinos
                W1, coupleW1 = self.BosonRecon(Leptons1, 24, (50,110), chooseNeutrinoByPt)
            
            ret = (W1, Z1)
            if useResidues : 
                Leptons2 = ChLeptons + Neutrinos
                W2, coupleW2 = self.BosonRecon(Leptons2, 24, (50,110), chooseNeutrinoByPt)
                if coupleW2 != None : 
                    RemainingChLeptons2 = [p for i, p in enumerate(ChLeptons) if i not in coupleW2]
                    Z2, coupleZ2 = self.BosonRecon(RemainingChLeptons2, 23, (60,120), False)

                    WRecIsBest = self.getLowerRes(W1,Z1,W2,Z2)
                    if WRecIsBest :  ret = (W2, Z2)
                        
            return ret
                
        return None, None


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
    # Return True if the residue between the gen boson masses and the theoretically one is lower for the second pair
    def getLowerRes(self, Wboson1, Zboson1, Wboson2, Zboson2) :
        WMassPdg = TDatabasePDG.Instance().GetParticle(24).Mass()
        ZMassPdg = TDatabasePDG.Instance().GetParticle(23).Mass()
        
        res = lambda b, mass: abs(b.mass - mass) if b is not None else float('inf')
        res1 = res(Wboson1, WMassPdg) + res(Zboson1, ZMassPdg)
        res2 = res(Wboson2, WMassPdg) + res(Zboson2, ZMassPdg)
        if res1 > res2 : return True
        else : return False
        
    
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
    def getInvMass(self, *particles) :
        ptot = ROOT.TLorentzVector()
        for i, part in enumerate(particles):
            if i == 2:  
                p4_temp = ROOT.TLorentzVector()
                p4_temp.SetPtEtaPhiM(part.pt, part.eta, part.phi, 0.0)
                ptot += p4_temp
            else:
                ptot += part.p4()
                
        return ptot.M()


    # extrapolate from a list of leptons the charged ones plus the neutrino with the biggest pt 
    def getLeptonList(self, Leptons) :
        ChLeptons = [p for p in Leptons if abs(p.pdgId) in (11, 13)]       
        chosenNeutrino = max((p for p in Leptons if abs(p.pdgId) in (12, 14)),key=lambda p: p.pt,default=None,)

        GoodLeptons = ChLeptons.copy()
        if chosenNeutrino is not None:
            GoodLeptons.append(chosenNeutrino)
        return GoodLeptons

        

#---------------------------
# other CLASSES
#---------------------------
#per il momento non tengo conto della FINAL STATE RADIATION (fotoni emessi da leptoni) 
class GenZ:   
    def __init__(self, Leptons):
        self.leptons = Leptons

        self.p4 = Leptons[0].p4() + Leptons[1].p4()
        self.mass = self.p4.M()
        self.pt = self.p4.Pt()
        self.eta = self.p4.Eta()
        self.phi = self.p4.Phi()  


class GenW:
    def __init__(self, Leptons):
        self.leptons = Leptons

        self.p4 = Leptons[0].p4() + Leptons[1].p4()
        self.mass = self.p4.M()
        self.pt = self.p4.Pt()
        self.eta = self.p4.Eta()
        self.phi = self.p4.Phi() 
