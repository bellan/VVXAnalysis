from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection
from ZZAnalysis.NanoAnalysis.tools import getLeptons, get_genEventSumw

from VVXAnalysis.NanoAnalysis.EventAnalyzer import EventAnalyzer
from VVXAnalysis.NanoAnalysis.Histogrammer import *
from VVXAnalysis.NanoAnalysis.Regions import Flags as Regions

from enum import IntFlag

from ROOT import TDatabasePDG
import itertools

class VZGammaAnalyzer(EventAnalyzer, analysis_name="VZGammaAnalyzer"):

    
    def __init__(self, regions):
        super().__init__(regions)
                    
        
    def analyze(self):
        EtichetteTagliSignalDef=["Total Events","After GenStatusFlag","After NQuarks","LeptonCut","AfterPhotonPtEta", "AfterJetCuts"]
        EtichetteTagliSignalReg=["Total Events","Trigger","LeptonCut","After electronVeto and WP90","AfterPhotonPtEta", "JetIdLoose", "AfterKinAndMassJetCuts"]
        self.hEvent.fill1D_label("Eventi per taglio (SD)","Eventi per taglio (SD)",EtichetteTagliSignalDef,EtichetteTagliSignalDef[0],self.weight)
        self.hEvent.fill1D_label("Eventi per taglio (SR)","Eventi per taglio (SR)",EtichetteTagliSignalReg,EtichetteTagliSignalReg[0],self.weight)


      #===========
      #Collection
      #===========
        if self.analyzeMC:
            GenParts=Collection(self.event,'GenPart')
            GenPhotons=[p for p in GenParts if p.pdgId==22 and p.status==1 and (p.statusFlags & GenStatusFlag.SEL)==GenStatusFlag.SEL ]
            GenLeptons=[p for p in GenParts if (abs(p.pdgId)==11 or abs(p.pdgId)==13) and p.status == 1 and (p.statusFlags & GenStatusFlag.SEL) == GenStatusFlag.SEL]
            GenJets=Collection(self.event,'GenJet')
            GenQuark=[p for p in GenParts if ((0<abs(p.pdgId)<7) and p.status==23 and (p.statusFlags & GenStatusFlag.QUARK_SEL) == GenStatusFlag.QUARK_SEL)]

        Photons=Collection(self.event,'Photon')
        Leptons=Collection(self.event,'Lepton')
        Zs=Collection(self.event, 'ZCand')
        Jets=Collection(self.event, 'Jet')
        bestZIdx = self.event.bestZIdx    
            
        #========================
        #N. particelle Signal Def
        #========================
        self.hEvent.fill1D("nGenLeptonsPostStatusFlag", "nGenLeptonsPostStatusFlag", 93, 0., 10., len(GenLeptons), self.weight)
        self.hEvent.fill1D("nGenPhotonsPostStatusFlag", "nGenPhotonsPostStatusFlag", 93, 0., 10., len(GenPhotons), self.weight)
        self.hEvent.fill1D("nGenJetsPostStatusFlag", "nGenJetsPostStatusFlag", 93, 0., 10., len(GenJets), self.weight)

        #========================
        #N. particelle Signal Reg
        #========================
        self.hEvent.fill1D("nLeptons", "nLeptons", 93, 0., 10., len(Leptons), self.weight)
        self.hEvent.fill1D("nPhotons", "nPhotons", 93, 0., 10., len(Photons), self.weight)
        self.hEvent.fill1D("nJets", "nJets", 93, 0., 10., len(Jets), self.weight)

        

        #=================
        #SIGNAL DEFINITION
        #=================  

        def isInSignalDef():

              bestjetpair=[]
              if not self.analyzeMC: return False, bestjetpair
              
              self.hEvent.fill1D_label("Eventi per taglio (SD)","Eventi per taglio (SD)",EtichetteTagliSignalDef,EtichetteTagliSignalDef[1],self.weight)

              if len(GenQuark) != 2: return False, bestjetpair

              self.hEvent.fill1D_label("Eventi per taglio (SD)","Eventi per taglio (SD)",EtichetteTagliSignalDef,EtichetteTagliSignalDef[2],self.weight)
              
              ptEtaCutLeptons=self.leptonPtEtaCut(GenLeptons)
              if self.leptonIsSignal(ptEtaCutLeptons)==False: return False, bestjetpair
              if self.leptonsMassCut(ptEtaCutLeptons)==False: return False, bestjetpair

              #=====================================
              #Ricostruire Z a partire dai 2 leptoni
              #=====================================
      
              self.hEvent.fill1D_label("Eventi per taglio (SD)","Eventi per taglio (SD)",EtichetteTagliSignalDef,EtichetteTagliSignalDef[3],self.weight)
              Zmass=(ptEtaCutLeptons[0].p4()+ptEtaCutLeptons[1].p4()).M()
              self.hEvent.fill1D("ZMass","Zmass",45,30,155,Zmass,self.weight)

              #=======
              # Fotoni
              # ======
              
              
              ptEtaCutPhoton=self.photonPtEtaCut(GenPhotons)
              if len(ptEtaCutPhoton)==0: return False, bestjetpair
              self.hEvent.fill1D_label("Eventi per taglio (SD)","Eventi per taglio (SD)",EtichetteTagliSignalDef,EtichetteTagliSignalDef[4],self.weight) 

              #======
              #Adroni
              #======
              ptEtaCutJets=self.jetPtEtaCut(GenJets)
              
              #=========================
              #Ricostruire m_inv dei jet
              #=========================
              #senza la selezione
              '''
              for i, jet1 in enumerate(finalCutJets):
                     for n in range(i+1,len(finalCutJets)):
                        jet2=finalCutJets[n]
                        jetMass=(jet1.p4()+jet2.p4()).M()
                        self.hEvent.fill1D("JetMass","Jetmass",110,30,150,jetMass,self.weight)
              
                  '''
              deltaRCutJets=self.jetDeltaRCut(ptEtaCutJets,ptEtaCutPhoton,ptEtaCutLeptons) 
              if (self.bestJetPair(deltaRCutJets)==False): return False, bestjetpair
              bestjetpair=self.bestJetPair(deltaRCutJets)
              #if (len(bestjetpair)!=2):  return False
              finalJetmass=(bestjetpair[0].p4()+bestjetpair[1].p4()).M()
              self.hEvent.fill1D("JetMassAfterSelection","JetmassAfterSelection",110,30,150,finalJetmass,self.weight)
              self.hEvent.fill1D_label("Eventi per taglio (SD)","Eventi per taglio (SD)",EtichetteTagliSignalDef,EtichetteTagliSignalDef[5],self.weight)
              
              #ricostruzione mlljj

              mlljj=(ptEtaCutLeptons[0].p4()+ptEtaCutLeptons[1].p4()+bestjetpair[0].p4()+bestjetpair[1].p4()).M()
              self.hEvent.fill1D("mlljj","mlljj",110,120,650,mlljj,self.weight)

              #ricostruzione mlljjgamma

              leadingPhoton=max(ptEtaCutPhoton, key=lambda ph:ph.pt)
              mlljjg=(ptEtaCutLeptons[0].p4()+ptEtaCutLeptons[1].p4()+bestjetpair[0].p4()+bestjetpair[1].p4()+leadingPhoton.p4()).M()
              self.hEvent.fill1D("mlljjg","mlljjg",110,120,850,mlljjg,self.weight)

              return True, bestjetpair

        def isInSignalReg():
              bestPairSR=[]  
              if not self.event.HLT_passZZ4l : return False,bestPairSR
              self.hEvent.fill1D_label("Eventi per taglio (SR)","Eventi per taglio (SR)",EtichetteTagliSignalReg,EtichetteTagliSignalReg[1],self.weight)
              #=========
              #Leptoni
              #=========

              ptEtaCutLeptonsSR=self.leptonPtEtaCut(Leptons)
              if self.leptonIsSignal(ptEtaCutLeptonsSR)==False: return False,bestPairSR
              if self.leptonsMassCut(ptEtaCutLeptonsSR)==False: return False,bestPairSR
              self.hEvent.fill1D_label("Eventi per taglio (SR)","Eventi per taglio (SR)",EtichetteTagliSignalReg,EtichetteTagliSignalReg[2],self.weight)

              #==========
              #Fotoni
              #==========
              truePhotons=self.photonIdCut(Photons)
              passedVeto,vetoPh=self.electronVeto(truePhotons)
              if (passedVeto==False): return False,bestPairSR
              self.hEvent.fill1D_label("Eventi per taglio (SR)","Eventi per taglio (SR)",EtichetteTagliSignalReg,EtichetteTagliSignalReg[3],self.weight)
              ptEtaCutPhotonSR=self.photonPtEtaCut(vetoPh)
              if len(ptEtaCutPhotonSR)==0: return False,bestPairSR

              self.hEvent.fill1D_label("Eventi per taglio (SR)","Eventi per taglio (SR)",EtichetteTagliSignalReg,EtichetteTagliSignalReg[4],self.weight)

              #=========
              #Adroni
              #=========
 
              JetIdSuccess,JetsId=self.JetIdCut(Jets)
              if (JetIdSuccess== False): return False,bestPairSR
              self.hEvent.fill1D_label("Eventi per taglio (SR)","Eventi per taglio (SR)",EtichetteTagliSignalReg,EtichetteTagliSignalReg[5],self.weight)
              
              ptEtaCutJetsSR=self.jetPtEtaCut(JetsId)
              deltaRCutJetsSR=self.jetDeltaRCut(ptEtaCutJetsSR,ptEtaCutPhotonSR,ptEtaCutLeptonsSR)

              if (self.bestJetPair(deltaRCutJetsSR)==False): return False,bestPairSR

              bestPairSR=self.bestJetPair(deltaRCutJetsSR)
              finalJetmassSR=(bestPairSR[0].p4()+bestPairSR[1].p4()).M()
              self.hEvent.fill1D("JetMassAfterSelection(SR)","JetmassAfterSelection(SR)",110,30,150,finalJetmassSR,self.weight)
              self.hEvent.fill1D_label("Eventi per taglio (SR)","Eventi per taglio (SR)",EtichetteTagliSignalReg,EtichetteTagliSignalReg[6],self.weight)
              
              return True,bestPairSR

        
        

        res,bestPairSR=isInSignalDef()
        if res==True:res1,bestjetpair=isInSignalReg()
       

        #check ricostruzione
        if(res==True and res1==True):
              success=0
              if (bestPairSR[0].DeltaR(bestjetpair[0])<0.4):
                    if(bestPairSR[1].DeltaR(bestjetpair[1])<0.4): success=1
              if (bestPairSR[0].DeltaR(bestjetpair[1])<0.4):
                    if(bestPairSR[1].DeltaR(bestjetpair[0])<0.4): success=1
              self.hEvent.fill1D("RicostruzioneJet", "RicostruzioneJet",2 , -0.5, 1.5, success, self.weight)     


        
            
    

            
              
          

    #=====================
    # SIGNAL DEF FUNCTIONS
    #=====================

    def photonPtEtaCut(self,Photons):
            CutPhotons=[p for p in Photons if abs(p.eta)<2.4 and p.pt>5 ]   
            return CutPhotons
              
    def leptonPtEtaCut(self,Leptons):
           firstCutLeptons=[p for p in Leptons if p.pt>10 and abs(p.eta)<2.5]
           return firstCutLeptons           

    def leptonIsSignal(self,Leptons):
           if len(Leptons)!=2: return False
           if (Leptons[0].pdgId!=-Leptons[1].pdgId): return False
           if not (Leptons[0].pt>20 or Leptons[1].pt>20): return False
           return True

    def leptonsMassCut(self,Leptons):
          Z_inv_mass=(Leptons[0].p4()+Leptons[1].p4()).M()
          if (50<Z_inv_mass<120 and Leptons[0].DeltaR(Leptons[1])>0.5): return True
          else: return False
         

                 
    def jetPtEtaCut(self,Jets):
          cutJets=[j for j in Jets if abs(j.eta)<4.7 and j.pt>30]
          return cutJets
    
    def jetMassCut(self,Jets):
          jetMassCutPassed=[]
          if (len(Jets)<2): return False
          for jet1,jet2 in itertools.combinations(Jets,2):
            jetMass=(jet1.p4()+jet2.p4()).M()
            if (50<jetMass<120): 
                  jetMassCutPassed.append(jet1)
                  jetMassCutPassed.append(jet2)
                                  
          jetMassCutPassed=list(set(jetMassCutPassed))  
          if len(jetMassCutPassed)==0: return False  
          return jetMassCutPassed

        
    def jetDeltaRCut(self,Jets,Photons,Leptons):
          jetDeltaRPassed=[p for p in Jets if all(p.DeltaR(l)>0.5 for l in Leptons) and all(p.DeltaR(g)>0.5 for g in Photons)]
          return jetDeltaRPassed

    def bestJetPair(self,Jets):
          if (len(Jets)<2): return False
          bestPair=[]
          best_min=50
          WMassPDG = TDatabasePDG.Instance().GetParticle(24).Mass()
          ZMassPDG = TDatabasePDG.Instance().GetParticle(23).Mass()
          for j1,j2 in itertools.combinations(Jets,2):
                jetMass=(j1.p4()+j2.p4()).M()
                if (50<jetMass<120):
                  W_diff=abs(WMassPDG-jetMass)
                  Z_diff=abs(ZMassPDG-jetMass)
                  current_min=min(W_diff,Z_diff)
                  if (current_min<best_min):
                        best_min=current_min
                        bestPair=[j1,j2]
          if len(bestPair)==0: return False         
          return bestPair

                      
   #=====================
   # SIGNAL REG FUNCTIONS
   # ====================
                
    def electronVeto(self, Photons):
         if len(Photons)<1: return False, None
         PhotonsVetoPass=[ph for ph in Photons if (ph.cutBased >= 1 and ph.electronVeto)]
         if len(PhotonsVetoPass)<1: return False, None
         return True, PhotonsVetoPass

    def JetIdCut(self,Jets):
         if len(Jets)<1: return False,None
         JetsIdPass=[j for j in Jets if j.jetId>=2]
         if (len(JetsIdPass)==0):return False,None
         return True,JetsIdPass

    def photonIdCut(self,Photons, Id="MVA"):
         if Id=="MVA":
              truePhotons=[ph for ph in Photons if (ph.mvaID_WP90)]
         return truePhotons
         
         

          
          
                      
            
                
class GenStatusFlag(IntFlag):
      IS_PROMPT = 1 << 0
      FROM_HARD_PROCESS = 1 << 8
      IS_FIRST_COPY = 1 << 12

      SEL = IS_PROMPT | FROM_HARD_PROCESS
      QUARK_SEL = IS_PROMPT | FROM_HARD_PROCESS | IS_FIRST_COPY
      #PH_SEL = IS_PROMPT                               


    
          


           
           
                      
                  


 
#ORIDINE TAGLI
#hlt trigger
#parte leptonica:
#tagli cinematici, 2 leptoni
#ricostruire la Z, regione di massa 60-120
#Selezione sui fotoni: pt,eta,(electronVeto),Id
#Disambiguazione sui jet (rimuovo i jet sovrapposti con i fotoni che hanno un grande parte del pt del jet) (signal region)
#separazione dei fotoni dagli altri oggetti
#più di due jet con tagli cinematici, massa inv più vicina a W o Z
#ricostr mllg
#Introduzione
#grafici cutflow in Sr e Sd
#tagli e varibili (descrivere gli algoritmi come ho ricostruito i jet, pulizia con il pt)
#significanza
#disambiguazione in DeltaR
#electronVeto, jetId, PhotonId anche gli eventi che non hanno una coppia di jet nella def di segnale vanno messi a denominatore
#controllare se mancano dei tagli 
           
   





           
           
                      
                  

