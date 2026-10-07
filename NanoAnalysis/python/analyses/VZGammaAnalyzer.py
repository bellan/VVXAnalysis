from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection
from ZZAnalysis.NanoAnalysis.tools import getLeptons, get_genEventSumw

from VVXAnalysis.NanoAnalysis.EventAnalyzer import EventAnalyzer
from VVXAnalysis.NanoAnalysis.Histogrammer import *
from VVXAnalysis.NanoAnalysis.Regions import Flags as Regions

from enum import IntFlag

class VZGammaAnalyzer(EventAnalyzer, analysis_name="VZGammaAnalyzer"):

    
    def __init__(self, regions):
        super().__init__(regions)
                    
        
    def analyze(self):
        EtichetteTagliSignalDef=["Total Events","After Trigger","AfterPhotonKin","AfterLepKin","AfterJetKin","AfterLepInvMassandDeltaR"]
        self.hEvent.fill1D_label("Eventi per taglio","Eventi per taglio",EtichetteTagliSignalDef,EtichetteTagliSignalDef[0],self.weight)
        if(self.event.HLT_passZZ4l): 

            #===========
            #Collection
            #===========

            GenParts=Collection(self.event,'GenPart')
            GenPhotons=[p for p in GenParts if p.pdgId==22 and p.status==1 and (p.statusFlags & GenStatusFlag.SEL)==GenStatusFlag.SEL ]
            GenLeptons=[p for p in GenParts if (abs(p.pdgId==11) or abs(p.pdgId==13)) and p.status == 1 and (p.statusFlags & GenStatusFlag.SEL) == GenStatusFlag.SEL]
            GenJets=Collection(self.event,'GenJet')

            Photons=Collection(self.event,'Photon')
            Leptons=Collection(self.event,'Lepton')
            Zs=Collection(self.event, 'ZCand')
            Jets=Collection(self.event, 'Jet')
            
            
            
            self.hEvent.fill1D("nGenLeptonsPostStatusFlag", "nGenLeptonsPostStatusFlag", 93, 0., 10., len(GenLeptons), self.weight)
            self.hEvent.fill1D("nGenPhotonsPostStatusFlag", "nGenPhotonsPostStatusFlag", 93, 0., 10., len(GenPhotons), self.weight)
            self.hEvent.fill1D("nGenJetsPostStatusFlag", "nGenJetsPostStatusFlag", 93, 0., 10., len(GenJets), self.weight)

            self.hEvent.fill1D_label("Eventi per taglio","Eventi per taglio",EtichetteTagliSignalDef,EtichetteTagliSignalDef[1],self.weight)

        #=================
        #SIGNAL DEFINITION
        #=================  

        def IsSignal():
              kinCutph=self.photonKinCut(GenPhotons)
              if len(kinCutph)==0: return False
              self.hEvent.fill1D_label("Eventi per taglio","Eventi per taglio",EtichetteTagliSignalDef,EtichetteTagliSignalDef[2],self.weight)
              kinCutlep=self.leptonKinCut(GenLeptons)
              if len(kinCutlep)==0: return False
              self.hEvent.fill1D_label("Eventi per taglio","Eventi per taglio",EtichetteTagliSignalDef,EtichetteTagliSignalDef[3],self.weight)
              kinCutJets=self.JetkinCut(GenJets)
              if len(kinCutJets)==0: return False
              self.hEvent.fill1D_label("Eventi per taglio","Eventi per taglio",EtichetteTagliSignalDef,EtichetteTagliSignalDef[4],self.weight)
              if self.leptonIsSignal(kinCutlep)==False: return False
              self.hEvent.fill1D_label("Eventi per taglio","Eventi per taglio",EtichetteTagliSignalDef,EtichetteTagliSignalDef[5],self.weight)
              
        

        res=IsSignal()  
            
    

            
              
          

    #=====================
    # SIGNAL DEF FUNCTIONS
    #=====================

    def photonKinCut(self,Photons):
            CutPhotons=[p for p in Photons if abs(p.eta)<2.4 and p.pt>15 ]   
            return CutPhotons

    def photonDeltaRCut(self,Photons,Leptons,Jets):
           secondCutPhotons=[p for p in Photons if all(p.DeltaR(l)>0.5 for l in Leptons) and all(p.DeltaR(j)>0.5 for j in Jets)]
           return secondCutPhotons                 

    def leptonKinCut(self,Leptons):
           firstCutLeptons=[p for p in Leptons if p.pt>5 and abs(p.eta)<2.5 and not 1.444<abs(p.eta)<1.566]
           return firstCutLeptons           

    def leptonIsSignal(self,Leptons):
           if len(Leptons)!=2: return False #mi servono 2 leptoni
           #if (sum(l.pdgId == 11 for l in Leptons) != sum(l.pdgId == -11 for l in Leptons)
           #or sum(l.pdgId == 13 for l in Leptons) != sum(l.pdgId == -13 for l in Leptons)): return False #di segno opposto
           Z_inv_mass=(Leptons[0].p4()+Leptons[1].p4()).M()
           if (50<Z_inv_mass<120 and Leptons[0].DeltaR(Leptons[1])>0.02): return True #taglio su massa inv e delta R
           return True
           
    def JetkinCut(self,Jets):
          cutJets=[j for j in Jets if abs(j.eta)<4.7 and j.pt>30]
          return cutJets





           
           
                      
                  
class GenStatusFlag(IntFlag):
                IS_PROMPT = 1 << 0
                FROM_HARD_PROCESS = 1 << 8

                SEL = IS_PROMPT | FROM_HARD_PROCESS
                #PH_SEL = IS_PROMPT
