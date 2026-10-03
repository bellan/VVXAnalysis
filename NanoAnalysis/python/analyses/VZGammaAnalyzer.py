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
        EtichetteTagliSignalDef=["Total Events","After Trigger","AfterPhotonKin","AfterPhotonDelta","AfterLeptonKin"]
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
              DeltaRCutph=self.photonDeltaRCut(kinCutph,kinCutlep,GenJets)
              if len(DeltaRCutph)==0: return False
              self.hEvent.fill1D_label("Eventi per taglio","Eventi per taglio",EtichetteTagliSignalDef,EtichetteTagliSignalDef[3],self.weight)
        

        res=IsSignal()  
            
    

            
              
          

 

    def photonKinCut(self,Photons):
            CutPhotons=[p for p in Photons if abs(p.eta)<2.4 and p.pt>20 ]   
            return CutPhotons

    def photonDeltaRCut(self,Photons,Leptons,Jets):
           secondCutPhotons=[p for p in Photons if all(p.DeltaR(l)>0.5 for l in Leptons) and all(p.DeltaR(j)>0.5 for j in Jets)]
           return secondCutPhotons                 

    def leptonKinCut(self,Leptons):
           firstCutLeptons=[p for p in Leptons if p.pt>10 and abs(p.eta)<2.5]
           return firstCutLeptons           
    
    def ZLeptons(self,Leptons):
        leptons=[p for p in Leptons if p.pt>20]
        lepton1=[p for p in leptons if all(leptons.DeltaR(l)>0.5 for l in leptons)]
        lepton2=[p for p in lepton1 if all(lepton1.DeltaR(l)>0.5 for l in lepton1)]
    #def invMassCut(self,Leptons,Photons):
           
                      
                  
class GenStatusFlag(IntFlag):
                IS_PROMPT = 1 << 0
                FROM_HARD_PROCESS = 1 << 8

                SEL = IS_PROMPT | FROM_HARD_PROCESS
                #PH_SEL = IS_PROMPT
