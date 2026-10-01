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
        if(self.event.HLT_passZZ4l): 
            weight = 1.

            GenParts=Collection(self.event,'GenPart')
            GenPhotons=[p for p in GenParts if p.pdgId==22 and p.status==1 and (p.statusFlags & GenStatusFlag.SEL)==GenStatusFlag.SEL ]
            GenLeptons=[p for p in GenParts if (abs(p.pdgId==11) or abs(p.pdgId==13)) and p.status == 1 and (p.statusFlags & GenStatusFlag.SEL) == GenStatusFlag.SEL]
            GenJets=Collection(self.event,'GenJet')

            Photons=Collection(self.event,'Photon')
            Leptons=Collection(self.event,'Lepton')
            Zs=Collection(self.event, 'ZCand')
            
            if self.analyzeMC: self.weight = (self.event.overallEventWeight/self.genEventSumw)      #da rivedere
            BestGamma=max(self.photonCut(GenPhotons,GenLeptons,GenJets),key=lambda ph: ph.pt, default=None)   
            if BestGamma!=None and len(BestGamma)==1:
                   self.hEvent.fill1D("BestGammaPt", "BestGammaPt", 50, 20., 220., BestGamma.pt, self.weight)
                
            
        
            
            
    

    def SignalDefinition(self,GenPhotons,GenLeptons):
                BestGamma=max(self.photonCut(GenPhotons,GenLeptons),key=lambda ph: ph.pt, default=None)   

    
    def photonCut(self,Photons,Leptons,Jets):
            firstCutPhotons=[p for p in Photons if abs(p.eta)<2.4 and p.pt>20 ]
            secondCutPhotons=[p for p in firstCutPhotons if all(firstCutPhotons.DeltaR(l)>0.5 for l in Leptons) and all(firstCutPhotons.DeltaR(j)>0.5 for j in Jets)]    
            return secondCutPhotons              

    def leptonCut(self,Leptons):
           firstCutLeptons=[p for p in Leptons if p.pt>10]
           secondCutLeptons=[p for p in firstCutLeptons if abs(p.eta)<2.5]
    def ZLeptons(self,Leptons):
        leptons=[p for p in Leptons if p.pt>20]
        lepton1=[p for p in leptons if all(leptons.DeltaR(l)>0.5 for l in leptons)]
        lepton2=[p for p in lepton1 if all(lepton1.DeltaR(l)>0.5 for l in lepton1)]
               
                  
class GenStatusFlag(IntFlag):
                IS_PROMPT = 1 << 0
                FROM_HARD_PROCESS = 1 << 8

                SEL = IS_PROMPT | FROM_HARD_PROCESS
                #PH_SEL = IS_PROMPT

