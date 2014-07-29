#include "VVXAnalysis/TreeAnalysis/interface/VVXAnalyzer.h"
#include "VVXAnalysis/Commons/interface/SignalDefinitions.h"
#include "VVXAnalysis/Commons/interface/Utilities.h"

#include <boost/foreach.hpp>
#define foreach BOOST_FOREACH

using std::cout;
using std::endl;


using namespace phys;

Int_t VVXAnalyzer::cut() {
  
  if(!ZZ->passTrigger())          theHistograms.fill("CutCheck",5,0,5,0);
  if(! 60 < ZZ->first().mass())   theHistograms.fill("CutCheck",5,0,5,1);
  if(! ZZ->first().mass() < 120)  theHistograms.fill("CutCheck",5,0,5,2);
  if(! 60 < ZZ->second().mass())  theHistograms.fill("CutCheck",5,0,5,3);
  if(! ZZ->second().mass() < 120) theHistograms.fill("CutCheck",5,0,5,4);
  
  if(ZZ->passTrigger() &&
     60 < ZZ->first().mass()  && ZZ->first().mass()  < 120 &&
     60 < ZZ->second().mass() && ZZ->second().mass() < 120)
    return 1;
  
  theHistograms.fill("mZ1_nonPassing", "Invariant mass of Z_{1} for non passing events",  50, 0, 200, ZZ->first().mass() , theWeight); 
  theHistograms.fill("mZ2_nonPassing", "Invariant mass of Z_{2} for non passing events",  50, 0, 200, ZZ->second().mass(), theWeight); 

  return -1;
}

void VVXAnalyzer::ZZplots(int id){

  if(ZZ->id() != id && id != -1) return; // -1 here means generic 4l final state

  std::string decay  = "4l";
  std::string decay1 = "l";
  std::string decay2 = "l";
  if      (id == 52) {decay = "4m"  ; decay1 = "m"; decay1 = "m";}
  else if (id == 48) {decay = "2e2m"; decay1 = "e"; decay2 = "m";}
  else if (id == 44) {decay = "4e"  ; decay1 = "e"; decay2 = "e";}
  theHistograms.fill(std::string("ZZTo")+decay+std::string("_mZ1To2")+decay1, std::string("Invariant mass of Z_{1}#rightarrow 2")+decay1,  15, 60,  120, ZZ->first().mass() , theWeight); 
  theHistograms.fill(std::string("ZZTo")+decay+std::string("_mZ2To2")+decay2, std::string("Invariant mass of Z_{2}#rightarrow 2")+decay2,  15, 60,  120, ZZ->second().mass(), theWeight); 
  theHistograms.fill(std::string("ZZTo")+decay+std::string("_mZZTo") +decay , std::string("Invariant mass of ZZ#rightarrow ")    +decay ,  40,  0, 1000, ZZ->mass()         , theWeight); 

  theHistograms.fill(std::string("ZZTo")+decay+"_nJets"       , "Number of jets (|#eta|<4.7 and p_T > 30 GeV)"        , 10, 0, 10, jets->size(), theWeight); 
  theHistograms.fill(std::string("ZZTo")+decay+"_nCentralJets", "Number of central jets (|#eta|<2.5 and p_T > 30 GeV)", 10, 0, 10, centralJets->size(), theWeight); 

  if(jets->size() >= 2)
    theHistograms.fill(std::string("ZZTo")+decay+"_deltaEtaJJ", "#Delta #eta(j,j) between the two most energetic jets",  10, 0, 8, fabs(jets->at(0).eta() - jets->at(1).eta()), theWeight); 
  
  
  if(centralJets->size() >= 2){
    theHistograms.fill(std::string("ZZTo")+decay+"_deltaEtaJJcentral", "#Delta #eta(j,j) between the two most energetyc central jets",  10, 0, 8, fabs(centralJets->at(0).eta() - centralJets->at(1).eta()), theWeight); 
    theHistograms.fill(std::string("ZZTo")+decay+"_mJJ", "m_{jj}",  20, 0, 1000, (centralJets->at(0).p4() + centralJets->at(1).p4()).M(), theWeight); 
  }

  theHistograms.fill(std::string("ZZTo")+decay+"_nExtraMuons"    , "Number of extra muons in the event"    , 10, 0, 10, muons->size(), theWeight); 
  theHistograms.fill(std::string("ZZTo")+decay+"_nExtraElectrons", "Number of extra electrons in the event", 10, 0, 10, electrons->size(), theWeight); 
  theHistograms.fill(std::string("ZZTo")+decay+"_nExtraLeptons"  , "Number of extra leptons in the event"  , 10, 0, 10, muons->size()+electrons->size(), theWeight); 

}


void VVXAnalyzer::analyze(){
  
  // Some basic plots on ZZ
  ZZplots();   // ZZ --> 4l
  ZZplots(52); // ZZ --> 4m
  ZZplots(48); // ZZ --> 2e2m
  ZZplots(44); // ZZ --> 4e
}



