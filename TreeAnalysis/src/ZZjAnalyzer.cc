#include "VVXAnalysis/TreeAnalysis/interface/ZZjAnalyzer.h"
#include "VVXAnalysis/Commons/interface/SignalDefinitions.h"
#include "VVXAnalysis/Commons/interface/Utilities.h"

#include <boost/foreach.hpp>
#define foreach BOOST_FOREACH

using std::cout;
using std::endl;


using namespace phys;

// Int_t VVXAnalyzer::cut() {
  
//   return 1;
// }

void ZZjAnalyzer::ZZplots(int id){

  if(ZZ->id() != id && id != -1) return; // -1 here means generic 4l final state

  std::string decay  = "4l";
  
  if      (id == 52) {decay = "4m";}
  else if (id == 48) {decay = "2e2m";}
  else if (id == 44) {decay = "4e";}

  theHistograms.fill(std::string("ZZTo")+decay+"_Mass"         , std::string("Invariant mass of ZZ_{1}#rightarrow ")+decay            ,  200, 50,  1000, ZZ->mass(),theWeight);
  theHistograms.fill(std::string("ZZTo")+decay+"_Mass"+"_FRVar", std::string("Var From FR Invariant mass of ZZ_{1}#rightarrow ")+decay,  200, 50,  1000, ZZ->mass(),ZZ->fakeRateSFVar());
 
 
  theHistograms.fill(std::string("ZZTo")+decay+"_nJets"      , "Number of jets (|#eta|<4.7 and p_T > 30 GeV)"            , 10, 0, 10, jets->size(), theWeight); 
  theHistograms.fill(std::string("ZZTo")+decay+"_nJets_FRVar", "Var From FR Number of jets (|#eta|<4.7 and p_T > 30 GeV)", 10, 0, 10, jets->size(), ZZ->fakeRateSFVar());
  

  theHistograms.fill(std::string("ZZTo")+decay+"_nCentralJets"      , "Number of central jets (|#eta|<2.5 and p_T > 30 GeV)"            , 10, 0, 10, centralJets->size(), theWeight); 
  theHistograms.fill(std::string("ZZTo")+decay+"_nCentralJets_FRVar", "Var from FR Number of central jets (|#eta|<2.5 and p_T > 30 GeV)", 10, 0, 10, centralJets->size(), ZZ->fakeRateSFVar());

 
  if(jets->size() >= 2) {
    theHistograms.fill(std::string("ZZTo")+decay+"_deltaEtaJJ"      , "#Delta #eta(j,j) between the two most energetic jets"            ,  10, 0, 8, fabs(jets->at(0).eta() - jets->at(1).eta()), theWeight); 
    theHistograms.fill(std::string("ZZTo")+decay+"_deltaEtaJJ_FRVar", "Var From FR #Delta #eta(j,j) between the two most energetic jets",  10, 0, 8, fabs(jets->at(0).eta() - jets->at(1).eta()),ZZ->fakeRateSFVar());

    theHistograms.fill(std::string("ZZTo")+decay+"_deltaYJJ"      , "#Delta Y(j,j) between the two most energetic jets"            ,  10, 0, 8, fabs(jets->at(0).rapidity() - jets->at(1).rapidity()), theWeight); 
    theHistograms.fill(std::string("ZZTo")+decay+"_deltaYJJ_FRVar", "Var From FR #Delta Y(j,j) between the two most energetic jets",  10, 0, 8, fabs(jets->at(0).rapidity() - jets->at(1).rapidity()),ZZ->fakeRateSFVar());
  
}


  if(centralJets->size() >= 2){
    theHistograms.fill(std::string("ZZTo")+decay+"_deltaEtaJJcentral"      , "#Delta #eta(j,j) between the two most energetyc central jets"            ,  10, 0, 8, fabs(centralJets->at(0).eta() - centralJets->at(1).eta()), theWeight); 
    theHistograms.fill(std::string("ZZTo")+decay+"_deltaEtaJJcentral_FRVar", "Var From FR #Delta #eta(j,j) between the two most energetyc central jets",  10, 0, 8, fabs(centralJets->at(0).eta() - centralJets->at(1).eta()), ZZ->fakeRateSFVar());

    theHistograms.fill(std::string("ZZTo")+decay+"_deltaYJJcentral"      , "#Delta Y(j,j) between the two most energetyc central jets"            ,  10, 0, 8, fabs(centralJets->at(0).rapidity() - centralJets->at(1).rapidity()), theWeight); 
    theHistograms.fill(std::string("ZZTo")+decay+"_deltaYJJcentral_FRVar", "Var From FR #Delta Y(j,j) between the two most energetyc central jets",  10, 0, 8, fabs(centralJets->at(0).rapidity() - centralJets->at(1).rapidity()), ZZ->fakeRateSFVar());


    theHistograms.fill(std::string("ZZTo")+decay+"_mJJ"      , "m_{jj}"            ,  20, 0, 1000, (centralJets->at(0).p4() + centralJets->at(1).p4()).M(), theWeight);
    theHistograms.fill(std::string("ZZTo")+decay+"_mJJ_FRVar", "Var From FR m_{jj}",  20, 0, 1000, (centralJets->at(0).p4() + centralJets->at(1).p4()).M(), ZZ->fakeRateSFVar());
 
  }

  theHistograms.fill(std::string("ZZTo")+decay+"_nExtraMuons"    , "Number of extra muons in the event"    , 10, 0, 10, muons->size(), theWeight); 
  theHistograms.fill(std::string("ZZTo")+decay+"_nExtraElectrons", "Number of extra electrons in the event", 10, 0, 10, electrons->size(), theWeight); 
  theHistograms.fill(std::string("ZZTo")+decay+"_nExtraLeptons"  , "Number of extra leptons in the event"  , 10, 0, 10, muons->size()+electrons->size(), theWeight); 

}


void ZZjAnalyzer::analyze(){

  // theHistograms.fill("fakeRateWeight", "FakeRate",  100, -2, 2, ZZ->fakeRateSF() , 1); 
  // if(ZZ->numberOfGoodGrandDaughters() == 3) theHistograms.fill("fakeRateWeight_3p1F", "FakeRate 3P1F",  100, -2, 2, ZZ->fakeRateSF() , 1); 
  // if(ZZ->numberOfGoodGrandDaughters() == 2) theHistograms.fill("fakeRateWeight_2p2F", "FakeRate 2P2F",  100, -2, 2, ZZ->fakeRateSF() , 1); 

  theHistograms.fill("ZZMass","Invariant Mass ",200,55,1000,ZZ->mass(),theWeight);
  theHistograms.fill("ZZMass_FRVar","Var From FR Invariant Mass ",200,55,1000,ZZ->mass(), ZZ->fakeRateSFVar());

  // Some basic plots on ZZ
  ZZplots();   // ZZ --> 4l
  ZZplots(52); // ZZ --> 4m
  ZZplots(48); // ZZ --> 2e2m
  ZZplots(44); // ZZ --> 4e

  if (topology.test(0)) theHistograms.fill("PassDef", "Number of events passing the signal definition", 200, 50,  1000, ZZ->mass(),theWeight);

  else theHistograms.fill("NoPassDef", "Number of events not passing the signal definition", 200, 50,  1000, ZZ->mass(),theWeight);

}



  
