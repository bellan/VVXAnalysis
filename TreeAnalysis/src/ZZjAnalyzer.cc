#include "VVXAnalysis/TreeAnalysis/interface/ZZjAnalyzer.h"
#include "VVXAnalysis/Commons/interface/SignalDefinitions.h"
#include "VVXAnalysis/Commons/interface/Utilities.h"

#include <boost/foreach.hpp>
#define foreach BOOST_FOREACH

using std::cout;
using std::endl;

using namespace phys;

void ZZjAnalyzer::ZZplots(int id){

  if(ZZ->id() != id && id != -1) return; // -1 here means generic 4l final state

  std::string eventstr=std::to_string(run)+":"+std::to_string(lumiBlock)+":"+std::to_string(event);

  std::string decay  = "4l";
  
  if      (id == 52) {
    decay = "4m";
    events4mu.push_back(eventstr);
  }
  
  else if (id == 48) {
    decay = "2e2m";
    events2e2mu.push_back(eventstr);
  }
  else if (id == 44) {
    decay = "4e";
    events4e.push_back(eventstr);
}


  //theWeight = theMCInfo.weight(); //To have plot without scale factors and fake rate weight"

 theHistograms.fill(std::string("ZZTo")+decay+"_Z0lep0_sip"         , std::string("sip of  Z0 lep0 of ZZ_{1}#rightarrow ")+decay            , 200, 0,5,ZZ->first().daughterPtr(0)->sip(),theWeight);

 theHistograms.fill(std::string("ZZTo")+decay+"_Z0lep1_sip"         , std::string("sip of  Z0 lep1 of ZZ_{1}#rightarrow ")+decay            , 200, 0,5,ZZ->first().daughterPtr(1)->sip(),theWeight);


  theHistograms.fill(std::string("ZZTo")+decay+"_Z1lep0_sip"         , std::string("sip of  Z1 lep0 of ZZ_{1}#rightarrow ")+decay            , 200, 0,5,ZZ->second().daughterPtr(0)->sip(),theWeight);

  theHistograms.fill(std::string("ZZTo")+decay+"_Z1lep1_sip"         , std::string("sip of  Z1 lep1 of ZZ_{1}#rightarrow ")+decay            , 200, 0,5,ZZ->second().daughterPtr(1)->sip(),theWeight);



  theHistograms.fill(std::string("ZZTo")+decay+"_Z0lep0_pt"         , std::string("pt of  Z0 lep0 of ZZ_{1}#rightarrow ")+decay            ,  300, 0,  300,ZZ->first().daughterPtr(0)->pt(),theWeight);

 theHistograms.fill(std::string("ZZTo")+decay+"_Z0lep1_pt"         , std::string("pt of  Z0 lep1 of ZZ_{1}#rightarrow ")+decay            ,  300, 0,  300,ZZ->first().daughterPtr(1)->pt(),theWeight);


  theHistograms.fill(std::string("ZZTo")+decay+"_Z1lep0_pt"         , std::string("pt of  Z1 lep0 of ZZ_{1}#rightarrow ")+decay            ,  300, 0,  300,ZZ->second().daughterPtr(0)->pt(),theWeight);

  theHistograms.fill(std::string("ZZTo")+decay+"_Z1lep1_pt"         , std::string("pt of  Z1 lep1 of ZZ_{1}#rightarrow ")+decay            ,  300, 0,  300,ZZ->second().daughterPtr(1)->pt(),theWeight);




  theHistograms.fill(std::string("ZZTo")+decay+"_Z0lep0_iso"         , std::string("iso of  Z0 lep0 of ZZ_{1}#rightarrow ")+decay            ,  100, 0, 2,ZZ->first().daughterPtr(0)->pfCombRelIsoFSRCorr(),theWeight);

 theHistograms.fill(std::string("ZZTo")+decay+"_Z0lep1_iso"         , std::string("iso of  Z0 lep1 of ZZ_{1}#rightarrow ")+decay            ,  100, 0, 2,ZZ->first().daughterPtr(1)->pfCombRelIsoFSRCorr(),theWeight);


  theHistograms.fill(std::string("ZZTo")+decay+"_Z1lep0_iso"         , std::string("iso of  Z1 lep0 of ZZ_{1}#rightarrow ")+decay            ,  100, 0, 2,ZZ->second().daughterPtr(0)->pfCombRelIsoFSRCorr(),theWeight);

  theHistograms.fill(std::string("ZZTo")+decay+"_Z1lep1_iso"         , std::string("iso of  Z1 lep1 of ZZ_{1}#rightarrow ")+decay            ,  100, 0, 2,ZZ->second().daughterPtr(1)->pfCombRelIsoFSRCorr(),theWeight);

 

  theHistograms.fill(std::string("ZZTo")+decay+"_Met"         , std::string("Invariant mass of ZZ_{1}#rightarrow ")+decay            ,  200, 0,  800, met->pt(),theWeight);

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


  std::string eventstr=std::to_string(run)+":"+std::to_string(lumiBlock)+":"+std::to_string(event);

  // Uncomment to look into a en event details

  // if(eventstr=="258448:311:374222263") 
  //   std::cout<<"\n"<<eventstr<<" ZZMass "<<ZZ->mass()<<" Z0m "<<ZZ->first().mass()<<" Z1m "<<ZZ->second().mass()<<" id Z0 "<<ZZ->first().daughterPtr(0)->id()<<" id Z1 "<<ZZ->second().daughterPtr(0)->id()<<"\n"

  // 	     <<" pt z0l0 "<<ZZ->first().daughterPtr(0)->pt()<<" pt z0l1 "<<ZZ->first().daughterPtr(1)->pt()<<" pt z1l0 "<<ZZ->second().daughterPtr(0)->pt()<<" pt z1l1 "<<ZZ->second().daughterPtr(1)->pt()<<"\n"
  // 	     <<" eta z0l0 "<<ZZ->first().daughterPtr(0)->eta()<<" eta z0l1 "<<ZZ->first().daughterPtr(1)->eta()<<" eta z1l0 "<<ZZ->second().daughterPtr(0)->eta()<<" eta z1l1 "<<ZZ->second().daughterPtr(1)->eta()<<"\n"
  // 	  <<" sip z0l0 "<<ZZ->first().daughterPtr(0)->sip()<<" sip z0l1 "<<ZZ->first().daughterPtr(1)->sip()<<" sip z1l0 "<<ZZ->second().daughterPtr(0)->sip()<<" sip z1l1 "<<ZZ->second().daughterPtr(1)->sip()<<"\n"
  // 	  <<" iso z0l0 "<<ZZ->first().daughterPtr(0)->pfCombRelIsoFSRCorr()<<" iso z0l1 "<<ZZ->first().daughterPtr(1)->pfCombRelIsoFSRCorr()<<" iso z1l0 "<<ZZ->second().daughterPtr(0)->pfCombRelIsoFSRCorr()<<" iso z1l1 "<<ZZ->second().daughterPtr(1)->pfCombRelIsoFSRCorr()<<"\n"
  // 	  <<" iso noFSR z0l0 "<<ZZ->first().daughterPtr(0)->pfCombRelIso()<<" iso noFSR z0l1 "<<ZZ->first().daughterPtr(1)->pfCombRelIso()<<" iso noFSR z1l0 "<<ZZ->second().daughterPtr(0)->pfCombRelIso()<<" iso noFSR z1l1 "<<ZZ->second().daughterPtr(1)->pfCombRelIso()<<"\n"
  // 	     <<" isgood z0l0 "<<ZZ->first().daughterPtr(0)->isGood()<<" isgood z0l1 "<<ZZ->first().daughterPtr(1)->isGood()<<" isgood z1l0 "<<ZZ->second().daughterPtr(0)->isGood()<<" isgood z1l1 "<<ZZ->second().daughterPtr(1)->isGood()<<"\n"
  // 	  <<"FSR\n"<<"fsr pt z0l0 "<<ZZ->first().fsrPhoton(0).pt()<<" fsr pt z0l1 "<<ZZ->first().fsrPhoton(1).pt()<<" fsr pt z1l0 "<<ZZ->second().fsrPhoton(0).pt()<<" fsr pt z1l1 "<<ZZ->second().fsrPhoton(1).pt()<<"\n"
  // 	     <<"Region Word "<<regionWord<<" is2p2f "<<regionWord.test(24)<<" is3p1f "<<regionWord.test(25)<<" is ZZSR "<<regionWord.test(26)
  // 	     <<std::endl;
  // }


  bool isZZRegion = 0;
 
 

  //  if(passSRZZOnShell){
  // if((region_ == phys::SR)&&(ZZ->first().daughterPtr(0)->pt() >10) && (ZZ->first().daughterPtr(1)->pt() >10)&& (ZZ->second().daughterPtr(0)->pt()>10) && (ZZ->second().daughterPtr(1)->pt()>10)) isZZRegion =1;
  // else if((region_ == phys::CR3P1F) && (((ZZ->second().daughterPtr(0)->pt()>10)&&(ZZ->second().daughterPtr(0)->passFullSel())) || ((ZZ->second().daughterPtr(1)->pt()>10)&&(ZZ->second().daughterPtr(1)->passFullSel()) ))) isZZRegion =1;
  // else if(region_ == phys::CR2P2F) isZZRegion=1;


  // To use the MC region
  // bool isIn = 0;
  //if(regionWord.test(25)) isIn=1; //CR3P1F
  //if(regionWord.test(24)) isIn=1; //CR2P2F
  //if(regionWord.test(26)) isIn=1; //SR

  // if(isZZRegion & isIn){


  // new Pt cuts
  if ((ZZ->first().daughterPtr(0)->pt() >10) && (ZZ->first().daughterPtr(1)->pt() >10)&& (ZZ->second().daughterPtr(0)->pt()>10) && (ZZ->second().daughterPtr(1)->pt()>10) && ZZ->mass()>100.) isZZRegion =1;


  if(isZZRegion) {
 
    //std::cout<<"lep1 "<<ZZ->second().daughterPtr(0)->passFullSel()<<" "<<ZZ->second().daughterPtr(0)->pfCombRelIsoFSRCorr()<<std::endl;
    //std::cout<<"lep2 "<<ZZ->second().daughterPtr(1)->passFullSel()<<" "<<ZZ->second().daughterPtr(1)->pfCombRelIsoFSRCorr()<<std::endl;
    
    theHistograms.fill("ZZMass","Invariant Mass ",200,55,1000,ZZ->mass(),theWeight);
    theHistograms.fill("ZZMass_FRVar","Var From FR Invariant Mass ",200,55,1000,ZZ->mass(), ZZ->fakeRateSFVar());
    

    // std::cout<<"ev "<<eventstr<<std::endl;
    
  if((region_ == phys::CR3P1F) && 
     ((((abs(ZZ->second().daughterPtr(0)->id())==11) &&  ZZ->second().daughterPtr(0)->pfCombRelIsoFSRCorr() >0.5) || 
       ((abs(ZZ->second().daughterPtr(0)->id())==13) &&  ZZ->second().daughterPtr(0)->pfCombRelIsoFSRCorr() >0.4) ) && 
      (((abs(ZZ->second().daughterPtr(1)->id())==11) &&  ZZ->second().daughterPtr(1)->pfCombRelIsoFSRCorr() >0.5) || 
       ((abs(ZZ->second().daughterPtr(1)->id())==13) &&  ZZ->second().daughterPtr(1)->pfCombRelIsoFSRCorr() >0.4) )))
     std::cout<<" 2 lep out od iso  ev "<<eventstr<<std::endl; 

  
  if((region_ == phys::CR3P1F) && ( ( !(ZZ->second().daughterPtr(0)->passFullSel()) ) && ( !(ZZ->second().daughterPtr(1)->passFullSel()) )))  std::cout<<" 2 FullSell "<< " ev "<<eventstr<<std::endl;
  if((region_ == phys::CR2P2F) && (  (ZZ->second().daughterPtr(0)->passFullSel())  || (ZZ->second().daughterPtr(1)->passFullSel() )))  std::cout<<" some FullSell "<< " ev "<<eventstr<<std::endl;
  
  
  if((abs(ZZ->first().daughterPtr(0)->id())!=(abs(ZZ->second().daughterPtr(0)->id())))  && ( 
 (deltaR(ZZ->first().daughterPtr(0)->p4().Rapidity(), ZZ->first().daughterPtr(0)->p4().Phi(), ZZ->second().daughterPtr(0)->p4().Rapidity(), ZZ->second().daughterPtr(0)->p4().Phi())<0.05) || 
 (deltaR(ZZ->first().daughterPtr(0)->p4().Rapidity(), ZZ->first().daughterPtr(0)->p4().Phi(), ZZ->second().daughterPtr(1)->p4().Rapidity(), ZZ->second().daughterPtr(1)->p4().Phi())<0.05) || 
 (deltaR(ZZ->first().daughterPtr(1)->p4().Rapidity(), ZZ->first().daughterPtr(1)->p4().Phi(), ZZ->second().daughterPtr(0)->p4().Rapidity(), ZZ->second().daughterPtr(0)->p4().Phi())<0.05) || 
 (deltaR(ZZ->first().daughterPtr(1)->p4().Rapidity(), ZZ->first().daughterPtr(1)->p4().Phi(), ZZ->second().daughterPtr(1)->p4().Rapidity(), ZZ->second().daughterPtr(1)->p4().Phi())<0.05)  )) 
    std::cout<<" delta R  <0.05 "<< " ev "<<eventstr<<std::endl;

if( 
 (deltaR(ZZ->first().daughterPtr(0)->p4().Rapidity(), ZZ->first().daughterPtr(0)->p4().Phi(), ZZ->second().daughterPtr(0)->p4().Rapidity(), ZZ->second().daughterPtr(0)->p4().Phi())<0.02) || 
 (deltaR(ZZ->first().daughterPtr(0)->p4().Rapidity(), ZZ->first().daughterPtr(0)->p4().Phi(), ZZ->second().daughterPtr(1)->p4().Rapidity(), ZZ->second().daughterPtr(1)->p4().Phi())<0.02) || 
 (deltaR(ZZ->first().daughterPtr(1)->p4().Rapidity(), ZZ->first().daughterPtr(1)->p4().Phi(), ZZ->second().daughterPtr(0)->p4().Rapidity(), ZZ->second().daughterPtr(0)->p4().Phi())<0.02) || 
 (deltaR(ZZ->first().daughterPtr(1)->p4().Rapidity(), ZZ->first().daughterPtr(1)->p4().Phi(), ZZ->second().daughterPtr(1)->p4().Rapidity(), ZZ->second().daughterPtr(1)->p4().Phi())<0.02)  ) 
    std::cout<<" delta R  < 0.02 "<<" ev "<<eventstr<<std::endl;

 if( ZZ->mass() < 100.)  std::cout<<" mass < 100 "<<" ev "<<eventstr<<" ZZ mass "<<ZZ->mass()<<std::endl;


    // Some basic plots on ZZ    
    ZZplots();   // ZZ --> 4l
    ZZplots(52); // ZZ --> 4m
    ZZplots(48); // ZZ --> 2e2m
    ZZplots(44); // ZZ --> 4e
    
    if (topology.test(0)) theHistograms.fill("PassDef", "Number of events passing the signal definition", 200, 50,  1000, ZZ->mass(),theWeight); 
    else theHistograms.fill("NoPassDef", "Number of events not passing the signal definition", 200, 50,  1000, ZZ->mass(),theWeight);
    
  }
}

void ZZjAnalyzer::end( TFile &) {
  
  // std::cout<<"eeee:"<<std::endl;
  //   for (std::vector<std::string>::iterator it = events4e.begin() ; it != events4e.end(); ++it) std::cout<<"    "<<*it<<std::endl;
  //   std::cout<<"eemm:"<<std::endl;
  //   for (std::vector<std::string>::iterator it = events2e2mu.begin() ; it != events2e2mu.end(); ++it) std::cout<<"    "<<*it<<std::endl;
  //   std::cout<<"mmmm:"<<std::endl;
  //   for (std::vector<std::string>::iterator it = events4mu.begin() ; it != events4mu.end(); ++it) std::cout<<"    "<<*it<<std::endl;

    std::cout<<"Final \neeee "<<events4e.size()<<std::endl<<"eemm "<<events2e2mu.size()<<std::endl<<"mmmm "<<events4mu.size()<<std::endl;
    
}

  
