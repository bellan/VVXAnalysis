#include "VVXAnalysis/TreeAnalysis/interface/ZZjAnalyzer.h"
#include "VVXAnalysis/Commons/interface/SignalDefinitions.h"
#include "VVXAnalysis/Commons/interface/Utilities.h"
#include "VVXAnalysis/Commons/interface/Comparators.h"
#include <boost/foreach.hpp>
#define foreach BOOST_FOREACH

using std::cout;
using std::endl;

using namespace phys;

void ZZjAnalyzer::ZZplots(int id){

  if(ZZ->id() != id && id != -1) return; // -1 here means generic 4l final state

  std::string eventstr=std::to_string(run)+":"+std::to_string(lumiBlock)+":"+std::to_string(event);

  std::string decay  = "4l";
  
  if  (id == 52) {
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
  else if ( decay == "4l") {
    events4l.push_back(eventstr);
  }


  theHistograms.fill(std::string("ZZTo")+decay+"_Z0Mass"         , std::string("Z0  mass  of ZZ_{1}#rightarrow ")+decay            , 200, 0,200,ZZ->first().mass(),theWeight);
  theHistograms.fill(std::string("ZZTo")+decay+"_Z1Mass"         , std::string("Z1  mass  of ZZ_{1}#rightarrow ")+decay            , 200, 0,200,ZZ->second().mass(),theWeight);

  theHistograms.fill(std::string("ZZTo")+decay+"_Z0pt"         , std::string("Z0  p_{T}  of ZZ_{1}#rightarrow ")+decay            , 200, 0,200,ZZ->first().pt(),theWeight);
  theHistograms.fill(std::string("ZZTo")+decay+"_Z1pt"         , std::string("Z1  p_{T}  of ZZ_{1}#rightarrow ")+decay            , 200, 0,200,ZZ->second().pt(),theWeight);

  theHistograms.fill(std::string("ZZTo")+decay+"_Mass"         , std::string("Invariant mass of ZZ_{1}#rightarrow ")+decay            ,  200, 50,  1000, ZZ->mass(), theWeight);
  
  theHistograms.fill(std::string("ZZTo")+decay+"_Mass"+"_FRVarHigh", std::string("Var From FR Invariant mass of ZZ_{1}#rightarrow ")+decay,  200, 50,  1000, ZZ->mass(),ZZ->fakeRateSFVarHigh());

  theHistograms.fill(std::string("ZZTo")+decay+"_Mass"+"_FRVarLow", std::string("Var From FR Invariant mass of ZZ_{1}#rightarrow ")+decay,  200, 50,  1000, ZZ->mass(),ZZ->fakeRateSFVarLow());

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



  Float_t   nJets = 0;
  Float_t   nJetsNoJer = 0;


  if(theMCInfo.isMC()){
    //JER

    //Smearing jet pt (Up and Down distributions just for systematic uncertainty estimate, Central for the standard analysis) without requiring it is signal (1D distributions made of ALL reco events)
    CentralJER_jets  = new std::vector<phys::Jet>();
    CentralJER_jetPt = new std::vector<Double_t>();
    
    CentralJER_jets->clear();
    CentralJER_jetPt->clear();
    
    
    foreach(const phys::Jet &jet, *pjets){
      double jetPt = 0;
      jetPt = jet.pt();
      
      //JER correction (applied only on MC reco). Up and down velues used to assess systematic uncertainty on jet energy resolution
      double newJetPtJER =0; 
      double width = 0;

      width = jet.jer_width(phys::Jet::central);
      
      newJetPtJER = JER_PtSmear(jetPt, width);
       
      if(newJetPtJER > 30){
	CentralJER_jets->push_back(jet);
	CentralJER_jetPt->push_back(newJetPtJER);
      }      
    }
    
 
    Int_t nCentralJERjets = CentralJER_jets->size();

    if (nCentralJERjets>3) nCentralJERjets=3;
    
    stable_sort(jets->begin(), jets->end(), PtComparator());
    stable_sort(CentralJER_jetPt->begin(), CentralJER_jetPt->end(), std::greater<float>());
    stable_sort(CentralJER_jets->begin(), CentralJER_jets->end(), PtComparator());
   
    nJetsNoJer = jets->size();
    jets = CentralJER_jets;
    nJets = CentralJER_jets->size();
  }

  else  nJets = jets->size();

  if(nJets>3) nJets=3;

  theHistograms.fill(std::string("ZZTo")+decay+"_nJets"      , "Number of jets (|#eta|<4.7 and p_T > 30 GeV)"            , 4, 0, 4, nJets, theWeight); 
  theHistograms.fill(std::string("ZZTo")+decay+"_nJetsNoJer"      , "Number of jets (|#eta|<4.7 and p_T > 30 GeV)"            , 4, 0, 4, nJetsNoJer, theWeight); 

  theHistograms.fill(std::string("ZZTo")+decay+"_nJets_FRVarHigh", " Var High From FR Number of jets (|#eta|<4.7 and p_T > 30 GeV)", 4, 0, 4, nJets, ZZ->fakeRateSFVarHigh());
  theHistograms.fill(std::string("ZZTo")+decay+"_nJets_FRVarLow", " Var Low From FR Number of jets (|#eta|<4.7 and p_T > 30 GeV)", 4, 0, 4, nJets, ZZ->fakeRateSFVarLow());
  

  if(jets->size()>0){
    theHistograms.fill(std::string("ZZTo")+decay+"_ptJet"      , "#Delta Y(j,j) between the two most energetyc central jets"            ,  100, 20, 400,jets->at(0).pt(), theWeight);     
    theHistograms.fill(std::string("ZZTo")+decay+"_etaJet"      , "#Delta Y(j,j) between the two most energetyc central jets"            ,  100, -5, 5,jets->at(0).eta(), theWeight);     
  }
  

  if(jets->size() >= 2) {
  

    //  Float_t zZ1 = (jets->at(2).eta()-(jets->at(0).eta() + jets->at(1).eta()))/fabs(jets->at(0).eta() - jets ->at(1).eta());
    Float_t zZ1 = ZZ->first().eta()-(jets->at(0).eta() + jets->at(1).eta())/2;
    Float_t zZ2 = ZZ->second().eta()-(jets->at(0).eta() + jets->at(1).eta())/2;


    //    Float_t SumVecPt = (ZZ->first().p4()+ZZ->second().p4()+jets->at(0).p4()+jets->at(1).p4()).Pt();
    //Float_t SumPt = (ZZ->first().pt()+ZZ->second().pt()+jets->at(0).pt()+jets->at(1).pt());

    Float_t PtRatio = ((ZZ->first().p4()+ZZ->second().p4()+jets->at(0).p4()+jets->at(1).p4()).Pt())/(ZZ->first().pt()+ZZ->second().pt()+jets->at(0).pt()+jets->at(1).pt());
    Float_t PtJRatio = ((jets->at(0).p4()+jets->at(1).p4()).Pt())/(jets->at(0).pt()+jets->at(1).pt());


    theHistograms.fill(std::string("ZZTo")+decay+"_ptRatio"      , "The ratio of transverse momentum of the vector sum of Z1, Z2, tj1, tj2 to the sum of pTs"        ,  50, 0, 0.8, PtRatio, theWeight); 
    theHistograms.fill(std::string("ZZTo")+decay+"_ptJRatio"      , "The ratio of transverse momentum of the vector sum of tj1, tj2 to the sum of pTs"        ,  50, 0, 1., PtJRatio, theWeight); 


    theHistograms.fill(std::string("ZZTo")+decay+"_Z1z"      , "Zeppenfeld Variable fo Z1 wrt the two leading jets"            ,  49, -6, 6, zZ1, theWeight); 
    theHistograms.fill(std::string("ZZTo")+decay+"_Z2z"      , "Zeppenfeld Variable fo Z2 wrt the two leading jets"            ,  49, -6, 6, zZ2, theWeight); 

    theHistograms.fill(std::string("ZZTo")+decay+"_deltaEtaJJ"      , "#Delta #eta(j,j) between the two most energetic jets"            ,  10, 0, 8, fabs(jets->at(0).eta() - jets->at(1).eta()), theWeight); 

    theHistograms.fill(std::string("ZZTo")+decay+"_deltaEtaJJ_FRVarHigh", "VarHigh From FR #Delta #eta(j,j) between the two most energetic jets",  10, 0, 8, fabs(jets->at(0).eta() - jets->at(1).eta()),ZZ->fakeRateSFVarHigh());

    theHistograms.fill(std::string("ZZTo")+decay+"_deltaEtaJJ_FRVarLow", "VarLow From FR #Delta #eta(j,j) between the two most energetic jets",  10, 0, 8, fabs(jets->at(0).eta() - jets->at(1).eta()),ZZ->fakeRateSFVarLow());

    theHistograms.fill(std::string("ZZTo")+decay+"_deltaYJJ"           , "#Delta Y(j,j) between the two most energetic jets"            ,  10, 0, 8, fabs(jets->at(0).rapidity() - jets->at(1).rapidity()), theWeight); 

    theHistograms.fill(std::string("ZZTo")+decay+"_deltaYJJ_FRVarHigh"  ,"VarHigh From FR #Delta Y(j,j) between the two most energetic jets",  10, 0, 8, fabs(jets->at(0).rapidity() - jets->at(1).rapidity()),ZZ->fakeRateSFVarHigh());

    theHistograms.fill(std::string("ZZTo")+decay+"_deltaYJJ_FRVarLow", "VarLow From FR #Delta Y(j,j) between the two most energetic jets",  10, 0, 8, fabs(jets->at(0).rapidity() - jets->at(1).rapidity()),ZZ->fakeRateSFVarLow());
  
}


  if(nJets >= 3){
    Float_t z =  (jets->at(2).eta()-(jets->at(0).eta() + jets->at(1).eta()))/fabs(jets->at(0).eta() - jets ->at(1).eta());

    //      std::cout<<<z1<<std::endl;
    theHistograms.fill(std::string("ZZTo")+decay+"_z"      , "z between the two most energetic jets"            ,  21, -10, 10, z, theWeight); 
  }

  theHistograms.fill(std::string("ZZTo")+decay+"_nExtraMuons"    , "Number of extra muons in the event"    , 10, 0, 10, muons->size(), theWeight); 
  theHistograms.fill(std::string("ZZTo")+decay+"_nExtraElectrons", "Number of extra electrons in the event", 10, 0, 10, electrons->size(), theWeight); 
  theHistograms.fill(std::string("ZZTo")+decay+"_nExtraLeptons"  , "Number of extra leptons in the event"  , 10, 0, 10, muons->size()+electrons->size(), theWeight);  

}

void ZZjAnalyzer::analyze(){

  std::string eventstr=std::to_string(run)+":"+std::to_string(lumiBlock)+":"+std::to_string(event);

  // Uncomment to look into a en event details

  // if(ZZ->second().mass()>60) {
  //   //  if(eventstr=="258448:311:374222263") 


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

 
 

  //  if(passSRZZOnShell){
  // if((region_ == phys::SR)&&(ZZ->first().daughterPtr(0)->pt() >10) && (ZZ->first().daughterPtr(1)->pt() >10)&& (ZZ->second().daughterPtr(0)->pt()>10) && (ZZ->second().daughterPtr(1)->pt()>10)) isZZRegion =1;
  // else if((region_ == phys::CR3P1F) && (((ZZ->second().daughterPtr(0)->pt()>10)&&(ZZ->second().daughterPtr(0)->passFullSel())) || ((ZZ->second().daughterPtr(1)->pt()>10)&&(ZZ->second().daughterPtr(1)->passFullSel()) ))) isZZRegion =1;
  // else if(region_ == phys::CR2P2F) isZZRegion=1;

  // To use the MC region
  // bool isIn = 0;
  //if(regionWord.test(25)) isIn=1; //CR3P1F
  //if(regionWord.test(24)) isIn=1; //CR2P2F
  //if(regionWord.test(26)) isIn=1; //SR

  
  //std::cout<<"lep1 "<<ZZ->second().daughterPtr(0)->passFullSel()<<" "<<ZZ->second().daughterPtr(0)->pfCombRelIsoFSRCorr()<<std::endl;
  //std::cout<<"lep2 "<<ZZ->second().daughterPtr(1)->passFullSel()<<" "<<ZZ->second().daughterPtr(1)->pfCombRelIsoFSRCorr()<<std::endl;
  
  Float_t scaleFacErrSq = ZZ->efficiencySFUnc(); 
  theHistograms.fill("SFErr","Invariant Mass ",200,0,5.,scaleFacErrSq,1);

  theHistograms.fill("ZZMass","Invariant Mass ",200,55,1000,ZZ->mass(),theWeight);
  theHistograms.fill("ZZMass_FRVarHigh","VarHigh From FR Invariant Mass ",200,55,1000,ZZ->mass(), ZZ->fakeRateSFVarHigh());
  theHistograms.fill("ZZMass_FRVarLow","VarLow From FR Invariant Mass ",200,55,1000,ZZ->mass(), ZZ->fakeRateSFVarLow());
  
  
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
 
 
 Float_t w_kf = 1.;
 
 if((theMCInfo.fileName()=="ggZZ2e2mu") || (theMCInfo.fileName()=="ggZZ4e") || (theMCInfo.fileName()=="ggZZ4mu") || (theMCInfo.fileName()=="ggH125") )   w_kf = kFactor_ggZZ ; 
 else if((theMCInfo.fileName()=="ZZTo4l") || (theMCInfo.fileName()=="ZZTo4lamcatnlo")) w_kf = kFactor_qqZZM * kFactor_EWKqqZZ ; 

 theWeight*=w_kf;
 
 
    // Some basic plots on ZZ    
    ZZplots();   // ZZ --> 4l
    ZZplots(52); // ZZ --> 4m
    ZZplots(48); // ZZ --> 2e2m
    ZZplots(44); // ZZ --> 4e
    
    if (topology.test(2)) theHistograms.fill("PassDef", "Number of events passing the signal definition", 200, 50,  1000, ZZ->mass(),theWeight); 
    else theHistograms.fill("NoPassDef", "Number of events not passing the signal definition", 200, 50,  1000, ZZ->mass(),theWeight);
    
  }


void ZZjAnalyzer::end( TFile &) {
  
  // std::cout<<"eeee:"<<std::endl;
  //   for (std::vector<std::string>::iterator it = events4e.begin() ; it != events4e.end(); ++it) std::cout<<"    "<<*it<<std::endl;
  //   std::cout<<"eemm:"<<std::endl;
  //   for (std::vector<std::string>::iterator it = events2e2mu.begin() ; it != events2e2mu.end(); ++it) std::cout<<"    "<<*it<<std::endl;
  //   std::cout<<"mmmm:"<<std::endl;
  //   for (std::vector<std::string>::iterator it = events4mu.begin() ; it != events4mu.end(); ++it) std::cout<<"    "<<*it<<std::endl;

  std::cout<<"Final \neeee "<<events4e.size()<<std::endl<<"eemm "<<events2e2mu.size()<<std::endl<<"mmmm "<<events4mu.size()<<"\nsum "<<events4e.size()+events2e2mu.size()+events4mu.size()<<" total "<<events4l.size()<<std::endl;
    
}

  
double ZZjAnalyzer::JER_PtSmear(double pt, double width)
{
  double ptsmear= gRandom->Gaus(pt,width);    
  return ptsmear;
}
