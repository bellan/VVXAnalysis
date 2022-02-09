#include "VVXAnalysis/TreeAnalysis/interface/ZZjAnalyzer.h"
#include "VVXAnalysis/Commons/interface/SignalDefinitions.h"
#include "VVXAnalysis/Commons/interface/Comparators.h"
#include <boost/foreach.hpp>
#include <boost/assign/std/vector.hpp> 
#include <boost/range/algorithm/remove_if.hpp>
#include <functional>   // std::bind
#include <algorithm>
#include <boost/bind.hpp>
#include <iterator>

#include "VVXAnalysis/Commons/interface/StringTools.h"

#define foreach BOOST_FOREACH

using std::cout;
using std::endl;
using namespace boost::assign; // bring 'operator+=()' into scope
using namespace phys;
//using namespace physmath;

void ZZjAnalyzer::ZZplots(int id){

  if(ZZ->id() != id && id != -1) return; // -1 here means generic 4l final state

  std::string eventstr=std::to_string(run)+":"+std::to_string(lumiBlock)+":"+std::to_string(event);
  std::string decay  = "4l";
  std::string channel  = "";
 

  if  (id == 52) {
    decay    = "4m";
    channel  = "mmmm";
    events4mu.push_back(eventstr);
  }
  
  else if (id == 48) {
    decay    = "2e2m";
    channel  = "eemm";
    events2e2mu.push_back(eventstr);
  }
  else if (id == 44) {
    decay    = "4e";
    channel  = "eeee";
    events4e.push_back(eventstr);
  }
  else if ( decay == "4l") {
    events4l.push_back(eventstr);
  }

  std::string eventStr ="";
  
  if(decay!="4l") eventStr = std::to_string(run)+":"+std::to_string(lumiBlock)+":"+std::to_string(event)+":"+channel+":"+strtool::sRound(ZZ->mass())+":"+strtool::sRound(ZZ->first().mass())+":"+strtool::sRound(ZZ->second().mass())+":"+std::to_string(jets->size());
 
  theHistograms->fill(std::string("ZZTo")+decay+"_Z0Mass"         , std::string("Z0  mass  of ZZ_{1}#rightarrow ")+decay            , 200, 0,200,ZZ->first().mass(),theWeight);
  theHistograms->fill(std::string("ZZTo")+decay+"_Z1Mass"         , std::string("Z1  mass  of ZZ_{1}#rightarrow ")+decay            , 200, 0,200,ZZ->second().mass(),theWeight);
  theHistograms->fill(std::string("ZZTo")+decay+"_Z0pt"           , std::string("Z0  p_{T}  of ZZ_{1}#rightarrow ")+decay            , 200, 0,200,ZZ->first().pt(),theWeight);  
  theHistograms->fill(std::string("ZZTo")+decay+"_Z1pt"           , std::string("Z1  p_{T}  of ZZ_{1}#rightarrow ")+decay            , 200, 0,200,ZZ->second().pt(),theWeight);
  theHistograms->fill(std::string("ZZTo")+decay+"_Mass"           , std::string("Invariant mass of ZZ_{1}#rightarrow ")+decay            ,  200, 50,  1000, ZZ->mass(), theWeight);
  theHistograms->fill(std::string("ZZTo")+decay+"_Z0lep0_sip"         , std::string("sip of  Z0 lep0 of ZZ_{1}#rightarrow ")+decay            ,100, 0,4,ZZ->first().daughterPtr(0)->sip(),theWeight);
  theHistograms->fill(std::string("ZZTo")+decay+"_Z0lep1_sip"         , std::string("sip of  Z0 lep1 of ZZ_{1}#rightarrow ")+decay            , 100, 0,4,ZZ->first().daughterPtr(1)->sip(),theWeight);
  theHistograms->fill(std::string("ZZTo")+decay+"_Z1lep0_sip"         , std::string("sip of  Z1 lep0 of ZZ_{1}#rightarrow ")+decay            , 100, 0,4,ZZ->second().daughterPtr(0)->sip(),theWeight);
  theHistograms->fill(std::string("ZZTo")+decay+"_Z1lep1_sip"         , std::string("sip of  Z1 lep1 of ZZ_{1}#rightarrow ")+decay            , 100, 0,4,ZZ->second().daughterPtr(1)->sip(),theWeight);
  theHistograms->fill(std::string("ZZTo")+decay+"_Z0lep0_pt"         , std::string("pt of  Z0 lep0 of ZZ_{1}#rightarrow ")+decay            ,  150, 0,  300,ZZ->first().daughterPtr(0)->pt(),theWeight);
  theHistograms->fill(std::string("ZZTo")+decay+"_Z0lep1_pt"         , std::string("pt of  Z0 lep1 of ZZ_{1}#rightarrow ")+decay            ,  150, 0,  300,ZZ->first().daughterPtr(1)->pt(),theWeight);
  theHistograms->fill(std::string("ZZTo")+decay+"_Z1lep0_pt"         , std::string("pt of  Z1 lep0 of ZZ_{1}#rightarrow ")+decay            ,  150, 0,  300,ZZ->second().daughterPtr(0)->pt(),theWeight);
  theHistograms->fill(std::string("ZZTo")+decay+"_Z1lep1_pt"         , std::string("pt of  Z1 lep1 of ZZ_{1}#rightarrow ")+decay            ,  150, 0,  300,ZZ->second().daughterPtr(1)->pt(),theWeight);
  theHistograms->fill(std::string("ZZTo")+decay+"_Z0lep0_iso"         , std::string("iso of  Z0 lep0 of ZZ_{1}#rightarrow ")+decay            ,  60, 0, 0.6,ZZ->first().daughterPtr(0)->pfCombRelIsoFSRCorr(),theWeight);
  theHistograms->fill(std::string("ZZTo")+decay+"_Z0lep1_iso"         , std::string("iso of  Z0 lep1 of ZZ_{1}#rightarrow ")+decay            ,  60, 0, 0.6,ZZ->first().daughterPtr(1)->pfCombRelIsoFSRCorr(),theWeight);
  theHistograms->fill(std::string("ZZTo")+decay+"_Z1lep0_iso"         , std::string("iso of  Z1 lep0 of ZZ_{1}#rightarrow ")+decay            ,  60, 0, 0.6,ZZ->second().daughterPtr(0)->pfCombRelIsoFSRCorr(),theWeight);
  theHistograms->fill(std::string("ZZTo")+decay+"_Z1lep1_iso"         , std::string("iso of  Z1 lep1 of ZZ_{1}#rightarrow ")+decay            ,  60, 0, 0.6,ZZ->second().daughterPtr(1)->pfCombRelIsoFSRCorr(),theWeight);
  theHistograms->fill(std::string("ZZTo")+decay+"_Met"         , std::string("Invariant mass of ZZ_{1}#rightarrow ")+decay            ,  200, 0,  800, met->pt(),theWeight);
  theHistograms->fill(std::string("ZZTo")+decay+"_PtZZ"         , std::string("p_{T} of ZZ_{1}#rightarrow ")+decay            ,  300, 0,  300, ZZ->pt(),theWeight);

  //JETS

  Float_t   nJets = 0;
  Float_t   nJets_central = 0;

  nJets = jets->size();
  nJets_central = centralJets->size();
    
  Int_t  njets_gen = genJets->size();
 
  if(nJets>3)   nJets=4;
  if(nJets_central>3)   nJets=4;

  theHistograms->fill("ResMat_ZZTo"+decay+"_nJets", "", 5,0,5,5,0,5,nJets,njets_gen, theWeight);   

  Int_t  njets = nJets;
  if(nJets>2)      njets=3;

  theHistograms->fill("ResMat_ZZTo"+decay+"_nJets_4", "", 4,0,4,4,0,4,njets,njets_gen, theWeight);   
  theHistograms->fill(std::string("ZZTo")+decay+"_nJets"      , "Number of jets (|#eta|<4.7 and p_T > 30 GeV)"            , 5, 0, 5, nJets, theWeight); 
  theHistograms->fill(std::string("ZZTo")+decay+"_nJets_central"      , "Number of jets (|#eta|<4.7 and p_T > 30 GeV)"            , 5, 0, 5, nJets_central, theWeight); 

  
  if(nJets>0){
    eventStr+=":";
    eventStr+=strtool::sRound(jets->at(0).pt());
    theHistograms->fill(std::string("ZZTo")+decay+"_PtJet1"      , "#p_{T} of the most energetic jet"            ,  100, 20, 400,jets->at(0).pt(), theWeight);
    theHistograms->fill(std::string("ZZTo")+decay+"_EtaJet1"      , "#eta  of the most energetic jet"            ,  100, -5, 5,jets->at(0).eta(), theWeight); 

    if(abs(jets->at(0).eta()) > 2.4){
    theHistograms->fill(std::string("ZZTo")+decay+"_PtJet1_noCentral"      , "#p_{T} of the most energetic jet"            ,  100, 20, 400,jets->at(0).pt(), theWeight);
    theHistograms->fill(std::string("ZZTo")+decay+"_EtaJet1_noCentral"      , "#eta  of the most energetic jet"            ,  100, -5, 5,jets->at(0).eta(), theWeight); 
    }
  }

  else{
    eventStr+=":-1.00";
  }  

 
  if(jets->size() >= 2) {
    Float_t  mjj =  (jets->at(0).p4() + jets->at(1).p4()).M();
    eventStr+=":";
    eventStr+=strtool::sRound(jets->at(1).pt());
    eventStr+=":";
    eventStr+=strtool::sRound(mjj); 
  }
  else{
    eventStr+=":-1.00:-1.00";
  }
  //  cout<<" theweight "<<theWeight<<" SF "<<ZZ->efficiencySF()<<" PU "<<theMCInfo.puWeight()<<endl;
  eventStr+=":";
  eventStr+=strtool::sRound(ZZ->efficiencySF()*theMCInfo.puWeight(),".4");

   if(nJets>1) {
    theHistograms->fill(std::string("ZZTo")+decay+"_PtJet2"      , "#p_{T} of the second  most energetic jet"            ,  100, 20, 400,jets->at(1).pt(), theWeight);     
    theHistograms->fill(std::string("ZZTo")+decay+"_EtaJet2"      , "#eta  of the second most energetic jet"            ,  100, -5, 5,jets->at(1).eta(), theWeight);     
    Float_t  mjj =  (jets->at(0).p4() + jets->at(1).p4()).M();
    Float_t  deta = fabs(jets->at(0).eta() - jets->at(1).eta());

    //  Float_t zZ1 = (jets->at(2).eta()-(jets->at(0).eta() + jets->at(1).eta()))/fabs(jets->at(0).eta() - jets ->at(1).eta());
    Float_t zZ1 = ZZ->first().eta()-(jets->at(0).eta() + jets->at(1).eta())/2;
    Float_t zZ2 = ZZ->second().eta()-(jets->at(0).eta() + jets->at(1).eta())/2;        
    Float_t PtRatio = ((ZZ->first().p4()+ZZ->second().p4()+jets->at(0).p4()+jets->at(1).p4()).Pt())/(ZZ->first().pt()+ZZ->second().pt()+jets->at(0).pt()+jets->at(1).pt());
    Float_t PtJRatio = ((jets->at(0).p4()+jets->at(1).p4()).Pt())/(jets->at(0).pt()+jets->at(1).pt());
    Float_t Dphi = physmath::deltaPhi(jets->at(0).phi(),jets->at(1).phi());    
  
    theHistograms->fill(std::string("ZZTo")+decay+"_Dphi"     , "Delta phi of the two leading jets",  100, 0, 3.20, Dphi, theWeight); 
    theHistograms->fill("ZZTo"+decay+"_Mjj","Invariant mass of the two leading jet",100,0,900,mjj,theWeight);
    theHistograms->fill("ZZTo"+decay+"_Deta","Delta eta of the two leading jet",30,0,9,deta,theWeight);

    theHistograms->fill(std::string("ZZTo")+decay+"_ptRatio"      , "The ratio of transverse momentum of the vector sum of Z1, Z2, tj1, tj2 to the sum of pTs"        ,  50, 0, 0.8, PtRatio, theWeight); 
    theHistograms->fill(std::string("ZZTo")+decay+"_ptJRatio"      , "The ratio of transverse momentum of the vector sum of tj1, tj2 to the sum of pTs"        ,  50, 0, 1., PtJRatio, theWeight); 
    theHistograms->fill(std::string("ZZTo")+decay+"_Z1z"      , "Zeppenfeld Variable fo Z1 wrt the two leading jets"            ,  49, -6, 6, zZ1, theWeight);
    theHistograms->fill(std::string("ZZTo")+decay+"_Z2z"      , "Zeppenfeld Variable fo Z2 wrt the two leading jets"            ,  49, -6, 6, zZ2, theWeight); 

    theHistograms->fill(std::string("ZZTo")+decay+"_deltaYJJ"           , "#Delta Y(j,j) between the two most energetic jets"            ,  10, 0, 8, fabs(jets->at(0).rapidity() - jets->at(1).rapidity()), theWeight); 

    if(nJets==2) theHistograms->fill("ZZTo"+decay+"_Deta2Jet","Delta eta of the two leading jet",30,0,9,deta,theWeight);
    if(nJets>2)  theHistograms->fill("ZZTo"+decay+"_Deta3Jet","Delta eta of the two leading jet",30,0,9,deta,theWeight);

    if( (abs(jets->at(0).eta()) > 2.4) || (abs(jets->at(1).eta()) > 2.4) )    theHistograms->fill("ZZTo"+decay+"_Deta_1noCentral","Delta eta of the two leading jet",30,0,9,deta,theWeight); 
    if( (abs(jets->at(0).eta()) > 2.4) && (abs(jets->at(1).eta()) > 2.4) )    theHistograms->fill("ZZTo"+decay+"_Deta_noCentral","Delta eta of the two leading jet",30,0,9,deta,theWeight); 
  }
    

  if(nJets_central>1) {
    Float_t  mjj =  (centralJets->at(0).p4() + centralJets->at(1).p4()).M();
    Float_t  deta = fabs(centralJets->at(0).eta() - centralJets->at(1).eta());
    theHistograms->fill("ZZTo"+decay+"_Mjj_Central","Invariant mass of the two central leading jet",100,0,900,mjj,theWeight);
    theHistograms->fill("ZZTo"+decay+"_Deta_Central","Delta eta of the two leading central jet",30,0,6,deta,theWeight);
  }

  if (decay != "4l")  eventsFull.push_back(eventStr);
  if(nJets >= 3){ 
    Float_t z =  (jets->at(2).eta()-(jets->at(0).eta() + jets->at(1).eta()))/2.  ;
    theHistograms->fill(std::string("ZZTo")+decay+"_z"      , "z between the two most energetic jets"            ,  21, -10, 10, z, theWeight); 
  }
  
  theHistograms->fill(std::string("ZZTo")+decay+"_nExtraMuons"    , "Number of extra muons in the event"    , 10, 0, 10, muons->size(), theWeight); 
  theHistograms->fill(std::string("ZZTo")+decay+"_nExtraElectrons", "Number of extra electrons in the event", 10, 0, 10, electrons->size(), theWeight); 
  theHistograms->fill(std::string("ZZTo")+decay+"_nExtraLeptons"  , "Number of extra leptons in the event"  , 10, 0, 10, muons->size()+electrons->size(), theWeight);    
}

void ZZjAnalyzer::analyze(){

  std::string eventstr=std::to_string(run)+":"+std::to_string(lumiBlock)+":"+std::to_string(event);


  // Uncomment to look into a en event details
  
  bool printDet = kFALSE;
  
  // if(eventstr=="274388:1503:2681886984"){  
  //  foreach(const Int_t  &i, *pileUpIds) {
  //   cout<<
  //   pjets->erase(pjets->begin()+i);
  // }
  
  // for (std::vector<phys::Jet>::const_iterator jet = pjets->begin(); jet != pjets->end();) {
  //   //    if(jet->pt()<30 || !(jet->fullPuId1))) // 0 loose, 1 medium, 2 tight 
  //   if(jet->pt()<30) // 0 loose, 1 medium, 2 tight 
  //     jet = pjets->erase(jet); 
  //   else ++jet;
  // } 
  
  
 if(printDet && ((jets->size()==1 && genJets->size()==0) || (jets->size()==1 && genJets->size()==1))){
    if(jets->size()==1 && genJets->size()==0) cout<<"\n\n##FAKE##\n"<<eventstr<<"\n"<<std::endl;
    else                                      cout<<"\n\n##TRUE##\n"<<eventstr<<"\n"<<std::endl;
    
    std::cout<<"Region Word "<<regionWord<<"\n"
	     <<"ZZMass "<<ZZ->mass()<<"\n"
	     <<"Z0) Mass = "<<ZZ->first().mass()<<" pt = "<<ZZ->first().pt()<<"\n"
	     <<"Z1) Msss = "<<ZZ->second().mass()<<" pt = "<<ZZ->second().pt()<<"\n"

             <<"Z0 l0) pt = "<<ZZ->first().daughterPtr(0)->pt()<<" eta = "<<ZZ->first().daughterPtr(0)->eta()<<" sip = "<<ZZ->first().daughterPtr(0)->sip()<<" iso = "<<ZZ->first().daughterPtr(0)->pfCombRelIsoFSRCorr()
	     <<" iso noFSR = "<<ZZ->first().daughterPtr(0)->pfCombRelIso()<<"\n"
	     <<"       is good "<<ZZ->first().daughterPtr(0)->isGood()<<" FSR photon pt = "<<ZZ->first().fsrPhoton(0).pt()<<"\n"     

             <<"Z0 l1) pt = "<<ZZ->first().daughterPtr(1)->pt()<<" eta = "<<ZZ->first().daughterPtr(1)->eta()<<" sip = "<<ZZ->first().daughterPtr(1)->sip()<<" iso = "<<ZZ->first().daughterPtr(1)->pfCombRelIsoFSRCorr()
	     <<" iso noFSR = "<<ZZ->first().daughterPtr(1)->pfCombRelIso()<<"\n"
	     <<"       is good "<<ZZ->first().daughterPtr(1)->isGood()<<" FSR photon pt = "<<ZZ->first().fsrPhoton(0).pt()<<"\n"     

             <<"Z1 l0) pt = "<<ZZ->second().daughterPtr(0)->pt()<<" eta = "<<ZZ->second().daughterPtr(0)->eta()<<" sip = "<<ZZ->second().daughterPtr(0)->sip()<<" iso = "<<ZZ->second().daughterPtr(0)->pfCombRelIsoFSRCorr()
	     <<" iso noFSR = "<<ZZ->second().daughterPtr(0)->pfCombRelIso()<<"\n"
	     <<"       is good "<<ZZ->second().daughterPtr(0)->isGood()<<" FSR photon pt = "<<ZZ->second().fsrPhoton(0).pt()<<"\n"     

             <<"Z1 l1) pt = "<<ZZ->second().daughterPtr(1)->pt()<<" eta = "<<ZZ->second().daughterPtr(1)->eta()<<" sip = "<<ZZ->second().daughterPtr(1)->sip()<<" iso = "<<ZZ->second().daughterPtr(1)->pfCombRelIsoFSRCorr()
	     <<" iso noFSR = "<<ZZ->second().daughterPtr(1)->pfCombRelIso()<<"\n"
	     <<"       is good "<<ZZ->second().daughterPtr(1)->isGood()<<" FSR photon pt = "<<ZZ->second().fsrPhoton(0).pt()<<"\n"     

	     <<"is2p2f "<<regionWord.test(phys::CR2P2F)<<" is3p1f "<<regionWord.test(phys::CR3P1F)<<" is ZZSR "<<regionWord.test(phys::SR4P)
	     <<std::endl; 
    int i = 0;
    foreach(const phys::Jet &jet, *jets){
      cout<<"jet "<<i<<") pt "<<jet.pt()<<"    eta "<<jet.eta()
	  <<"  sum pt^2 = "<<jet.ptd()  
	  <<"  jetArea = "<<jet.jetArea()<<std::endl;
    }  
  }
  
  Float_t scaleFacErrSq = ZZ->efficiencySFUnc(); 
  theHistograms->fill("SFErr","Invariant Mass ",200,0,5.,scaleFacErrSq,1);
  theHistograms->fill("ZZMass","Invariant Mass ",200,55,1000,ZZ->mass(),theWeight);
  
  if((region_ == phys::CR3P1F) && 
     ((((abs(ZZ->second().daughterPtr(0)->id())==11) &&  ZZ->second().daughterPtr(0)->pfCombRelIsoFSRCorr() >0.5) || 
       ((abs(ZZ->second().daughterPtr(0)->id())==13) &&  ZZ->second().daughterPtr(0)->pfCombRelIsoFSRCorr() >0.4) ) && 
      (((abs(ZZ->second().daughterPtr(1)->id())==11) &&  ZZ->second().daughterPtr(1)->pfCombRelIsoFSRCorr() >0.5) || 
       ((abs(ZZ->second().daughterPtr(1)->id())==13) &&  ZZ->second().daughterPtr(1)->pfCombRelIsoFSRCorr() >0.4) )))
    std::cout<<" 2 lep out od iso  ev "<<eventstr<<std::endl; 

  
  if((region_ == phys::CR3P1F) && ( ( !(ZZ->second().daughterPtr(0)->passFullSel()) ) && ( !(ZZ->second().daughterPtr(1)->passFullSel()) )))  std::cout<<" 2 FullSell "<< " ev "<<eventstr<<std::endl;
  if((region_ == phys::CR2P2F) && (  (ZZ->second().daughterPtr(0)->passFullSel())  || (ZZ->second().daughterPtr(1)->passFullSel() )))  std::cout<<" some FullSell "<< " ev "<<eventstr<<std::endl;
  
  
  if((abs(ZZ->first().daughterPtr(0)->id())!=(abs(ZZ->second().daughterPtr(0)->id())))  && ((physmath::deltaR(*ZZ->first().daughterPtr(0), *ZZ->second().daughterPtr(0))<0.05) || 
											    (physmath::deltaR(*ZZ->first().daughterPtr(0), *ZZ->second().daughterPtr(1))<0.05) || 
											    (physmath::deltaR(*ZZ->first().daughterPtr(1), *ZZ->second().daughterPtr(0))<0.05) || 
											    (physmath::deltaR(*ZZ->first().daughterPtr(1), *ZZ->second().daughterPtr(1))<0.05)  )) 
    std::cout<<" delta R  <0.05 "<< " ev "<<eventstr<<std::endl;

if( 
   (physmath::deltaR(*ZZ->first().daughterPtr(0), *ZZ->second().daughterPtr(0))<0.02) || 
   (physmath::deltaR(*ZZ->first().daughterPtr(0), *ZZ->second().daughterPtr(1))<0.02) || 
   (physmath::deltaR(*ZZ->first().daughterPtr(1), *ZZ->second().daughterPtr(0))<0.02) || 
   (physmath::deltaR(*ZZ->first().daughterPtr(1), *ZZ->second().daughterPtr(1))<0.02)  ) 
    std::cout<<" delta R  < 0.02 "<<" ev "<<eventstr<<std::endl;

 if( ZZ->mass() < 70.)  std::cout<<" mass < 70 "<<" ev "<<eventstr<<" ZZ mass "<<ZZ->mass()<<std::endl;
  
 Float_t w_kf = 1.;
 
 // if((theMCInfo.fileName()=="ggZZ2e2mu") || (theMCInfo.fileName()=="ggZZ4e") || (theMCInfo.fileName()=="ggZZ4mu") || (theMCInfo.fileName()=="ggH125") )   w_kf = theMCInfo.kF_ggZZ() ; 
 //else if((theMCInfo.fileName()=="ZZTo4l") || (theMCInfo.fileName()=="ZZTo4lamcatnlo")) w_kf = theMCInfo.kF_qqZZM() * theMCInfo.kF_EWKqqZZ() ; 
 if((theMCInfo.fileName()=="ggZZ2e2mu") || (theMCInfo.fileName()=="ggZZ4e") || (theMCInfo.fileName()=="ggZZ4mu") || (theMCInfo.fileName()=="ggTo2e2mu_Contin_MCFM701") || (theMCInfo.fileName()=="ggTo4e_Contin_MCFM701") || (theMCInfo.fileName()=="ggTo4mu_Contin_MCFM701"))   w_kf = 1.7 ; 
 else if(theMCInfo.fileName()=="ZZTo4l") w_kf = 1.1; 

 theWeight*=w_kf;


 // std::bitset<16> trigger(triggerWord);
 // if( trigger.test(7) && !trigger.test(1) && !trigger.test(2) && !trigger.test(3) && !trigger.test(4) && !trigger.test(4) && !trigger.test(6) ) cout<<"hei ";
 // else{


 // bool  pass12ele = 0;

 // Float_t secMax = 0;
 // Int_t idSecMax= 0;

 // Float_t max = ZZ->first().daughterPtr(0)->pt();
 // Int_t idMax = abs(ZZ->first().daughterPtr(0)->id());
 
 // if(ZZ->first().daughterPtr(1)->pt() > max){ 
 //   secMax = max;
 //   secIdMax = idMax; 
 //   max = ZZ->first().daughterPtr(1)->pt();
 //   idMax = abs(ZZ->first().daughterPtr(1)->id())
 //     }
 

 std::vector<phys::Lepton> Leps;

 Leps.push_back(*ZZ->first().daughterPtr(0));
 Leps.push_back(*ZZ->first().daughterPtr(1));
 Leps.push_back(*ZZ->second().daughterPtr(0));
 Leps.push_back(*ZZ->second().daughterPtr(1));

 stable_sort(Leps.begin(),     Leps.end(),     phys::PtComparator());
 
 
 if((abs(Leps.at(1).id())==11) && (Leps.at(1).pt()<12)) return;


   ZZplots();   // ZZ --> 4l
   ZZplots(48); // ZZ --> 2e2m
   ZZplots(44); // ZZ --> 4e
   ZZplots(52); // ZZ --> 4m
   if (topology.test(2)) theHistograms->fill("PassDef", "Number of events passing the signal definition", 200, 50,  1000, ZZ->mass(),theWeight); 
   else theHistograms->fill("NoPassDef", "Number of events not passing the signal definition", 200, 50,  1000, ZZ->mass(),theWeight);
   if (topology.test(3)) theHistograms->fill("PassFidDef", "Number of events passing the signal definition", 200, 50,  1000, ZZ->mass(),theWeight); 
   else theHistograms->fill("NoPassFidDef", "Number of events not passing the signal definition", 200, 50,  1000, ZZ->mass(),theWeight);
 
}

void ZZjAnalyzer::begin() {

  Xbins_pt += 20,30,50,100;
  Xbins_eta += 0,2.5,2.75,3.0,5.0;
  myfile.open("example.txt",std::ios::app);
}

void ZZjAnalyzer::end( TFile &) {
  
  // std::cout<<"eeee:"<<std::endl;
  //   for (std::vector<std::string>::iterator it = events4e.begin() ; it != events4e.end(); ++it) std::cout<<"    "<<*it<<std::endl;
  //   std::cout<<"eemm:"<<std::endl;
  //   for (std::vector<std::string>::iterator it = events2e2mu.begin() ; it != events2e2mu.end(); ++it) std::cout<<"    "<<*it<<std::endl;
  //   std::cout<<"mmmm:"<<std::endl;
  //   for (std::vector<std::string>::iterator it = events4mu.begin() ; it != events4mu.end(); ++it) std::cout<<"    "<<*it<<std::endl;

  cout<<"size "<<eventsFull.size()<<endl;
  std::sort (eventsFull.begin(), eventsFull.end(), strtool::sortEvents);
  for (std::vector<std::string>::iterator it = eventsFull.begin() ; it != eventsFull.end(); ++it){
    //    std::cout<<*it<<std::endl;
    myfile<<*it<<"\n";
  }
  cout<<myfile.good()<<endl;
  myfile.close();
  cout<<myfile.good()<<endl;
  std::cout<<"Final \neeee "<<events4e.size()<<std::endl<<"eemm "<<events2e2mu.size()<<std::endl<<"mmmm "<<events4mu.size()<<"\nsum "<<events4e.size()+events2e2mu.size()+events4mu.size()<<" total "<<events4l.size()<<std::endl;
    
}

 
