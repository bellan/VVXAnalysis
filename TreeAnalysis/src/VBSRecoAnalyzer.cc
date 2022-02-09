#include "VVXAnalysis/TreeAnalysis/interface/VBSRecoAnalyzer.h"
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

void VBSRecoAnalyzer::analyze(){
  

  //  if((region_ == phys::SR && topology.test(3)) || (region_ == phys::SR_HZZ && topology.test(1))){ 

  UpJER_jets->clear();
  DownJER_jets->clear();
  UpJES_jets->clear();
  DownJES_jets->clear();
  UpJESData_jets->clear();
  DownJESData_jets->clear();

  std::string decay  = "";
  Int_t id = ZZ->id();
  
  if      (id == 52) {decay = "4m";}
  else if (id == 48) {decay = "2e2m";}
  else if (id == 44) {decay = "4e";}

  Float_t w_kf = 1.;
  nJets = jets->size();


  if((theSampleInfo.fileName()=="ggZZ2e2mu") || (theSampleInfo.fileName()=="ggZZ4e") || (theSampleInfo.fileName()=="ggZZ4mu"))   w_kf = 1.7 ; 
  else if(theSampleInfo.fileName()=="ZZTo4l") w_kf = 1.1; 
              
  theWeight*=w_kf;

  Float_t scaleFacErrSq = ZZ->efficiencySFUnc(); 

  FillHistosJets(decay,theWeight,jets,"01");   
  FillHistosJets(decay,theWeight*(1-scaleFacErrSq),jets,"SFErrSqMinus_01");
  FillHistosJets(decay,theWeight*(1+scaleFacErrSq),jets,"SFErrSqPlus_01");

  if(region_ == phys::CR3P1F || region_ == phys::CR2P2F) FillHistosJets(decay,ZZ->fakeRateSFVar(),jets,"FRVar");
   
  //Data  
  if(!theSampleInfo.isMC()){         
    foreach(const phys::Jet &dataJet, *pjets){
      if(!dataJet.fullPuId(1)) continue;
      double dataJetPt = 0;
      dataJetPt = dataJet.pt();
      
      //JES correction: Up and down velues used to assess systematic uncertainty on jet energy resolution
      double newJetPtJESData_up   =0;  
      double newJetPtJESData_down =0;
      
      newJetPtJESData_up   = dataJetPt*(1 + dataJet.jecUncertainty());
      newJetPtJESData_down = dataJetPt*(1 - dataJet.jecUncertainty());
      
      if(newJetPtJESData_up   < dataJetPt) {cout<<"Error "<<endl; abort();}
      if(newJetPtJESData_up   > 30)     UpJESData_jets->push_back(dataJet);       
      if(newJetPtJESData_down > 30)   DownJESData_jets->push_back(dataJet);
      
      stable_sort(UpJESData_jets->begin(), UpJESData_jets->end(), PtComparator());
      stable_sort(DownJESData_jets->begin(), DownJESData_jets->end(), PtComparator());
    }
     FillHistosJets(decay,theWeight,UpJESData_jets,"JESDataUpSmear_01");
     FillHistosJets(decay,theWeight,DownJESData_jets,"JESDataDownSmear_01");
  }

   else{

    foreach(const phys::Jet &jet, *pjets){
      if(!jet.fullPuId(1)) continue;

      if(jet.ptJerUp() > 30.) 	UpJER_jets->push_back(jet); 
      if(jet.ptJerDn() > 30.)	DownJER_jets->push_back(jet);
  
      //JES correction: Up and down velues used to assess systematic uncertainty on jet energy resolution
      double newJetPtJES_up   = 0;  
      double newJetPtJES_down = 0;
           
      newJetPtJES_up   = jet.pt()*(1+jet.jecUncertainty());
      newJetPtJES_down = jet.pt()*(1-jet.jecUncertainty());
      
      if(newJetPtJES_up > 30)  	  UpJES_jets->push_back(jet);
      if(newJetPtJES_down > 30)   DownJES_jets->push_back(jet);
    }
    
    FillHistosJets(decay,theWeight,UpJER_jets,"JERUpSmear_01");
    FillHistosJets(decay,theWeight,DownJER_jets,"JERDownSmear_01");

    FillHistosJets(decay,theWeight,UpJES_jets,"JESUpSmear_01");
    FillHistosJets(decay,theWeight,DownJES_jets,"JESDownSmear_01");

   }
  }
//}
void VBSRecoAnalyzer::begin() {  

  UpJESData_jets           = new std::vector<phys::Jet>();
  DownJESData_jets         = new std::vector<phys::Jet>();
  
  UpJER_jets               = new std::vector<phys::Jet>();
  DownJER_jets             = new std::vector<phys::Jet>();

  UpJES_jets               = new std::vector<phys::Jet>();
  DownJES_jets             = new std::vector<phys::Jet>();

  UpJESData_jets->reserve(2);           
  DownJESData_jets->reserve(2);         
  
  UpJER_jets->reserve(2);               
  DownJER_jets->reserve(2);             

  UpJES_jets->reserve(2);               
  DownJES_jets->reserve(2);             

  Xbins_mass += 100,200,250,300,350,400,500,600,800;  
  nJets=0;
  mjj=0;
  deltaEtajj=0;
}

void VBSRecoAnalyzer::end( TFile &) {
  cout<<"finish"<<endl;
}

void VBSRecoAnalyzer::FillHistosJets(std::string decay,float Wh,std::vector<phys::Jet> *jetsVec,std::string type){
  if(jetsVec->size()>1){  
    mjj        =  (jetsVec->at(0).p4() + jetsVec->at(1).p4()).M();
    deltaEtajj =  fabs(jetsVec->at(0).eta() - jetsVec->at(1).eta());
    //    if(jetsVec->at(0).pt()>100. && jetsVec->at(1).pt()>70. &&  mjj>300 && deltaEtajj>2.4 && met->pt()<60){
    if(mjj>400 && deltaEtajj>2.4){
      
      theHistograms->fill("ZZTo"+decay+"_Mass_"+type,"", Xbins_mass,  ZZ->mass(),Wh);
    }
  }  
  //  theHistograms->fill("ZZTo"+decay+"_nJets_"+type, "", Xbins_nJets,njets, Wh);     
  // if(njets>0){  
  
  //   stable_sort(jetsVec->begin(), jetsVec->end(), PtComparator());
  //   ptJet1 = jetsVec->at(0).pt();
  //   if (ptJet1>=500) ptJet1 = 499;    
  //   theHistograms->fill("ZZTo"+decay+"_PtJet1_"+type," ", Xbins_ptJet1, ptJet1, Wh); 
    
  //   etaJet1 = fabs(jetsVec->at(0).eta());
  //   if (etaJet1>=4.7) etaJet1 = 4.6;
  //   theHistograms->fill("ZZTo"+decay+"_EtaJet1_"+type, "", Xbins_etaJet1, etaJet1, Wh);
  // }
   
  //  if(njets>1){  
  //    deta = fabs(jetsVec->at(0).eta() - jetsVec->at(1).eta());
     
  //    if (deta>=4.7) deta = 4.6;
  //    mjj =  (jetsVec->at(0).p4() + jetsVec->at(1).p4()).M();
  //    if (mjj>=800) mjj = 799;
  //    ptJet2 = jetsVec->at(1).pt();
  //    if (ptJet2>=500) ptJet2 = 499;

  //    dphi = physmath::deltaPhi(jetsVec->at(0).phi(),jetsVec->at(1).phi());
  //    if (dphi>=6) dphi = 5;
     
  //    theHistograms->fill("ZZTo"+decay+"_PtJet2_"+type, "", Xbins_ptJet2, ptJet2, Wh);
     
  //    etaJet2 = fabs(jetsVec->at(1).eta());
  //    if (etaJet2>=4.7) etaJet2 = 4.6;
     
  //    theHistograms->fill("ZZTo"+decay+"_EtaJet2_"+type, "", Xbins_etaJet2, etaJet2, Wh);
  //    theHistograms->fill("ZZTo"+decay+"_Mjj_"+type,"",Xbins_mjj,mjj,Wh); 
  //    theHistograms->fill("ZZTo"+decay+"_Deta_"+type,"",Xbins_deta,deta,Wh); 
  //    theHistograms->fill("ZZTo"+decay+"_Phi_"+type,"",Xbins_dphi,dphi,Wh); 
  //  }
  //  if(njets>2){  
     
  //    ptJet3 = jetsVec->at(2).pt();
  //    if (ptJet3>=500) ptJet3 = 499;
  //    theHistograms->fill("ZZTo"+decay+"_PtJet3_"+type, "", Xbins_ptJet3, ptJet3, Wh);
  //  }

}






// void VBSRecoAnalyzer::ZZplots(int id){

//   if(ZZ->id() != id && id != -1) return; // -1 here means generic 4l final state

//   std::string eventstr=std::to_string(run)+":"+std::to_string(lumiBlock)+":"+std::to_string(event);
//   std::string decay  = "4l";
//   std::string channel  = "";
 

//   if  (id == 52) {
//     decay    = "4m";
//     channel  = "mmmm";
//     events4mu.push_back(eventstr);
//   }
  
//   else if (id == 48) {
//     decay    = "2e2m";
//     channel  = "eemm";
//     events2e2mu.push_back(eventstr);
//   }
//   else if (id == 44) {
//     decay    = "4e";
//     channel  = "eeee";
//     events4e.push_back(eventstr);
//   }
//   else if ( decay == "4l") {
//     events4l.push_back(eventstr);
//   }
  
//   theHistograms->fill(std::string("ZZTo")+decay+"_Z0Mass"         , std::string("Z0  mass  of ZZ_{1}#rightarrow ")+decay            , 200, 0,200,ZZ->first().mass(),theWeight);
//   theHistograms->fill(std::string("ZZTo")+decay+"_Z1Mass"         , std::string("Z1  mass  of ZZ_{1}#rightarrow ")+decay            , 200, 0,200,ZZ->second().mass(),theWeight);
//   theHistograms->fill(std::string("ZZTo")+decay+"_Z0pt"         , std::string("Z0  p_{T}  of ZZ_{1}#rightarrow ")+decay            , 200, 0,200,ZZ->first().pt(),theWeight);  
//   theHistograms->fill(std::string("ZZTo")+decay+"_Z1pt"         , std::string("Z1  p_{T}  of ZZ_{1}#rightarrow ")+decay            , 200, 0,200,ZZ->second().pt(),theWeight);
//   theHistograms->fill(std::string("ZZTo")+decay+"_Mass"         , std::string("Invariant mass of ZZ_{1}#rightarrow ")+decay            ,  200, 50,  1000, ZZ->mass(), theWeight);  
//   theHistograms->fill(std::string("ZZTo")+decay+"_Mass"+"_FRVar", std::string("Var From FR Invariant mass of ZZ_{1}#rightarrow ")+decay,  200, 50,  1000, ZZ->mass(),ZZ->fakeRateSFVar());

//   theHistograms->fill(std::string("ZZTo")+decay+"_Met"         , std::string("Invariant mass of ZZ_{1}#rightarrow ")+decay            ,  200, 0,  800, met->pt(),theWeight);
//   theHistograms->fill(std::string("ZZTo")+decay+"_PtZZ"         , std::string("p_{T} of ZZ_{1}#rightarrow ")+decay            ,  300, 0,  300, ZZ->pt(),theWeight);
  
//   theHistograms->fill(std::string("ZZTo")+decay+"_ptJet1"      , "#Delta Y(j,j) between the two most energetyc central jets"            ,  100, 100, 400,jets->at(0).pt(), theWeight);     
//   theHistograms->fill(std::string("ZZTo")+decay+"_etaJet1"      , "#Delta Y(j,j) between the two most energetyc central jets"            ,  100, -5, 5,jets->at(0).eta(), theWeight);     
  
//   theHistograms->fill(std::string("ZZTo")+decay+"_ptJet2"      , "#Delta Y(j,j) between the two most energetyc central jets"            ,  100, 70, 370,jets->at(1).pt(), theWeight);     
//   theHistograms->fill(std::string("ZZTo")+decay+"_etaJet2"      , "#Delta Y(j,j) between the two most energetyc central jets"            ,  100, -5, 5,jets->at(1).eta(), theWeight);     
  
//   //  Float_t zZ1 = (jets->at(2).eta()-(jets->at(0).eta() + jets->at(1).eta()))/fabs(jets->at(0).eta() - jets ->at(1).eta());
//   Float_t zZ1 = ZZ->first().eta()-(jets->at(0).eta() + jets->at(1).eta())/2;
//   Float_t zZ2 = ZZ->second().eta()-(jets->at(0).eta() + jets->at(1).eta())/2;
  
//   // Float_t SumVecPt = (ZZ->first().p4()+ZZ->second().p4()+jets->at(0).p4()+jets->at(1).p4()).Pt();
//   // Float_t SumPt = (ZZ->first().pt()+ZZ->second().pt()+jets->at(0).pt()+jets->at(1).pt());
  
//   Float_t PtRatio = ((ZZ->first().p4()+ZZ->second().p4()+jets->at(0).p4()+jets->at(1).p4()).Pt())/(ZZ->first().pt()+ZZ->second().pt()+jets->at(0).pt()+jets->at(1).pt());
//     Float_t PtJRatio = ((jets->at(0).p4()+jets->at(1).p4()).Pt())/(jets->at(0).pt()+jets->at(1).pt());
//     Float_t Dphi = physmath::deltaPhi(jets->at(0).phi(),jets->at(1).phi());    
  
  
//     theHistograms->fill(std::string("ZZTo")+decay+"_Dphi"     , "Delta phi of the two leading jets",  100, 0, 3.20, Dphi, theWeight); 
//     theHistograms->fill("ZZTo"+decay+"_Mjj","Invariant mass of the two leading jet",100,300,1600,mjj,theWeight);
//     theHistograms->fill(std::string("ZZTo")+decay+"_ptRatio"      , "The ratio of transverse momentum of the vector sum of Z1, Z2, tj1, tj2 to the sum of pTs"        ,  50, 0, 0.8, PtRatio, theWeight); 
//     theHistograms->fill(std::string("ZZTo")+decay+"_ptJRatio"      , "The ratio of transverse momentum of the vector sum of tj1, tj2 to the sum of pTs"        ,  50, 0, 1., PtJRatio, theWeight); 
//     theHistograms->fill(std::string("ZZTo")+decay+"_Z1z"      , "Zeppenfeld Variable fo Z1 wrt the two leading jets"            ,  49, -6, 6, zZ1, theWeight); 
//     theHistograms->fill(std::string("ZZTo")+decay+"_Z2z"      , "Zeppenfeld Variable fo Z2 wrt the two leading jets"            ,  49, -6, 6, zZ2, theWeight); 
//     theHistograms->fill(std::string("ZZTo")+decay+"_deltaEtaJJ"      , "#Delta #eta(j,j) between the two most energetic jets"            ,  10, 0, 8, fabs(jets->at(0).eta() - jets->at(1).eta()), theWeight); 
//     theHistograms->fill(std::string("ZZTo")+decay+"_deltaEtaJJ_FRVarLow", "VarLow From FR #Delta #eta(j,j) between the two most energetic jets",  10, 0, 8, fabs(jets->at(0).eta() - jets->at(1).eta()),ZZ->fakeRateSFVar());
//     theHistograms->fill(std::string("ZZTo")+decay+"_deltaYJJ"           , "#Delta Y(j,j) between the two most energetic jets"            ,  10, 0, 8, fabs(jets->at(0).rapidity() - jets->at(1).rapidity()), theWeight); 
//     theHistograms->fill(std::string("ZZTo")+decay+"_deltaYJJ_FRVarLow", "VarLow From FR #Delta Y(j,j) between the two most energetic jets",  10, 0, 8, fabs(jets->at(0).rapidity() - jets->at(1).rapidity()),ZZ->fakeRateSFVar());

  
//  if(nJets >= 3){ 
//     Float_t z =  (jets->at(2).eta()-(jets->at(0).eta() + jets->at(1).eta()))/2.  ;
//     theHistograms->fill(std::string("ZZTo")+decay+"_z"      , "z between the two most energetic jets"            ,  21, -10, 10, z, theWeight); 
//   }
  
//   theHistograms->fill(std::string("ZZTo")+decay+"_nExtraMuons"    , "Number of extra muons in the event"    , 10, 0, 10, muons->size(), theWeight); 
//   theHistograms->fill(std::string("ZZTo")+decay+"_nExtraElectrons", "Number of extra electrons in the event", 10, 0, 10, electrons->size(), theWeight); 
//   theHistograms->fill(std::string("ZZTo")+decay+"_nExtraLeptons"  , "Number of extra leptons in the event"  , 10, 0, 10, muons->size()+electrons->size(), theWeight);    
// }



// }


 
