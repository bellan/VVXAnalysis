#include "VVXAnalysis/TreeAnalysis/interface/ZZSAnalyzer.h"
#include "VVXAnalysis/Commons/interface/SignalDefinitions.h"
#include "VVXAnalysis/Commons/interface/Comparators.h"

#include <boost/foreach.hpp>
#define foreach BOOST_FOREACH

#include <TVector.h>
//#include "TPair.h"
#include <TChain.h>
#include <TGraph.h>
#include <TVectorD.h>
#include <sstream> 
#include <string> 
// #include <vector>
// #include <boost/assign/std/vector.hpp> // for 'operator+=()'
// #include <boost/assert.hpp>
using namespace std;
//using namespace boost::assign; // bring 'operator+=()' into scope
using namespace phys;

// Int_t ZZSAnalyzer::cut() {
  
//   if(!ZZ->passTrigger())          theHistograms.fill("CutCheck",5,0,5,0);
//   if(! 60 < ZZ->first().mass())   theHistograms.fill("CutCheck",5,0,5,1);
//   if(! ZZ->first().mass() < 120)  theHistograms.fill("CutCheck",5,0,5,2);
//   if(! 60 < ZZ->second().mass())  theHistograms.fill("CutCheck",5,0,5,3);
//   if(! ZZ->second().mass() < 120) theHistograms.fill("CutCheck",5,0,5,4);

//   if(ZZ->passTrigger())
//     return 1;
  
//   theHistograms.fill("mZ1_nonPassing", "Invariant mass of Z_{1} for non passing events",  50, 0, 200, ZZ->first().mass() , theWeight); 
//   theHistograms.fill("mZ2_nonPassing", "Invariant mass of Z_{2} for non passing events",  50, 0, 200, ZZ->second().mass(), theWeight); 

//   return -1;
// }

// void ZZSAnalyzer::ZZplots(int id){

//    if(ZZ->id() != id && id != -1) return; // -1 here means generic 4l final state

//   std::string decay  = "4l";
//   std::string decay1 = "l";
//   std::string decay2 = "l";

//   if      (id == 52) {decay = "4m"  ; decay1 = "m"; decay1 = "m";}
//   else if (id == 48) {decay = "2e2m"; decay1 = "e"; decay2 = "m";}
//   else if (id == 44) {decay = "4e"  ; decay1 = "e"; decay2 = "e";}

//   theHistograms.fill(std::string("mZ1To2")+decay1, std::string("Invariant mass of Z_{1}#rightarrow 2")+decay1,  15, 60,  120, ZZ->first().mass() , theWeight); 
//   theHistograms.fill(std::string("mZ2To2")+decay2, std::string("Invariant mass of Z_{2}#rightarrow 2")+decay2,  15, 60,  120, ZZ->second().mass(), theWeight); 
//   theHistograms.fill(std::string("mZZTo") +decay , std::string("Invariant mass of ZZ#rightarrow ")    +decay ,  40,  0, 1000, ZZ->mass()         , theWeight);
 
//  theHistograms.fill(std::string("nJets_ZZTo") +decay      , std::string("Number of jets (|#eta|<4.7 and p_T > 30 GeV) ZZ#rightarrow ") +decay        , 10, 0, 10, jets->size(), theWeight); 
//   theHistograms.fill(std::string("nCentralJets_ZZTo") +decay ,std::string("Number of jets (|#eta|<2.5 and p_T > 30 GeV) ZZ#rightarrow ") +decay, 10, 0, 10, centralJets->size(), theWeight);

//  theHistograms.fill(std::string("met_ZZTo") +decay      , std::string("Missing transverse energy (|#eta|<4.7 and p_T > 30 GeV) ZZ#rightarrow ") +decay        , 200, 0, 800, met->pt(), theWeight); 

//  if(jets->size() >= 2){
//    theHistograms.fill(std::string("DeltaYJJ_ZZTo") +decay, std::string("#Delta Y(j,j) between the two most energetic jets ZZ#rightarrow ") +decay,  50, 0, 8, fabs(jets->at(0).rapidity() - jets->at(1).rapidity()), theWeight); 
//    theHistograms.fill(std::string("DeltaEtaJJ_ZZTo") +decay, std::string("#Delta #eta(j,j) between the two most energetic jets ZZ#rightarrow ") +decay,  50, 0, 8, fabs(jets->at(0).eta() - jets->at(1).eta()), theWeight); 
//   }  

//   if(centralJets->size() >= 2){
//     theHistograms.fill(std::string("DeltaYJJcentral_ZZTo") +decay, std::string("#Delta Y(j,j) between the two most energetic central jets ZZ#rightarrow ") +decay ,  50, 0, 8, fabs(centralJets->at(0).rapidity() - centralJets->at(1).rapidity()), theWeight);
//     theHistograms.fill(std::string("DeltaEtaJJcentral_ZZTo") +decay,  std::string("#Delta #eta(j,j) between the two most energetic central jets ZZ#rightarrow ") +decay,  50, 0, 8, fabs(centralJets->at(0).eta() - centralJets->at(1).eta()), theWeight);  
//  		       theHistograms.fill(std::string("mJJ")+decay, "m_{jj}",  100, 0, 3000, (centralJets->at(0).p4() + centralJets->at(1).p4()).M(), theWeight); 
//   }

// }



void ZZSAnalyzer::analyze(){
  
  e++;
  
  Long64_t nentries = tree()->GetEntries(); 
  
  // cout << nentries <<endl;
  bins = {100.,200.,250.,300.,350.,400.,500.,600.,800.};  
  
  
  m4L_r = ZZ->mass();  
  
  if(topology.test(0)){
    
    m4L_g = sqrt((genVBParticles->at(0).p4()+genVBParticles->at(1).p4())*(genVBParticles->at(0).p4()+genVBParticles->at(1).p4()));
    
    sig++;
    
    theHistograms.fill("ResponseMatrix_m4L_fb_4l_01", "Response Matrix m_{4l}" , 40, 0, 800, 40, 0, 800, m4L_r ,m4L_g , theWeight);
    theHistograms.fill("ResponseMatrix_m4L_vb_4l_01", "Response Matrix m_{4l}" , bins, bins, m4L_r ,m4L_g , theWeight);
    
    if(ZZ->id()==52){
      theHistograms.fill("ResponseMatrix_m4L_fb_4mu_01", "Response Matrix m_{4#mu}" ,  40, 0, 800, 40, 0, 800, m4L_r ,m4L_g , theWeight);
      theHistograms.fill("ResponseMatrix_m4L_vb_4mu_01", "Response Matrix m_{4#mu}" , bins, bins, m4L_r ,m4L_g , theWeight);
    }
    
    if(ZZ->id()==44){
      theHistograms.fill("ResponseMatrix_m4L_fb_4e_01", "Response Matrix m_{4e}" ,  40, 0, 800, 40, 0, 800, m4L_r ,m4L_g , theWeight);
      theHistograms.fill("ResponseMatrix_m4L_vb_4e_01", "Response Matrix m_{4e}" , bins, bins, m4L_r ,m4L_g , theWeight);
    }
    
    if(ZZ->id()==48){
      theHistograms.fill("ResponseMatrix_m4L_fb_2e2mu_01", "Response Matrix m_{2e2#mu}" ,  40, 0, 800, 40, 0, 800, m4L_r ,m4L_g , theWeight);
      theHistograms.fill("ResponseMatrix_m4L_vb_2e2mu_01", "Response Matrix m_{2e2#mu}" , bins, bins, m4L_r ,m4L_g , theWeight);
    }
    
    
    if(e < nentries/2){
      theHistograms.fill("ResponseMatrix_m4L_fb_4l_0", "Response Matrix m_{4l}" , 40, 0, 800, 40, 0, 800, m4L_r ,m4L_g , theWeight);
      theHistograms.fill("ResponseMatrix_m4L_vb_4l_0", "Response Matrix m_{4l}" , bins, bins, m4L_r ,m4L_g , theWeight);

      if(ZZ->id()==52){
      theHistograms.fill("ResponseMatrix_m4L_fb_4mu_0", "Response Matrix m_{4#mu}" ,  40, 0, 800, 40, 0, 800, m4L_r ,m4L_g , theWeight);
      theHistograms.fill("ResponseMatrix_m4L_vb_4mu_0", "Response Matrix m_{4#mu}" , bins, bins, m4L_r ,m4L_g , theWeight);
      }

      if(ZZ->id()==44){
  	theHistograms.fill("ResponseMatrix_m4L_fb_4e_0", "Response Matrix m_{4e}" ,  40, 0, 800, 40, 0, 800, m4L_r ,m4L_g , theWeight);
  	theHistograms.fill("ResponseMatrix_m4L_vb_4e_0", "Response Matrix m_{4e}" , bins, bins, m4L_r ,m4L_g , theWeight);
      }

      if(ZZ->id()==48){
  	theHistograms.fill("ResponseMatrix_m4L_fb_2e2mu_0", "Response Matrix m_{2e2#mu}" ,  40, 0, 800, 40, 0, 800, m4L_r ,m4L_g , theWeight);
  	theHistograms.fill("ResponseMatrix_m4L_vb_2e2mu_0", "Response Matrix m_{2e2#mu}" , bins, bins, m4L_r ,m4L_g , theWeight);
      }
    }
    
    if(e >= nentries/2){
      theHistograms.fill("ResponseMatrix_m4L_fb_4l_1", "Response Matrix m_{4l}" , 40, 0, 800, 40, 0, 800, m4L_r ,m4L_g , theWeight);
      theHistograms.fill("ResponseMatrix_m4L_vb_4l_1", "Response Matrix m_{4l}" , bins, bins, m4L_r ,m4L_g , theWeight);
      
      if(ZZ->id()==52){
  	theHistograms.fill("ResponseMatrix_m4L_fb_4mu_1", "Response Matrix m_{4#mu}" ,  40, 0, 800, 40, 0, 800, m4L_r ,m4L_g , theWeight);
  	theHistograms.fill("ResponseMatrix_m4L_vb_4mu_1", "Response Matrix m_{4#mu}" , bins, bins, m4L_r ,m4L_g , theWeight);
      }
      
      if(ZZ->id()==44){
  	theHistograms.fill("ResponseMatrix_m4L_fb_4e_1", "Response Matrix m_{4e}" ,  40, 0, 800, 40, 0, 800, m4L_r ,m4L_g , theWeight);
  	theHistograms.fill("ResponseMatrix_m4L_vb_4e_1", "Response Matrix m_{4e}" , bins, bins, m4L_r ,m4L_g , theWeight);
      }
      
      if(ZZ->id()==48){
  	theHistograms.fill("ResponseMatrix_m4L_fb_2e2mu_1", "Response Matrix m_{2e2#mu}" ,  40, 0, 800, 40, 0, 800, m4L_r ,m4L_g , theWeight);
  	theHistograms.fill("ResponseMatrix_m4L_vb_2e2mu_1", "Response Matrix m_{2e2#mu}" , bins, bins, m4L_r ,m4L_g , theWeight);
      } 
    }
        
  }
  
  else {
    bkg++;
    
    //std::cout << "sig = " << sig << " bkg = " << bkg << std::endl;
  }
  
  /////////////Histograms with all events recostructed in the signal region///////////////////
  
  theHistograms.fill("m4L_allSR_fb_4l_01"       , "m_{4l} distribution - reco level - All SR", 40, 0, 800, m4L_r, theWeight);
  theHistograms.fill("m4L_allSR_vb_4l_01"       , "m_{4l} distribution - reco level - All SR",bins, m4L_r, theWeight);
  
  if(ZZ->id()==52){
    theHistograms.fill("m4L_allSR_fb_4mu_01"       , "m_{4#mu} distribution - reco level - All SR", 40, 0, 800, m4L_r, theWeight);
    theHistograms.fill("m4L_allSR_vb_4mu_01"       , "m_{4#mu} distribution - reco level - All SR",bins, m4L_r, theWeight);
  }
  
  if(ZZ->id()==44){
    theHistograms.fill("m4L_allSR_fb_4e_01"       , "m_{4e} distribution - reco level - All SR", 40, 0, 800, m4L_r, theWeight);
    theHistograms.fill("m4L_allSR_vb_4e_01"       , "m_{4e} distribution - reco level - All SR",bins, m4L_r, theWeight);
  }
  
  if(ZZ->id()==48){
    theHistograms.fill("m4L_allSR_fb_2e2mu_01"       , "m_{2e2#mu} distribution - reco level - All SR", 40, 0, 800, m4L_r, theWeight);
    theHistograms.fill("m4L_allSR_vb_2e2mu_01"       , "m_{2e2#mu} distribution - reco level - All SR",bins, m4L_r, theWeight);
  }
  
  
  if(e < nentries/2){
    
    theHistograms.fill("m4L_allSR_fb_4l_0"       , "m_{4l} distribution - reco level - All SR", 40, 0, 800, m4L_r, theWeight);
    theHistograms.fill("m4L_allSR_vb_4l_0"       , "m_{4l} distribution - reco level - All SR",bins, m4L_r, theWeight);
    
    if(ZZ->id()==52){
      theHistograms.fill("m4L_allSR_fb_4mu_0"       , "m_{4#mu} distribution - reco level - All SR", 40, 0, 800, m4L_r, theWeight);
      theHistograms.fill("m4L_allSR_vb_4mu_0"       , "m_{4#mu} distribution - reco level - All SR",bins, m4L_r, theWeight);
    }
    
    if(ZZ->id()==44){
      theHistograms.fill("m4L_allSR_fb_4e_0"       , "m_{4e} distribution - reco level - All SR", 40, 0, 800, m4L_r, theWeight);
      theHistograms.fill("m4L_allSR_vb_4e_0"       , "m_{4e} distribution - reco level - All SR",bins, m4L_r, theWeight);
    }
    
    if(ZZ->id()==48){
      theHistograms.fill("m4L_allSR_fb_2e2mu_0"       , "m_{2e2#mu} distribution - reco level - All SR", 40, 0, 800, m4L_r, theWeight);
      theHistograms.fill("m4L_allSR_vb_2e2mu_0"       , "m_{2e2#mu} distribution - reco level - All SR",bins, m4L_r, theWeight);
    }
    
  }
  
  if(e >= nentries/2){

    theHistograms.fill("m4L_allSR_fb_4l_1"       , "m_{4l} distribution - reco level - All SR", 40, 0, 800, m4L_r, theWeight);
    theHistograms.fill("m4L_allSR_vb_4l_1"       , "m_{4l} distribution - reco level - All SR",bins, m4L_r, theWeight);
    
    if(ZZ->id()==52){
      theHistograms.fill("m4L_allSR_fb_4mu_1"       , "m_{4#mu} distribution - reco level - All SR", 40, 0, 800, m4L_r, theWeight);
      theHistograms.fill("m4L_allSR_vb_4mu_1"       , "m_{4#mu} distribution - reco level - All SR",bins, m4L_r, theWeight);
    }
    
    if(ZZ->id()==44){
      theHistograms.fill("m4L_allSR_fb_4e_1"       , "m_{4e} distribution - reco level - All SR", 40, 0, 800, m4L_r, theWeight);
      theHistograms.fill("m4L_allSR_vb_4e_1"       , "m_{4e} distribution - reco level - All SR",bins, m4L_r, theWeight);
    }
    
    if(ZZ->id()==48){
      theHistograms.fill("m4L_allSR_fb_2e2mu_1"       , "m_{2e2#mu} distribution - reco level - All SR", 40, 0, 800, m4L_r, theWeight);
      theHistograms.fill("m4L_allSR_vb_2e2mu_1"       , "m_{2e2#mu} distribution - reco level - All SR",bins, m4L_r, theWeight);
    }
  }
}

void ZZSAnalyzer::end( TFile &) {
  cout <<nentries << endl;
}



  



  

