#include "VVXAnalysis/TreeAnalysis/interface/ZZSAnalyzer.h"
#include "VVXAnalysis/Commons/interface/SignalDefinitions.h"
#include "VVXAnalysis/Commons/interface/Comparators.h"
//#include "VVXAnalysis/Commons/interface/Colours.h"
#include <boost/foreach.hpp>
#define foreach BOOST_FOREACH

#include <TVector.h>
//#include "TPair.h"
#include <TChain.h>
#include <TGraph.h>
#include <TVectorD.h>
#include <vector>
#include <utility>
#include <iostream>
#include <algorithm>


using namespace phys;
using namespace std;
//using namespace colour;

Int_t ZZSAnalyzer::cut() {
  
  if(!ZZ->passTrigger())          theHistograms.fill("CutCheck",5,0,5,0);
  if(! 60 < ZZ->first().mass())   theHistograms.fill("CutCheck",5,0,5,1);
  if(! ZZ->first().mass() < 120)  theHistograms.fill("CutCheck",5,0,5,2);
  if(! 60 < ZZ->second().mass())  theHistograms.fill("CutCheck",5,0,5,3);
  if(! ZZ->second().mass() < 120) theHistograms.fill("CutCheck",5,0,5,4);

  if(ZZ->passTrigger())
    return 1;
  
  theHistograms.fill("mZ1_nonPassing", "Invariant mass of Z_{1} for non passing events",  50, 0, 200, ZZ->first().mass() , theWeight); 
  theHistograms.fill("mZ2_nonPassing", "Invariant mass of Z_{2} for non passing events",  50, 0, 200, ZZ->second().mass(), theWeight); 

  return -1;
}

void ZZSAnalyzer::ZZplots(int id){

   if(ZZ->id() != id && id != -1) return; // -1 here means generic 4l final state

  std::string decay  = "4l";
  std::string decay1 = "l";
  std::string decay2 = "l";

  if      (id == 52) {decay = "4m"  ; decay1 = "m"; decay1 = "m";}
  else if (id == 48) {decay = "2e2m"; decay1 = "e"; decay2 = "m";}
  else if (id == 44) {decay = "4e"  ; decay1 = "e"; decay2 = "e";}

  theHistograms.fill(std::string("mZ1To2")+decay1, std::string("Invariant mass of Z_{1}#rightarrow 2")+decay1,  15, 60,  120, ZZ->first().mass() , theWeight); 
  theHistograms.fill(std::string("mZ2To2")+decay2, std::string("Invariant mass of Z_{2}#rightarrow 2")+decay2,  15, 60,  120, ZZ->second().mass(), theWeight); 
  theHistograms.fill(std::string("mZZTo") +decay , std::string("Invariant mass of ZZ#rightarrow ")    +decay ,  40,  0, 1000, ZZ->mass()         , theWeight);
 
 theHistograms.fill(std::string("nJets_ZZTo") +decay      , std::string("Number of jets (|#eta|<4.7 and p_T > 30 GeV) ZZ#rightarrow ") +decay        , 10, 0, 10, jets->size(), theWeight); 
  theHistograms.fill(std::string("nCentralJets_ZZTo") +decay ,std::string("Number of jets (|#eta|<2.5 and p_T > 30 GeV) ZZ#rightarrow ") +decay, 10, 0, 10, centralJets->size(), theWeight);

 theHistograms.fill(std::string("met_ZZTo") +decay      , std::string("Missing transverse energy (|#eta|<4.7 and p_T > 30 GeV) ZZ#rightarrow ") +decay        , 200, 0, 800, met->pt(), theWeight); 

 if(jets->size() >= 2){
   theHistograms.fill(std::string("DeltaYJJ_ZZTo") +decay, std::string("#Delta Y(j,j) between the two most energetic jets ZZ#rightarrow ") +decay,  50, 0, 8, fabs(jets->at(0).rapidity() - jets->at(1).rapidity()), theWeight); 
   theHistograms.fill(std::string("DeltaEtaJJ_ZZTo") +decay, std::string("#Delta #eta(j,j) between the two most energetic jets ZZ#rightarrow ") +decay,  50, 0, 8, fabs(jets->at(0).eta() - jets->at(1).eta()), theWeight); 
  }  

  if(centralJets->size() >= 2){
    theHistograms.fill(std::string("DeltaYJJcentral_ZZTo") +decay, std::string("#Delta Y(j,j) between the two most energetic central jets ZZ#rightarrow ") +decay ,  50, 0, 8, fabs(centralJets->at(0).rapidity() - centralJets->at(1).rapidity()), theWeight);
    theHistograms.fill(std::string("DeltaEtaJJcentral_ZZTo") +decay,  std::string("#Delta #eta(j,j) between the two most energetic central jets ZZ#rightarrow ") +decay,  50, 0, 8, fabs(centralJets->at(0).eta() - centralJets->at(1).eta()), theWeight);  
 		       theHistograms.fill(std::string("mJJ")+decay, "m_{jj}",  100, 0, 3000, (centralJets->at(0).p4() + centralJets->at(1).p4()).M(), theWeight); 
  }

}

int e =0;
int genCat = 0;
int top_0 = 0;
int top_2 = 0;
int top_3 = 0;
int top_4 = 0;
int top_5 = 0;
int sig = 0;
int inDeltaEta_4l = 0;
int inDeltaY_4l = 0;
int inDeltaEta_ZZ = 0;
int inDeltaY_ZZ = 0;

int q = 0;
int g = 0;
int hq = 0;
int p = 0;
int l =0;
int o = 0;

void ZZSAnalyzer::analyze(){
  
  e++;
  
  if(genCategory>0){
    
    // Particle gjet_1;
    // Jet jet_1;
    
    // std::vector< std::pair<const Particle , const Jet > > JetcomparatorVector;
    std::vector< std::pair<const Particle* , const Jet*> > GenRecoJetMatchedVector;
    std::vector< std::pair<const Particle* , const Particle*> > RecoJetGenParticleMatchedVector;

    foreach(const Jet &jet, *jets){
      double bestDeltaR = 1e30;
      const Particle* matchedGenJet = 0; 
      foreach(const Particle &gjet, *pgenJets){
	double myDeltaR=physmath::deltaR(jet, gjet);
	if (myDeltaR<bestDeltaR) {
	  matchedGenJet=&gjet;
	  bestDeltaR=myDeltaR;
	}
      }
      if(bestDeltaR < 0.4) GenRecoJetMatchedVector.push_back(make_pair(matchedGenJet,&jet));
    }


    for(std::vector<std::pair<const Particle*, const Jet*> >::iterator it = GenRecoJetMatchedVector.begin(); it != GenRecoJetMatchedVector.end(); it++){
     
      std::cout << "DeltaR = " << physmath::deltaR(*(it->first), *(it->second)) << std::endl;
 
      double bestDeltaR = 1e30;
      const Particle* matchedGenParticle =0; 
      foreach(const Particle &gparticle, *genParticles){
     	double myDeltaR=physmath::deltaR(*(it->first), gparticle);
     	if (myDeltaR<bestDeltaR) {
     	  matchedGenParticle=&gparticle;
     	  bestDeltaR=myDeltaR;
     	}
      }
     if(bestDeltaR < 0.4) RecoJetGenParticleMatchedVector.push_back(make_pair(matchedGenParticle,it->second));

    }

   
    
    for(std::vector<std::pair<const Particle*, const Particle*> >::iterator j = RecoJetGenParticleMatchedVector.begin(); j != RecoJetGenParticleMatchedVector.end(); j++){
      
      std::cout << "gen particle id = " << j->first->id() << std::endl; 

      if(abs(j->first->id()) == 1 || abs(j->first->id()) == 2 )	q++;
      
      else if(abs(j->first->id()) == 21)  g++;

      else if(abs(j->first->id()) == 3 ||abs(j->first->id()) == 4 || abs(j->first->id()) == 5 || abs(j->first->id()) == 6 ) hq++;
      

      else if(abs(j->first->id()) == 22) p++;
      
      else if(abs(j->first->id()) == 11 || abs(j->first->id()) == 13) l++;

      else o++;
    }
    
    std::cout << "quark jets = " << q  << " gluon jets = " << g << " heavy jets = " << hq << " photon jets = " << p << " lepton jets = " << l << " other jets = " << o << std::endl;

     
      
    if(topology.test(2)) {
      
      top_2++; 
      theHistograms.fill("ngenJets_top2_sigVVS"       , "Number of generated jets (|#eta|<4.7 and p_T > 30 GeV) - VVS signal condictions required"        , 10, 0, 10, genJets->size(), theWeight);
      
    }

    if(topology.test(3)) top_3++; 
    if(topology.test(4)) top_4++; 
    if(topology.test(5)) top_5++; 
    
    if(topology.test(2)==1 && topology.test(3)==0 && topology.test(4)==0 && topology.test(5)==0){
      
      sig ++;   
      
      theHistograms.fill("ngenJets_sigVVS"       , "Number of generated jets (|#eta|<4.7 and p_T > 30 GeV) - VVS signal condictions required"        , 10, 0, 10, genJets->size(), theWeight);
      theHistograms.fill("ngenCentralJets_sigVVS", "Number of generated central jets (|#eta|<2.5 and p_T > 30 GeV) - VVS signal condictions required", 10, 0, 10,  centralGenJets->size(), theWeight);
      
    }

  }


  if(ZZ->passTrigger()){//baseline selection required? FIXME
    if(topology.test(2)==1 && topology.test(3)==0 && topology.test(4)==0 && topology.test(5)==0){
     	if(jets->size() >= 2){

      
	} 
     }
    else{
    }
  }
}

 

  



  

