#include "VVXAnalysis/TreeAnalysis/interface/ZZRecoAnalyzer.h"
#include "VVXAnalysis/Commons/interface/SignalDefinitions.h"
#include "VVXAnalysis/Commons/interface/Utilities.h"

#include "TRandom.h"
#include "TTree.h"
#include <boost/assign/std/vector.hpp> 
#include <boost/assert.hpp> 
using namespace std;
using namespace boost::assign; // bring 'operator+=()' into scope

#include <boost/foreach.hpp>
#define foreach BOOST_FOREACH

using std::cout;
using std::endl;

using namespace phys;

void ZZRecoAnalyzer::ZZplots(int id, int e){

  if(ZZ->id() != id && id != -1) return; // -1 here means generic 4l final state
  
  std::string decay  = "4l";
  std::string sample;
  if      (id == 52) {decay = "4m";}
  else if (id == 48) {decay = "2e2m";}
  else if (id == 44) {decay = "4e";}
  
  if(e < nentries/2){sample = "0";}
  else {sample = "1";}
  
  Int_t njets = jets->size();
  if (njets>3) njets=3;
  theHistograms.fill(std::string("ZZTo")+decay+"_Jets_"+sample, std::string("Number of jets of ZZ_{1}#rightarrow ")+decay,4,0,4,njets,theWeight); 
  theHistograms.fill(std::string("ZZTo")+decay+"_Jets_01", std::string("Number of jets of ZZ_{1}#rightarrow ")+decay,4,0,4,njets,theWeight);  
  theHistograms.fill(std::string("ZZTo")+decay+"_Mass_01", std::string("Invariant mass of ZZ_{1}#rightarrow ")+decay , Xbins, ZZ->mass(),theWeight);
  theHistograms.fill(std::string("ZZTo")+decay+"_Mass_"+sample, std::string("Invariant mass of ZZ_{1}#rightarrow ")+decay , Xbins, ZZ->mass(),theWeight);
  
  if (genCategory !=-1 && topology.test(0)){
    Int_t ngenjets =  genJets->size(); 
    if (ngenjets>3) ngenjets=3;
    theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_Jets_"+sample, std::string("Response matrix number of jets of ZZ_{1}#rightarrow ")+decay, 4,0,4,4,0,4, njets,ngenjets, theWeight); 
    theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_Jets_01", std::string("Response matrix number of jets of ZZ_{1}#rightarrow ")+decay, 4,0,4,4,0,4, njets,ngenjets, theWeight); 
    
    m4L_gen = sqrt((genVBParticles->at(0).p4()+genVBParticles->at(1).p4())*(genVBParticles->at(0).p4()+genVBParticles->at(1).p4()));
    theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_Mass_"+sample, std::string("Response matrix invariant mass of ZZ_{1}#rightarrow ")+decay, Xbins, Xbins, ZZ->mass() ,m4L_gen , theWeight); 
    theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_Mass_01", std::string("Response matrix invariant mass of ZZ_{1}#rightarrow ")+decay, Xbins, Xbins, ZZ->mass() ,m4L_gen , theWeight);
  }
  
  if(region_ == phys::CR3P1F || region_ == phys::CR2P2F) {
    theHistograms.fill(std::string("ZZTo")+decay+"_Mass"+"_FRVar", std::string("Var From FR Invariant mass of ZZ_{1}#rightarrow ")+decay, Xbins, ZZ->mass(),ZZ->fakeRateSFVar());
    theHistograms.fill(std::string("ZZTo")+decay+"_Jets"+"_FRVar", std::string("Var From FR Invariant mass of ZZ_{1}#rightarrow ")+decay ,4,0,4,njets,ZZ->fakeRateSFVar());
  }
}

double ZZRecoAnalyzer::JER_PtSmear(double pt, double width)
{
  double ptsmear= gRandom->Gaus(pt,width);    
  return ptsmear;
}

void ZZRecoAnalyzer::analyze(){
  
  e++;
  
  ZZplots();   // ZZ --> 4l
  ZZplots(52,e); // ZZ --> 4m
  ZZplots(48,e); // ZZ --> 2e2m
  ZZplots(44,e); // ZZ --> 4e
  
  if (genCategory !=-1 && topology.test(0)){     
    centralJER_jets  = new std::vector<phys::Jet>();
    centralJER_jets_up  = new std::vector<phys::Jet>();
    centralJER_jets_down  = new std::vector<phys::Jet>();
    
    centralJER_jets->clear();
    centralJER_jets_up->clear();
    centralJER_jets_down->clear();
    
    foreach(const phys::Jet &jet, *pjets){
      double jetPt = 0;
      double newJetPt =0; 
      double newJetPt_up =0;  
      double newJetPt_down =0;
      double width = 0;
      double width_up = 0; 
      double width_down = 0;
      
      width = jet.jer_width(phys::Jet::central);
      width_up = jet.jer_width(phys::Jet::up); 
      width_down = jet.jer_width(phys::Jet::down);
      
      jetPt = jet.pt();
      newJetPt = JER_PtSmear(jetPt, width);
      newJetPt_up = JER_PtSmear(jetPt, width_up);  
      newJetPt_down = JER_PtSmear(jetPt, width_down);
      
      if(newJetPt > 30) centralJER_jets->push_back(jet);
      if(newJetPt_up > 30) centralJER_jets_up->push_back(jet); 
      if(newJetPt_down > 30) centralJER_jets_down->push_back(jet);
      
    }
    
    
    std::string decay  = "4l";
    std::string sample;
    if      (ZZ->id() == 52) {decay = "4m";}
    else if (ZZ->id() == 48) {decay = "2e2m";}
    else if (ZZ->id() == 44) {decay = "4e";}

    if(e < nentries/2){sample = "0";}
    else {sample = "1";}
    
    theHistograms.fill(std::string("ZZTo")+decay+"_Jets_JERCentralSmear_"+sample, "Number of reco jets (|#eta|<4.7 and p_T > 30", 4,0,4,centralJER_jets->size() , theWeight);
    theHistograms.fill(std::string("ZZTo")+decay+"_Jets_JERUpSmear_"+sample, "Number of reco jets (|#eta|<4.7 and p_T > 30", 4,0,4,centralJER_jets_up->size() , theWeight);
    theHistograms.fill(std::string("ZZTo")+decay+"_Jets_JERDownSmear_"+sample, "Number of reco jets (|#eta|<4.7 and p_T > 30", 4,0,4,centralJER_jets_down->size() , theWeight);
    theHistograms.fill(std::string("ZZTo")+decay+"_Jets_JERCentralSmear_01", "Number of reco jets (|#eta|<4.7 and p_T > 30", 4,0,4,centralJER_jets->size() , theWeight);
    theHistograms.fill(std::string("ZZTo")+decay+"_Jets_JERUpSmear_01", "Number of reco jets (|#eta|<4.7 and p_T > 30", 4,0,4,centralJER_jets_up->size() , theWeight);
    theHistograms.fill(std::string("ZZTo")+decay+"_Jets_JERDownSmear_01", "Number of reco jets (|#eta|<4.7 and p_T > 30", 4,0,4,centralJER_jets_down->size() , theWeight);
  }
}

void ZZRecoAnalyzer::begin() {
  nentries = tree()->GetEntries("ZZCand.passFullSel_");
  Xbins += 100,200,250,300,350,400,500,600,800,1000;
  m4L_gen =0;
  
}

void ZZRecoAnalyzer::end( TFile &) {
  cout<<theMCInfo.analyzedEvents()<< " " << e <<endl;
  
  if(region_ == phys::CR3P1F || region_ == phys::CR2P2F) {
    vector<std::string>  FinalState = {"4m","4e","2e2m"};
    
    for (std::vector<std::string>::iterator it = FinalState.begin() ; it != FinalState.end(); ++it){
      
      TH1 *hvar =  new TH1F();
      hvar =  theHistograms.get(("ZZTo"+*it+"_Mass_FRVar").c_str());
      
      TH1 *h =  new TH1F();
      h =  theHistograms.get(("ZZTo"+*it+"_Mass_01").c_str());
      
      if(!h) continue;
      for(int i = 1; i<=h->GetNbinsX();i++){
	
	Float_t Err = h->GetBinError(i);
	h->SetBinError(i,sqrt(Err*Err+hvar->GetBinContent(i)));
      }
    }
  }  
}
