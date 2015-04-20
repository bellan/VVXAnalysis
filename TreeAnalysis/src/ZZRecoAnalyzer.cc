#include "VVXAnalysis/TreeAnalysis/interface/ZZRecoAnalyzer.h"
#include "VVXAnalysis/Commons/interface/SignalDefinitions.h"
#include "VVXAnalysis/Commons/interface/Utilities.h"

#include <boost/assign/std/vector.hpp> 
#include <boost/assert.hpp> 
using namespace std;
using namespace boost::assign; // bring 'operator+=()' into scope

#include <boost/foreach.hpp>
#define foreach BOOST_FOREACH

using std::cout;
using std::endl;


using namespace phys;

// Int_t VVXAnalyzer::cut() {
  
//   return 1;
// }

void ZZRecoAnalyzer::ZZplots(int id){

  if(ZZ->id() != id && id != -1) return; // -1 here means generic 4l final state

  std::string decay  = "4l";
  
  if      (id == 52) {decay = "4m";}
  else if (id == 48) {decay = "2e2m";}
  else if (id == 44) {decay = "4e";}
    
  theHistograms.fill(std::string("ZZTo")+decay+"_Mass_01" , std::string("Invariant mass of ZZ_{1}#rightarrow ")+decay , Xbins, ZZ->mass(),theWeight);
  
  if (topology.test(0)){
    m4L_gen = sqrt((genVBParticles->at(0).p4()+genVBParticles->at(1).p4())*(genVBParticles->at(0).p4()+genVBParticles->at(1).p4()));
    theHistograms.fill(std::string("ResMat_ZZTo")+decay+"_Mass_01", std::string("Response matrix invariant mass of ZZ_{1}#rightarrow ")+decay, Xbins, Xbins, ZZ->mass() ,m4L_gen , theWeight);
  }
  
 theHistograms.fill(std::string("ZZTo")+decay+"_Mass"+"_FRVar", std::string("Var From FR Invariant mass of ZZ_{1}#rightarrow ")+decay, Xbins, ZZ->mass(),ZZ->fakeRateSFVar());

 

}


void ZZRecoAnalyzer::analyze(){

  ZZplots();   // ZZ --> 4l
  ZZplots(52); // ZZ --> 4m
  ZZplots(48); // ZZ --> 2e2m
  ZZplots(44); // ZZ --> 4e
  
  if (topology.test(0)) theHistograms.fill("PassDef", "Number of events passing the signal definition", 200, 50,  1000, ZZ->mass(),theWeight);

  else theHistograms.fill("NoPassDef", "Number of events not passing the signal definition", 200, 50,  1000, ZZ->mass(),theWeight);

}


void ZZRecoAnalyzer::begin() {

  Xbins += 100,200,250,300,350,400,500,600,800;
  m4L_gen =0;
}

void ZZRecoAnalyzer::end( TFile &) {
  
  if(region_ == phys::CR3P1F || region_ == phys::CR2P2F) {
    vector<std::string>  FinalState = {"4m","4e","2e2m"};
    
    for (std::vector<std::string>::iterator it = FinalState.begin() ; it != FinalState.end(); ++it){
      
      TH1 *hvar =  new TH1F();
      hvar =  theHistograms.GetHisto(("ZZTo"+*it+"_Mass_FRVar").c_str());
      
      TH1 *h =  new TH1F();
      h =  theHistograms.GetHisto(("ZZTo"+*it+"_Mass_01").c_str());
      
      if(!h) continue;
      for(int i = 1; i<=h->GetNbinsX();i++){
	
	Float_t Err = h->GetBinError(i);
	h->SetBinError(i,sqrt(Err*Err+hvar->GetBinContent(i)));
      }
    }
  }  
}
