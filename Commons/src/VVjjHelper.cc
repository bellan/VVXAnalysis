#include "VVXAnalysis/Commons/interface/VVjjHelper.h"
#include "VVXAnalysis/Commons/interface/Comparators.h"
#include "VVXAnalysis/Commons/interface/Constants.h"
#include "VVXAnalysis/Commons/interface/SignalDefinitions.h"
#include "VVXAnalysis/Commons/interface/Utilities.h"
#include "VVXAnalysis/DataFormats/interface/TypeDefs.h"
#include "VVXAnalysis/TreeAnalysis/interface/Histogrammer.h"

#include <boost/foreach.hpp>
#define foreach BOOST_FOREACH

#include <boost/assign/std/vector.hpp> 
using namespace boost::assign;

using namespace std;
using namespace phys;
using namespace physmath;

void VVjjHelper::test() {
  cout << "pippo!!" << endl;
}


bool VVjjHelper::FindDiBoson(vector<Particle> &genparticles, VVtype &VV, string eventtype){
  // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  // ~~~~~~~~~~~~~~~~~~ Function to find DiBoson ~~~~~~~~~~~~~~~~~~~
  // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  VVjjHelper helper;
  
  // Search for and get leptons in genParticles
  helper.LeptonSearch(genparticles, eventtype);

  // Check if there are enough leptons and neutrinos
  unsigned int leptonnumber = helper.GetAllLeptonsNumber();
  unsigned int neutrinonumber = helper.GetNeutrinosNumber();

  if(leptonnumber != 4 || neutrinonumber > 1){
    return false;
  }
  else{    
    VV = helper.BuildVV(eventtype);
  }

  if(VV.first().mass() == 0){
    return false;
  }
    
  return true;
}


void VVjjHelper::LeptonSearch(vector<Particle> &genparticles, string eventtype){
  // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  // ~~~~~~~~~~~~~~~~~~ Function to find Leptons ~~~~~~~~~~~~~~~~~~~
  // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  Histogrammer theHistograms;
  
  // Reset Data Member containers
  leptons_.clear();   leptons_.shrink_to_fit();
  neutrinos_.clear(); neutrinos_.shrink_to_fit();
  
  // Get leptons from genParticles
  foreach(const Particle &gen, genparticles){
    if((abs(gen.id()) != 11 && abs(gen.id()) != 13 && abs(gen.id()) != 12 && abs(gen.id()) != 14) || (!(gen.genStatusFlags().test(phys::GenStatusBit::isPrompt)) || !(gen.genStatusFlags().test(phys::GenStatusBit::fromHardProcess)))) 
      continue;

    // Leptons accepted are:
    // - electrons with pt >= 7. and |eta| <= 2.5
    if(abs(gen.id()) == 11 && gen.pt() >= 7. && abs(gen.eta()) <= 2.5){
      leptons_.push_back(gen);
    }
    // - muons with pt >= 5. and |eta| <= 2.4    
    else if(abs(gen.id()) == 13 && gen.pt() >= 5. && abs(gen.eta()) <= 2.4){
      leptons_.push_back(gen);
    }
    // - electronic or muonic neutrinos
    else if(abs(gen.id()) == 12 || abs(gen.id()) == 14){
      neutrinos_.push_back(gen);
    }
  }

  // sorting leptons by pt
  sort(leptons_.begin(), leptons_.end(), PtComparator());

  // different selection for differents event types: 1 for WZ and 0 for ZZ
  int switchcase = -1;
  if(eventtype == "WZ" && leptons_.size() == 3){
    switchcase = 1;
  }
  else if(eventtype == "ZZ" && leptons_.size() == 4){
    switchcase = 0;
  }

  switch(switchcase){
  case 0:{// ZZevent
    // things for ZZ event
  }
    break;
   
  case 1:{// WZ event
    // first lepton must have at least 20GeV pt
    if(leptons_[0].pt() < 20.){
      leptons_.clear();
      leptons_.shrink_to_fit();
    }

    // second lepton must have at least 12GeV pt if electron or 10GeV pt if muon
    switch(abs(leptons_[1].id())){
    case 11:{
      if(leptons_[1].pt() < 12.){
	leptons_.clear();
	leptons_.shrink_to_fit();
      }
    }
      break;

    case 13:{
      if(leptons_[1].pt() < 10.){
	leptons_.clear();
	leptons_.shrink_to_fit();
      }
    }
      break;

    default:{
      cout << "ERROR: second lepton's ID is: " << leptons_[1].id() << endl << endl;
      return;
    }
    }
  }
    break;
     
  default:{
    //cout << "ERROR: lepton's number is: " << leptons_.size() << endl << endl;
    return;
  }
  }

  TLorentzVector Ptot;

  foreach(const Particle neu, neutrinos_)
    Ptot += neu.p4();
  
  foreach(const Particle lep, leptons_)
    Ptot += lep.p4();
  
  //theHistograms.fill("AllGenlllnu_mass",   "m 3 leptons and #nu",     150, 0, 1500, Ptot.M());
  //theHistograms.fill("AllGenlllnu_trmass", "m_{T} 3 leptons and #nu", 150, 0, 1500, Ptot.Mt());
  
  if(Ptot.M() < 165.){
    leptons_.clear();
    leptons_.shrink_to_fit();    
  }
}


VVtype VVjjHelper::BuildVV(string eventtype){
  // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  // ~~~~~~~~~~~~~~~~~~ Function to build DiBoson ~~~~~~~~~~~~~~~~~~
  // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  Histogrammer theHistograms;  
  VVtype VV;
  VVtype null;

  // different ways to build up the VV couple if event is WZ or ZZ
  int switchcase = -1;
  if(eventtype == "WZ"){
    switchcase = 1;
  }
  else if(eventtype == "ZZ"){
    switchcase = 0;
  }
  
  switch(switchcase){
  case 0:{// ZZ event
    // things
  }
    break;

  case 1:{// WZ event
    // useful variables
    vector<Zltype> Zls;
    vector<VVtype> WZs;
    double differenceZ = 0;
    double differenceW = 0;
    unsigned int choice = 0;

    // fill Zl candidates
    for(int i = 0; i < (int)leptons_.size() -1; i++){
      for(int j = i + 1; j < (int)leptons_.size(); j++){
	for(int k = 0; k < (int)leptons_.size(); k++){
	  if(k != i && k != j){
	    if(leptons_[i].id() == -leptons_[j].id()){// same flavour, opposite charge
	      Zls.push_back(Zltype(BosonParticle(leptons_[i], leptons_[j], 23), leptons_[k]));
	    }
	  }
	}
      }
    }
    
    // Z is made up of the couple which gives a better Zmass 
    sort(Zls.begin(), Zls.end(), pairMassComparator(0, ZMASS));
    differenceZ = fabs(ZMASS - Zls[0].first.p4().M());

    for(int i = 0; i < (int)Zls.size(); i++){
      WZs.push_back(VVtype(BosonParticle(Zls[i].second, neutrinos_[0], copysign(24, Zls[i].second.charge())), Zls[i].first));
    }
    
    // W is made up of the couple which gives a better Wmass 
    sort(WZs.begin(), WZs.end(), pairMassComparator(0, WMASS));
    differenceW = fabs(WMASS - WZs[0].first().p4().M());
    
    // Best couple has less difference in mass from main boson
    if(differenceZ < differenceW){ // Z is better
      VV = VVtype(BosonParticle(Zls[0].second, neutrinos_[0], copysign(24, Zls[0].second.charge())), Zls[0].first);
      choice = 1;
      //theHistograms.fill("Helper_choosingWZ", "best WZ couple choice", 10, 0.5, 10.5, choice);
    }
    else{ // W is better
      VV = WZs[0];
      choice = 2;
      //theHistograms.fill("Helper_choosingWZ", "best WZ couple choice", 10, 0.5, 10.5, choice);
    }

    if(fabs(VV.first().mass() - WMASS) > rangeVmass || fabs(VV.second().mass() - ZMASS) > rangeVmass){
      VV = null;
      choice = 3;
      //theHistograms.fill("Helper_choosingWZ", "best WZ couple choice", 10, 0.5, 10.5, choice);
    }

    if(isTheSame(Zls[0].first, WZs[0].second())){
      choice = 4;
      //theHistograms.fill("Helper_choosingWZ", "best WZ couple choice", 10, 0.5, 10.5, choice);
    }

  }
    break;

  default:{
    cout << "ERROR: wrong initialisation for eventtype, " << eventtype << endl << endl;
  }
  }

  return VV;
}


  // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  // ~~~~~~~~~~~~~~~~~~~~~~ Getter functions ~~~~~~~~~~~~~~~~~~~~~~~
  // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
unsigned int VVjjHelper::GetAllLeptonsNumber(){
  return leptons_.size() + neutrinos_.size();
}

unsigned int VVjjHelper::GetNeutrinosNumber(){
  return neutrinos_.size();
}
