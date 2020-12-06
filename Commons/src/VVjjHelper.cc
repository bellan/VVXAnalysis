#include "VVXAnalysis/Commons/interface/VVjjHelper.h"
#include "VVXAnalysis/Commons/interface/Comparators.h"
#include "VVXAnalysis/Commons/interface/Constants.h"
#include "VVXAnalysis/Commons/interface/SignalDefinitions.h"
#include "VVXAnalysis/Commons/interface/Utilities.h"
#include "VVXAnalysis/DataFormats/interface/TypeDefs.h"

#include <boost/foreach.hpp>
#define foreach BOOST_FOREACH

#include <boost/assign/std/vector.hpp> 
using namespace boost::assign;

using namespace std;
using namespace phys;

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
  } else{
    VV = helper.BuildVV(eventtype);
  }
    
  return true;
}


void VVjjHelper::LeptonSearch(vector<Particle> &genparticles, string eventtype){
  // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  // ~~~~~~~~~~~~~~~~~~ Function to find Leptons ~~~~~~~~~~~~~~~~~~~
  // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
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
    cout << "ERROR: lepton's number is: " << leptons_.size() << endl << endl;
    return;
  }
  }
}


VVtype VVjjHelper::BuildVV(string eventtype){
  // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  // ~~~~~~~~~~~~~~~~~~ Function to build DiBoson ~~~~~~~~~~~~~~~~~~
  // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  VVtype VV;
  //things

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
