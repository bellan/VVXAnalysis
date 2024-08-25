#include "VVXAnalysis/TreeAnalysis/interface/WlllnuAnalyzer.h"
#include "VVXAnalysis/Commons/interface/SignalDefinitions.h"
#include "VVXAnalysis/Commons/interface/Utilities.h"
#include "VVXAnalysis/Commons/interface/Constants.h"
#include "VVXAnalysis/Commons/interface/Colours.h"
#include "VVXAnalysis/DataFormats/interface/TypeDefs.h"
#include "VVXAnalysis/DataFormats/interface/Boson.h"
#include "VVXAnalysis/DataFormats/interface/DiBoson.h"
#include "VVXAnalysis/DataFormats/interface/Lepton.h"
#include "VVXAnalysis/Commons/interface/Comparators.h"

#include <TCanvas.h>
#include <TH1F.h>

#include <boost/foreach.hpp>
#define foreach BOOST_FOREACH

#include <boost/assign/std/vector.hpp>
using namespace boost::assign;
using namespace colour;

using std::cout;
using std::endl;


using namespace phys;
using namespace physmath;

int start;

std::vector<phys::Particle>* genElectrons_;
std::vector<phys::Particle>* genElectronsHist_;
std::vector<phys::Particle>* electronsHist_;

std::vector<phys::Particle>* genMuons_;
//std::vector<phys::Particle>* muonsHist_;

double ptEffNum; double ptEffDen;
double deltaRMax;
int posI; int posJ;
int numWrongCharge;

void WlllnuAnalyzer::begin(){
  
}

Int_t WlllnuAnalyzer::cut() {
  genEventSetup();
  return 1;
}

void WlllnuAnalyzer::analyze(){

  // cout << "------------------------------------------------------------------------"<<endl;
  // cout << "Run: " << run << " event: " << event << endl;

  // -- Particle ID -- //
  foreach(const phys::Particle &gen, *genParticles){
    theHistograms->fill("gen_ID"      , "gen ID"  , 30 , 0, 30, gen.id(), theWeight);
  }
  
  // -- Gen electron pt -- //
  foreach(const phys::Particle &genEl, *genElectrons_){
    theHistograms->fill("GEN_el_pt"      , "Gen Electron pt;pt"  , 50 , 0, 300.,  genEl.pt(), theWeight);
  }
  
  // -- Gen photon pt -- //
  /*if(genPhotons_->size() > 0){
    theHistograms->fill("GEN_photon_pt", "GEN photon pt;pt;Counts" , 60 , 0., 2., genPhotons_->at(0).pt(), theWeight);
  }*/
  
  // -- Gen electrons & electrons pt match -- //
  /*for(int i=0;i<genElectrons_->size();i++){
    ptEffDen += 1;
    for(int j=0;j<electrons->size();j++){
      if(physmath::deltaR(genElectrons_->at(i),electrons->at(j)) < 0.1 ){ 
	ptEffNum += 1;
	theHistograms->fill("gen_rec_el_pt"      , "Gen electrons & electrons pt match;pt"  , 50 , 0., 300., genElectrons_->at(i).pt(), theWeight);
      }
    }
  }*/
  
  for(int i=0;i<genElectrons_->size();i++){
    deltaRMax = 0.;
    ptEffDen += 1;
    for(int j=0;j<electrons->size();j++){
      if( deltaR(genElectrons_->at(i),electrons->at(j)) < 0.1 ){ 
	if(deltaRMax==0.){
	  deltaRMax = deltaR(genElectrons_->at(i),electrons->at(j));
	  posJ = j;
	}
        else if(deltaR(genElectrons_->at(i),electrons->at(j))<deltaRMax){
	  deltaRMax = deltaR(genElectrons_->at(i),electrons->at(j));
	  posJ = j;
	}
      }
    }
    if(deltaRMax != 0){ 
      ptEffNum += 1;
      genElectronsHist_->push_back(genElectrons_->at(i));
      electronsHist_->push_back(electrons->at(posJ));
      theHistograms->fill("gen_rec_el_pt"      , "Gen electrons & electrons pt match;pt"  , 50 , 0., 300., genElectrons_->at(i).pt(), theWeight);
    }
  }
  start += 1;
  
  double ptEff = ptEffNum/ptEffDen;
  // cout << "-----------------------------------------------------------------" << endl;
  // cout << "genElectrons_ vs electrons: pt Efficiency = " << ptEff << endl;
  
  /*
  // -- Gen electrons & electrons charge match -- //  
  if(genElectronsHist_->size()>0 && electronsHist_->size()>0){
    for(int i=0;i<genElectronsHist_->size();i++){
      if( genElectrons_->at(i).charge() != electrons->at(i).charge() ){
	numWrongCharge += 1;
      }
    }
  }
  // cout << "--------------------------------------------------------------------" << endl;
  // cout << "Number of wrong charge = " << numWrongCharge << endl;
}
  */

void WlllnuAnalyzer::end(TFile &){
}


void WlllnuAnalyzer::genEventSetup(){
  
  //genQuarks_->clear();
  genChLeptons_->clear();

  genElectrons_ = new std::vector<phys::Particle>;
  electronsHist_ = new std::vector<phys::Particle>;
  genMuons_ = new std::vector<phys::Particle>;
  
  //genNeutrinos_->clear();
  genPhotons_->clear();
  genPhotonsPrompt_->clear();
  
 
  /*
  // -- Sort gen particles -- //
  for(auto p : *genParticles){
    unsigned int aPID = abs(p.id());
    if(aPID < 9){
      genQuarks_->push_back(p);
    }
    else if(aPID == 11 || aPID == 13){
      genChLeptons_->push_back(p);
    }
    else if(aPID == 12 || aPID == 14){
      genNeutrinos_->push_back(p);
    }
    else if(aPID == 22){
      genPhotons_->push_back(p);
      if(p.genStatusFlags().test(phys::isPrompt)){
	genPhotonsPrompt_->push_back(p);
      }
    }
  }
  */

  // -- Sort genChLeptons -- //
  for(auto p : *genParticles){
    unsigned int aPID = abs(p.id());
    if(aPID == 11 || aPID == 13){
      genChLeptons_->push_back(p);
      // cout << "------------------------------------------------------------------------"<<endl;
      // cout << "Charged Lepton pt: " << p.pt() << endl;
      if(aPID == 11){
	genElectrons_->push_back(p);
	// cout << "------------------------------------------------------------------------"<<endl;
        // cout << "Electron pt: " << p.pt() << endl;
      }
      else{
	genMuons_->push_back(p);
	// cout << "------------------------------------------------------------------------"<<endl;
        // cout << "Muon pt: " << p.pt() << endl;
      }
    }
    else if(aPID == 22){
      genPhotons_->push_back(p);
      if(p.genStatusFlags().test(phys::isPrompt)){
	genPhotonsPrompt_->push_back(p);
        // cout << "------------------------------------------------------------------------"<<endl;
        // cout << "Photon pt: " << p.pt() << endl;
      }
    }
  }

  
  // std::sort(genQuarks_       ->begin(), genQuarks_       ->end(), [](const Particle& a, const Particle& b){ return a.pt() < b.pt(); });
  std::sort(genChLeptons_    ->begin(), genChLeptons_    ->end(), [](const Particle& a, const Particle& b){ return a.pt() < b.pt(); });

  std::sort(genElectrons_    ->begin(), genElectrons_    ->end(), [](const Particle& a, const Particle& b){ return a.pt() < b.pt(); });
  std::sort(genMuons_    ->begin(), genMuons_    ->end(), [](const Particle& a, const Particle& b){ return a.pt() < b.pt(); });

  // std::sort(genNeutrinos_    ->begin(), genNeutrinos_    ->end(), [](const Particle& a, const Particle& b){ return a.pt() < b.pt(); });
  std::sort(genPhotons_      ->begin(), genPhotons_      ->end(), [](const Particle& a, const Particle& b){ return a.pt() < b.pt(); });
  std::sort(genPhotonsPrompt_->begin(), genPhotonsPrompt_->end(), [](const Particle& a, const Particle& b){ return a.pt() < b.pt(); });

  if(start==0){
    genElectronsHist_ = new std::vector<phys::Particle>;
    electronsHist_ = new std::vector<phys::Particle>;
  }
  //start +=1;


}

