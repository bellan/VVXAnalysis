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

std::vector<phys::Particle>* genElectrons_;
std::vector<phys::Particle>* genElectronsHist_;
std::vector<phys::Particle>* electronsHist_;
std::vector<phys::Particle>* genMuons_;
std::vector<phys::Particle>* genMuonsHist_;
std::vector<phys::Particle>* muonsHist_;
std::vector<phys::Particle>* genNeutrinosHist_;
std::vector<phys::Particle>* metHist_;

double deltaRMax;
double deltaPhiMax;
int posI; int posJ;
int start;
//TString histName; TString histTitle;

void WlllnuAnalyzer::begin(){
  /*genElectronsHist_ = new std::vector<phys::Particle>;
  electronsHist_ = new std::vector<phys::Particle>;
  genMuonsHist_ = new std::vector<phys::Particle>;
  muonsHist_ = new std::vector<phys::Particle>;
  genNeutrinosHist_ = new std::vector<phys::Particle>;
  metHist_ = new std::vector<phys::Particle>;*/
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

  // -- Gen muons pt -- //
  foreach(const phys::Particle &genMu, *genMuons_){
    theHistograms->fill("GEN_mu_pt"      , "Gen Muon pt;pt"  , 50 , 0, 300.,  genMu.pt(), theWeight);
  }
  
  // -- Gen photon pt -- //
  /*if(genPhotons_->size() > 0){
    theHistograms->fill("GEN_photon_pt", "GEN photon pt;pt;Counts" , 60 , 0., 2., genPhotons_->at(0).pt(), theWeight);
  }*/
  
  // -- Gen electrons & electrons pt match -- //
  for(int i=0;i<genElectrons_->size();i++){
    deltaRMax = 0.;
    bool usedElectrons_[electrons->size()] = {};
    for(int j=0;j<electrons->size();j++){
      usedElectrons_[j] = false;
    }
    
    for(int j=0;j<electrons->size();j++){
      if(usedElectrons_[j]==false){
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
    }
    if(deltaRMax != 0){ 
      genElectronsHist_->push_back(genElectrons_->at(i));
      electronsHist_->push_back(electrons->at(posJ));
      usedElectrons_[posJ] = true;
      theHistograms->fill("GEN_REC_el_pt"      , "Gen electrons & electrons pt match;pt"  , 50 , 0., 300., genElectrons_->at(i).pt(), theWeight);
      
    }
  }
  // -- Gen muons & muons pt match -- //
  for(int i=0;i<genMuons_->size();i++){
    deltaRMax = 0.;
    bool usedMuons_[muons->size()] = {};
    for(int j=0;j<muons->size();j++){
      usedMuons_[j] = false;
    }
    
    for(int j=0;j<muons->size();j++){
      if(usedMuons_[j]==false){
        if( deltaR(genMuons_->at(i),muons->at(j)) < 0.1 ){ 
	  if(deltaRMax==0.){
	    deltaRMax = deltaR(genMuons_->at(i),muons->at(j));
	    posJ = j;
	  }
          else if(deltaR(genMuons_->at(i),muons->at(j))<deltaRMax){
	    deltaRMax = deltaR(genMuons_->at(i),muons->at(j));
	    posJ = j;
	  }
        }
      }
    }
    if(deltaRMax != 0){ 
      genMuonsHist_->push_back(genMuons_->at(i));
      muonsHist_->push_back(muons->at(posJ));
      usedMuons_[posJ] = true;
      theHistograms->fill("GEN_REC_mu_pt"      , "Gen muons & muons pt match;pt"  , 50 , 0., 300., genMuons_->at(i).pt(), theWeight);
      /*if( genMuons_->at(i).eta() < 1.4 ){                         // Barrel region
	theHistograms->fill("GEN_REC_mu_pt_barrel_region"      , "Gen muons & muons pt match (BARREL REGION);pt"  , 50 , 0., 300., genMuons_->at(i).pt(), theWeight);
      }
      else if( 1.4 < genMuons_->at(i).eta() < 1.6 ){              // Endcap region         
	theHistograms->fill("GEN_REC_mu_pt_endcap_region"      , "Gen muons & muons pt match (ENDCAP REGION);pt"  , 50 , 0., 300., genMuons_->at(i).pt(), theWeight);
      }
      else{                                                             // Other regions external to endcap region
	theHistograms->fill("GEN_REC_mu_pt_other_regions"      , "Gen muons & muons pt match (OTHER REGIONS EXTERNAL TO ENDCAP);pt"  , 50 , 0., 300., genMuons_->at(i).pt(), theWeight);
      }*/
    }
  }

  // -- Gen neutrinos & MET pt match -- //
  deltaPhiMax = 0.;
  for(int i=0;i<genNeutrinos_->size();i++){
    if( deltaPhi(genNeutrinos_->at(i).phi(),met->phi()) < 0.1 ){ 
      if(deltaPhiMax==0.){
	deltaPhiMax = deltaPhi(genNeutrinos_->at(i).phi(),met->phi());
	posI = i;
      }
      else if(deltaPhi(genNeutrinos_->at(i).phi(),met->phi())<deltaPhiMax){
	deltaPhiMax = deltaPhi(genNeutrinos_->at(i).phi(),met->phi());
	posI = i;
      }
    }
  }
  if(deltaPhiMax != 0){ 
    //genNeutrinosHist_->push_back(genNeutrinos_->at(posI));
    //metHist_->push_back(met);
    theHistograms->fill("GEN_REC_neutrinos_pt"      , "Gen neutrinos & MET pt match;pt"  , 50 , 0., 300., genNeutrinos_->at(posI).pt(), theWeight);
  }
  
  // -- W boson possible configurations -- //
  if(genChLeptons_->size()==3){
    if(genElectrons_->size()==2 && genMuons_->size()==1){
      //histName = "GEN_W_decay_ch_leptons_pt_mode1";  // 2 elettroni, 1 muone
      //histTitle = "Gen Charged Leptons pt (2 electrons & 1 muon)";
      foreach(const phys::Particle &genCh, *genChLeptons_){
	theHistograms->fill("GEN_W_decay_ch_leptons_pt_mode1", "Gen Charged Leptons pt (2 electrons & 1 muon)" , 50 , 0, 300.,  genCh.pt(), theWeight);
      }
    }
    else if(genMuons_->size()==2 && genElectrons_->size()==1){
      //histName = "GEN_W_decay_ch_leptons_pt_mode2";  // 2 muoni, 1 elettrone
      //histTitle = "Gen Charged Leptons pt (2 muons & 1 electrons)";
      foreach(const phys::Particle &genCh, *genChLeptons_){
	theHistograms->fill("GEN_W_decay_ch_leptons_pt_mode2", "Gen Charged Leptons pt (2 muons & 1 electron)" , 50 , 0, 300.,  genCh.pt(), theWeight);
      }
    }
    else if(genElectrons_->size()==3){
      //histName = "GEN_W_decay_ch_leptons_pt_mode3";  // 3 elettroni
      //histTitle = "Gen Charged Leptons pt (3 electrons)";
      foreach(const phys::Particle &genCh, *genChLeptons_){
        theHistograms->fill("GEN_W_decay_ch_leptons_pt_mode3", "Gen Charged Leptons pt (3 electrons)" , 50 , 0, 300.,  genCh.pt(), theWeight);
      }
    }
    else if(genMuons_->size()==3){
      //histName = "GEN_W_decay_ch_leptons_pt_mode4";  // 3 muoni
      //histTitle = "Gen Charged Leptons pt (3 muons)";
      foreach(const phys::Particle &genCh, *genChLeptons_){
        theHistograms->fill("GEN_W_decay_ch_leptons_pt_mode4", "Gen Charged Leptons pt (3 muons)" , 50 , 0, 300.,  genCh.pt(), theWeight);
      }
    }
    //foreach(const phys::Particle &genCh, *genChLeptons_){
    //  theHistograms->fill(histName, histTitle , 50 , 0, 300.,  genCh.pt(), theWeight);
    //}
  }
  
  





}
  

void WlllnuAnalyzer::end(TFile &){
  
  // -- genElectrons pt Resolution -- // 
  for(int i=0;i<genElectronsHist_->size();i++){
    if( genElectronsHist_->at(i).eta() < 1.4 ){                         // Barrel region
      theHistograms->fill("GEN_el_pt_resolution_barrel_region" , "GEN electron pt resolution (BARREL REGION)"  , 30 , -10, 10, genElectronsHist_->at(i).pt() - electronsHist_->at(i).pt(), theWeight);
    }
    else if( 1.4 < genElectronsHist_->at(i).eta() < 3.0 ){              // Endcap region         
      theHistograms->fill("GEN_el_pt_resolution_endcap_region" , "GEN electron pt resolution (ENDCAP REGION)"  , 30 , -10, 10, genElectronsHist_->at(i).pt() - electronsHist_->at(i).pt(), theWeight);
    }
    else if( 3.0 < genElectronsHist_->at(i).eta() < 5.0 ){              // Other regions external to endcap region
      theHistograms->fill("GEN_el_pt_resolution_forward_region" , "GEN electron pt resolution (FORWARD REGION)"  , 30 , -10, 10, genElectronsHist_->at(i).pt() - electronsHist_->at(i).pt(), theWeight);
    }
  }
  // -- genMuons pt Resolution -- //
  for(int i=0;i<genMuonsHist_->size();i++){
    if( genMuonsHist_->at(i).eta() < 1.4 ){                         // Barrel region
      theHistograms->fill("GEN_mu_pt_resolution_barrel_region"      , "GEN muon pt resolution (BARREL REGION)"  , 30 , -10, 10, genMuonsHist_->at(i).pt() - muonsHist_->at(i).pt(), theWeight);
    }
    else if( 1.4 < genMuonsHist_->at(i).eta() < 3.0 ){              // Endcap region         
      theHistograms->fill("GEN_mu_pt_resolution_endcap_region"      , "GEN muon pt resolution (ENDCAP REGION)"  , 30 , -10, 10, genMuonsHist_->at(i).pt() - muonsHist_->at(i).pt(), theWeight);
    }
    else if( 3.0 < genMuonsHist_->at(i).eta() < 5.0 ){              // Other regions external to endcap region
      theHistograms->fill("GEN_mu_pt_resolution_forward_region"      , "GEN muon pt resolution (FORWARD REGION)"  , 30 , -10, 10, genMuonsHist_->at(i).pt() - muonsHist_->at(i).pt(), theWeight);
    }
  }
  








  
}


void WlllnuAnalyzer::genEventSetup(){
  
  //genQuarks_->clear();
  genChLeptons_->clear();
  genElectrons_ = new std::vector<phys::Particle>;
  genMuons_ = new std::vector<phys::Particle>;
  genNeutrinos_->clear();
  genPhotons_->clear();
  genPhotonsPrompt_->clear();

  
  if(start==0){
    genElectronsHist_ = new std::vector<phys::Particle>;
    electronsHist_ = new std::vector<phys::Particle>;
    genMuonsHist_ = new std::vector<phys::Particle>;
    muonsHist_ = new std::vector<phys::Particle>;
    genNeutrinosHist_ = new std::vector<phys::Particle>;
    metHist_ = new std::vector<phys::Particle>;
  }
 
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

  // -- Sort genChLeptons, genNeutrinos & genPhotons -- //
  for(auto p : *genParticles){
    unsigned int aPID = abs(p.id());
    if(aPID == 11 || aPID == 13){
      genChLeptons_->push_back(p);
      if(aPID == 11){
	genElectrons_->push_back(p);
      }
      else{
	genMuons_->push_back(p);
      }
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
  
  // std::sort(genQuarks_       ->begin(), genQuarks_       ->end(), [](const Particle& a, const Particle& b){ return a.pt() < b.pt(); });
  std::sort(genChLeptons_    ->begin(), genChLeptons_    ->end(), [](const Particle& a, const Particle& b){ return a.pt() < b.pt(); });
  std::sort(genElectrons_    ->begin(), genElectrons_    ->end(), [](const Particle& a, const Particle& b){ return a.pt() < b.pt(); });
  std::sort(genMuons_    ->begin(), genMuons_    ->end(), [](const Particle& a, const Particle& b){ return a.pt() < b.pt(); });
  std::sort(genNeutrinos_    ->begin(), genNeutrinos_    ->end(), [](const Particle& a, const Particle& b){ return a.pt() < b.pt(); });
  std::sort(genPhotons_      ->begin(), genPhotons_      ->end(), [](const Particle& a, const Particle& b){ return a.pt() < b.pt(); });
  std::sort(genPhotonsPrompt_->begin(), genPhotonsPrompt_->end(), [](const Particle& a, const Particle& b){ return a.pt() < b.pt(); });





  
}
