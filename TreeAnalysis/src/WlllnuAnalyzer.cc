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

//double deltaRMax;
double deltaPhiMax;
int posI; //int posJ;
//TString histName; TString histTitle;

void WlllnuAnalyzer::begin(){
  genElectronsHist_ = new std::vector<phys::Particle>;
  electronsHist_ = new std::vector<phys::Particle>;
  genMuonsHist_ = new std::vector<phys::Particle>;
  muonsHist_ = new std::vector<phys::Particle>;
  genNeutrinosHist_ = new std::vector<phys::Particle>;
  metHist_ = new std::vector<phys::Particle>;
}

Int_t WlllnuAnalyzer::cut() {
  genEventSetup();
  return 1;
}

void WlllnuAnalyzer::analyze(){

  // -- Particle ID -- //
  foreach(const phys::Particle &gen, *genParticles){
    theHistograms->fill("gen_ID"      , "gen ID"  , 30 , 0, 30, gen.id(), theWeight);
  }
  
  // -- Gen neutrino pt -- //
  foreach(const phys::Particle &genNu, *genNeutrinos_){
    theHistograms->fill("GEN_nu_pt", "Gen Neutrino pt;pt", 75 , 0, 300.,  genNu.pt(), theWeight);
  }
  
  /*
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
  */
  
  
  // -- W boson 3 leptons decay possible configurations -- //
  if(genChLeptons_->size()==3){
  
    if(genElectrons_->size()==2 && genMuons_->size()==1){
      // -- Gen leptons pt -- //
      foreach(const phys::Particle &genEl, *genElectrons_){
        theHistograms->fill("GEN_el_pt_mode1", "Gen Electron pt (e+e- mu);pt"  , 75 , 0, 300.,  genEl.pt(), theWeight);
      }
      foreach(const phys::Particle &genMu, *genMuons_){
        theHistograms->fill("GEN_mu_pt_mode1", "Gen Muons pt (e+e- mu);pt"  , 75 , 0, 300.,  genMu.pt(), theWeight);
      }
    
      // -- Invariant mass of the e+e- pair -- //
      phys::Particle el1 = genElectrons_->at(0);
      phys::Particle el2 = genElectrons_->at(1);
      double elPairInvMass = pow(el1.p4().Energy() + el2.p4().Energy(),2) - pow(el1.p4().Px() + el2.p4().Px(),2) - pow(el1.p4().Py() + el2.p4().Py(),2) - pow(el1.p4().Pz() + el2.p4().Pz(),2);
      theHistograms->fill("GEN_el_pair_invariant_mass_mode1", "Gen Electron pair e+e- invariant mass", 75, 0, 300., elPairInvMass, theWeight);
      
      // -- Invariant mass of the 4 leptons e+ e- mu nu -- // 
      phys::Particle mu = genMuons_->at(0);
      phys::Particle nu = genNeutrinos_->at(0);
      double fourLepInvMass = pow(el1.p4().Energy() + el2.p4().Energy() + mu.p4().Energy() + nu.p4().Energy(),2) - pow(el1.p4().Px() + el2.p4().Px() + mu.p4().Px() + nu.p4().Px(),2) - pow(el1.p4().Py() + el2.p4().Py() + mu.p4().Py() + nu.p4().Py(),2) - pow(el1.p4().Pz() + el2.p4().Pz() + mu.p4().Pz() + nu.p4().Pz(),2);
      theHistograms->fill("GEN_four_leptons_invariant_mass_mode1", "Gen Leptons e+e-mu nu invariant mass", 75, 0, 300., fourLepInvMass, theWeight);
      
      // -- Transverse mass of the 3 charged leptons e+ e- mu -- //
      double transverseMass = mT(el1,el2,mu);
      theHistograms->fill("GEN_charged_leptons_transverse_mass_mode1", "Gen Charged Leptons e+e-mu transverse mass", 75, 0, 300., transverseMass, theWeight);
      
      // -- Reconstructed leptons efficiency e+ e- mu -- //
      if(electrons->size()==2 && muons->size()==1){
        // -- Gen Leptons & Rec Leptons pt match -- //
        reconstructionLepCompatibility(genElectrons_, electrons, "REC_el_compatibility_mode1", "Reconstructed electron pair e+e- compatibility");
        reconstructionLepCompatibility(genMuons_, muons, "REC_mu_compatibility_mode1", "Reconstructed muon mu compatibility");
        
        
        
      }
      
      
      
    }
    
    else if(genMuons_->size()==2 && genElectrons_->size()==1){
      // -- Gen leptons pt -- //
      foreach(const phys::Particle &genEl, *genElectrons_){
        theHistograms->fill("GEN_el_pt_mode2", "Gen Electron pt (mu+mu- e);pt"  , 75 , 0, 300.,  genEl.pt(), theWeight);
      }
      foreach(const phys::Particle &genMu, *genMuons_){
        theHistograms->fill("GEN_mu_pt_mode2", "Gen Muons pt (mu+mu- e);pt"  , 75 , 0, 300.,  genMu.pt(), theWeight);
      }
      
      // -- Invariant mass of the mu+mu- pair -- //
      phys::Particle mu1 = genMuons_->at(0);
      phys::Particle mu2 = genMuons_->at(1);
      double muPairInvMass = pow(mu1.p4().Energy() + mu2.p4().Energy(),2) - pow(mu1.p4().Px() + mu2.p4().Px(),2) - pow(mu1.p4().Py() + mu2.p4().Py(),2) - pow(mu1.p4().Pz() + mu2.p4().Pz(),2);
      theHistograms->fill("GEN_mu_pair_invariant_mass", "Gen Muon pair mu+mu- invariant mass", 75, 0, 300., muPairInvMass, theWeight);
      
      // -- Invariant mass of the 4 leptons mu+ mu- e nu -- // 
      phys::Particle el = genElectrons_->at(0);
      phys::Particle nu = genNeutrinos_->at(0);
      double fourLepInvMass = pow(mu1.p4().Energy() + mu2.p4().Energy() + el.p4().Energy() + nu.p4().Energy(),2) - pow(mu1.p4().Px() + mu2.p4().Px() + el.p4().Px() + nu.p4().Px(),2) - pow(mu1.p4().Py() + mu2.p4().Py() + el.p4().Py() + nu.p4().Py(),2) - pow(mu1.p4().Pz() + mu2.p4().Pz() + el.p4().Pz() + nu.p4().Pz(),2);
      theHistograms->fill("GEN_four_leptons_invariant_mass_mode2", "Gen Leptons mu+mu-e nu invariant mass", 75, 0, 300., fourLepInvMass, theWeight);      
      
      // -- Transverse mass of the 3 charged leptons mu+ mu- e -- //
      double transverseMass = mT(mu1,mu2,el);
      theHistograms->fill("GEN_charged_leptons_transverse_mass_mode2", "Gen Charged Leptons mu+mu-e transverse mass", 75, 0, 300., transverseMass, theWeight);
      
      // -- Reconstructed leptons efficiency mu+ mu- e -- //
      if(muons->size()==2 && electrons->size()==1){
        // -- Gen Leptons & Rec Leptons pt match -- //
        reconstructionLepCompatibility(genElectrons_, electrons, "REC_el_compatibility_mode2", "Reconstructed electron e compatibility");
        reconstructionLepCompatibility(genMuons_, muons, "REC_mu_compatibility_mode2", "Reconstructed muon pair mu+mu- compatibility");
      } 
      
      
      
      
      
    }
    
    else if(genElectrons_->size()==3){
    /*
      // -- Gen electron pt -- //
      foreach(const phys::Particle &genEl, *genElectrons_){
        theHistograms->fill("GEN_el_pt"      , "Gen Electron pt;pt"  , 50 , 0, 300.,  genEl.pt(), theWeight);
      }
    */
      
      
      
      
    }
    
    else if(genMuons_->size()==3){
      /*
      // -- Gen muons pt -- //
      foreach(const phys::Particle &genMu, *genMuons_){
        theHistograms->fill("GEN_mu_pt"      , "Gen Muon pt;pt"  , 50 , 0, 300.,  genMu.pt(), theWeight);
      }
      */
      
      
      
      
    }
    
  }
  
  





}
  

void WlllnuAnalyzer::end(TFile &){
  /*
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
  */








  
}


void WlllnuAnalyzer::genEventSetup(){
  
  genChLeptons_->clear();
  genElectrons_ = new std::vector<phys::Particle>;
  genMuons_ = new std::vector<phys::Particle>;
  genNeutrinos_->clear();

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
  }
  
  std::sort(genChLeptons_    ->begin(), genChLeptons_    ->end(), [](const Particle& a, const Particle& b){ return a.pt() < b.pt(); });
  std::sort(genElectrons_    ->begin(), genElectrons_    ->end(), [](const Particle& a, const Particle& b){ return a.pt() < b.pt(); });
  std::sort(genMuons_    ->begin(), genMuons_    ->end(), [](const Particle& a, const Particle& b){ return a.pt() < b.pt(); });
  std::sort(genNeutrinos_    ->begin(), genNeutrinos_    ->end(), [](const Particle& a, const Particle& b){ return a.pt() < b.pt(); });
  
}


void WlllnuAnalyzer::reconstructionLepCompatibility(std::vector<phys::Particle>* genLep, std::vector<phys::Lepton>* recLep, string histName, string histTitle){
  
  for(int i=0;i<genLep->size();i++){
    	double deltaRMax = 0.; int posJ = 0;
    	bool usedLep[recLep->size()] = {};
    	for(int j=0;j<recLep->size();j++){
          usedLep[j] = false;
    	}
    	for(int j=0;j<recLep->size();j++){
      	  if(usedLep[j]==false){
            if( deltaR(genLep->at(i),recLep->at(j)) < 0.1 ){ 
	      if(deltaRMax==0.){
	        deltaRMax = deltaR(genLep->at(i),recLep->at(j));
	        posJ = j;
	      }
              else if(deltaR(genLep->at(i),recLep->at(j))<deltaRMax){
	        deltaRMax = deltaR(genLep->at(i),recLep->at(j));
	        posJ = j;
	      }
            }
          }
        }
        if(deltaRMax != 0){ 
          usedLep[posJ] = true;
          theHistograms->fill(histName, histTitle, 75, 0., 300., genLep->at(i).pt(), theWeight);
        }
  }
  
}
