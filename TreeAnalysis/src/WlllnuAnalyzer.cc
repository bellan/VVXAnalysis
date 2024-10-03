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
std::vector<phys::Particle>* genMuons_;

 

//double deltaPhiMax;
//int posI;


void WlllnuAnalyzer::begin(){
  genElectrons_ = new std::vector<phys::Particle>;
  genMuons_ = new std::vector<phys::Particle>;  
  
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
    theHistograms->fill("GEN_nu_pt", "Gen Neutrino pt;pt", 75 , 0, 150.,  genNu.pt(), theWeight);
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
    theHistograms->fill("GEN_REC_neutrinos_pt"      , "Gen neutrinos & MET pt match;pt"  , 50 , 0., 150., genNeutrinos_->at(posI).pt(), theWeight);
  }
  */
  
  
  // -- W boson 3 leptons decay possible configurations -- //
  if(genChLeptons_->size()==3){
    
    if(genElectrons_->size()==2 && genMuons_->size()==1){  // MODE 1 e+ e- mu nu
      phys::Particle genEl1 = genElectrons_->at(0);
      phys::Particle genEl2 = genElectrons_->at(1);
      phys::Particle genMu = genMuons_->at(0);
      phys::Particle genNu = genNeutrinos_->at(0);

      // -- Gen leptons pt -- //
      foreach(const phys::Particle &genEl, *genElectrons_){
        theHistograms->fill("GEN_el_pt_mode1", "Gen Electron pt (e+e- mu);pt"  , 75 , 0, 150.,  genEl.pt(), theWeight);
      }
      foreach(const phys::Particle &genMu, *genMuons_){
        theHistograms->fill("GEN_mu_pt_mode1", "Gen Muons pt (e+e- mu);pt"  , 75 , 0, 150.,  genMu.pt(), theWeight);
      }
    
      // -- Invariant mass of the e+e- pair -- //
      double genElPairInvMass = (genEl1.p4() + genEl2.p4()).M();
      //double genElPairInvMass = sqrt( pow(el1.p4().Energy() + el2.p4().Energy(),2) - pow(el1.p4().Px() + el2.p4().Px(),2) - pow(el1.p4().Py() + el2.p4().Py(),2) - pow(el1.p4().Pz() + el2.p4().Pz(),2) );
      theHistograms->fill("GEN_el_pair_invariant_mass_mode1", "Gen Electron pair e+e- invariant mass (mode 1)", 75, 0, 150., genElPairInvMass, theWeight);
      
      // -- Invariant mass of the 4 gen leptons e+ e- mu nu -- // 
      double genFourLepInvMass = (genEl1.p4() + genEl2.p4() + genMu.p4() + genNu.p4()).M();
      //double genFourLepInvMass = sqrt( pow(el1.p4().Energy() + el2.p4().Energy() + mu.p4().Energy() + nu.p4().Energy(),2) - pow(el1.p4().Px() + el2.p4().Px() + mu.p4().Px() + nu.p4().Px(),2) - pow(el1.p4().Py() + el2.p4().Py() + mu.p4().Py() + nu.p4().Py(),2) - pow(el1.p4().Pz() + el2.p4().Pz() + mu.p4().Pz() + nu.p4().Pz(),2) );
      theHistograms->fill("GEN_four_leptons_invariant_mass_mode1", "Gen Leptons e+e-mu nu invariant mass (mode 1)", 150, 0, 300., genFourLepInvMass, theWeight);
      
      // -- Transverse mass of the 4 gen leptons e+ e- mu nu -- //
      double genFourLepTransverseMass = (genEl1.p4() + genEl2.p4() + genMu.p4() + genNu.p4()).Mt();;
      theHistograms->fill("GEN_four_leptons_transverse_mass_mode1", "Gen Four Leptons e+e-mu nu transverse mass (mode 1)", 150, 0., 300., genFourLepTransverseMass, theWeight);
      
      // -- Invariant mass of the gen electron pair e+e- vs Invariant mass of the 4 gen leptons e+e-mu nu -- //
      theHistograms->fill("GEN_el_pair_four_lep_inv_mass_mode1", "Gen Electron pair vs GEN Four Leptons invariant mass (mode 1)", 60, 0., 300., 60, 0., 300., genFourLepInvMass, genElPairInvMass, theWeight);
      
      // -- Gen Neutrino pt vs Invariant mass of the 4 gen leptons e+e-mu nu -- //
      theHistograms->fill("GEN_nu_pt_four_lep_inv_mass_mode1", "Gen Neutrino pt vs GEN Four Leptons invariant mass (mode 1)", 150, 0., 300., 150, 0., 300., genFourLepInvMass, genNu.pt(), theWeight);
      
      // -- Gen Muon pt vs Invariant mass of the 4 gen leptons e+e-mu nu -- //
      theHistograms->fill("GEN_mu_pt_four_lep_inv_mass_mode1", "Gen Muon pt vs GEN Four Leptons invariant mass (mode 1)", 250, 0., 500., 150, 0., 300., genFourLepInvMass, genMu.pt(), theWeight);
      
      // -- Transverse mass of the gen lepton pair mu nu vs Invariant mass of the 4 gen leptons e+e-mu nu -- //
      double genLepPairTransverseMass = (genMu.p4() + genNu.p4()).Mt();
      theHistograms->fill("GEN_lep_pair_transv_mass_mode1", "Gen Lepton pair mu nu transverse mass (mode 1)", 75, 0., 150., genLepPairTransverseMass, theWeight);
      theHistograms->fill("GEN_lep_pair_transv_mass_four_lep_inv_mass_mode1", "Gen Lepton pair mu nu transverse mass vs Four Leptons invariant mass (mode 1)", 250, 0., 500., 75, 0., 150., genFourLepInvMass, genLepPairTransverseMass, theWeight);
      
      // -- Transverse vs Invariant mass of the 4 gen leptons e+e-mu nu -- // 
      theHistograms->fill("GEN_four_lep_transv_mass_four_lep_inv_mass_mode1", "Gen Four Lepton transverse vs invariant mass (mode 1); Inv mass; Transv mass", 100, 0., 500., 100, 0., 500., genFourLepInvMass, genFourLepTransverseMass, theWeight);
      
      // -- RECONSTRUCTED LEPTONS EFFICIENCY e+ e- mu -- //
      if(electrons->size()==2 && muons->size()==1){
        phys::Particle recEl1 = electrons->at(0);
        phys::Particle recEl2 = electrons->at(1);
        phys::Particle recMu = muons->at(0);
        phys::Particle recNu = *met;
      
      	// -- Gen Leptons & Rec Leptons match -- //
        foreach(const phys::Particle &genEl, *genElectrons_){
          theHistograms->fill("REC_el_compatibility_mode1", "Reconstructed electron pair e+e- compatibility", 75, 0, 150.,  genEl.pt(), theWeight);
        }
        foreach(const phys::Particle &genMu, *genMuons_){
          theHistograms->fill("REC_mu_compatibility_mode1", "Reconstructed muon mu compatibility", 75, 0, 150.,  genMu.pt(), theWeight);
        }
        
        // -- Gen Leptons & Rec Leptons pt match (deltaR < 0.01) -- //
        reconstructionLepCompatibility(genElectrons_, electrons, "REC_el_pt_compatibility_mode1", "Reconstructed electron pair e+e- pt compatibility (deltaR < 0.01)");
        reconstructionLepCompatibility(genMuons_, muons, "REC_mu_pt_compatibility_mode1", "Reconstructed muon mu pt compatibility (deltaR < 0.01)");
        
        // -- Gen Muon & Rec Muon charge compatibility -- //
        if(genMuons_->at(0).id() == muons->at(0).id()){
          foreach(const phys::Particle &genMu, *genMuons_){
            theHistograms->fill("REC_mu_charge_compatibility_mode1", "Reconstructed muon mu charge compatibility", 75, 0., 150., genMu.pt(), theWeight);
          }
        }
        
        // -- Transverse mass of the 4 REC leptons e+e-mu nu -- //
        double recFourLepTransverseMass = (recEl1.p4() + recEl2.p4() + recMu.p4() + recNu.p4()).Mt();
        theHistograms->fill("REC_four_lep_transv_mass_mode1", "Rec Four Leptons e+e-mu nu transverse mass (mode 1)", 250, 0., 500., recFourLepTransverseMass, theWeight);
        
        // -- Transverse mass of the 4 REC leptons e+e-mu nu vs Transverse mass of the 4 GEN leptons e+e-mu nu -- //
        theHistograms->fill("REC_four_lep_transv_mass_GEN_four_lep_transv_mass_mode1", "Rec Four Leptons transverse mass vs Gen Four Leptons transverse mass e+e-mu nu (mode 1)", 250, 0., 500., 250, 0., 500., genFourLepTransverseMass, recFourLepTransverseMass, theWeight);
        
        // -- Transverse mass of the 4 REC leptons e+e-mu nu vs Invariant mass of the 4 GEN leptons e+e-mu nu -- //
        theHistograms->fill("REC_four_lep_transv_mass_GEN_four_lep_inv_mass_mode1", "Rec Four Leptons transverse mass vs Gen Four Leptons invariant mass e+e-mu nu (mode 1)", 250, 0., 500., 250, 0., 500., genFourLepInvMass, recFourLepTransverseMass, theWeight);
        
        
        
        
      }
      
      
      
      
      
  
      
      
      
      
      
      
      
      
    }
    
    else if(genMuons_->size()==2 && genElectrons_->size()==1){    // MODE 2 mu+ mu- e nu
      phys::Particle genMu1 = genMuons_->at(0);
      phys::Particle genMu2 = genMuons_->at(1);
      phys::Particle genEl = genElectrons_->at(0);
      phys::Particle genNu = genNeutrinos_->at(0);
      
      // -- Gen leptons pt -- //
      foreach(const phys::Particle &genEl, *genElectrons_){
        theHistograms->fill("GEN_el_pt_mode2", "Gen Electron pt (mu+mu- e);pt"  , 75 , 0, 150.,  genEl.pt(), theWeight);
      }
      foreach(const phys::Particle &genMu, *genMuons_){
        theHistograms->fill("GEN_mu_pt_mode2", "Gen Muons pt (mu+mu- e);pt"  , 75 , 0, 150.,  genMu.pt(), theWeight);
      }
      
      // -- Invariant mass of the mu+mu- pair -- //
      double genMuPairInvMass = (genMu1.p4() + genMu2.p4()).M();
      //double genMuPairInvMass = sqrt( pow(mu1.p4().Energy() + mu2.p4().Energy(),2) - pow(mu1.p4().Px() + mu2.p4().Px(),2) - pow(mu1.p4().Py() + mu2.p4().Py(),2) - pow(mu1.p4().Pz() + mu2.p4().Pz(),2) );
      theHistograms->fill("GEN_mu_pair_invariant_mass_mode2", "Gen Muon pair mu+mu- invariant mass (mode 2)", 75, 0, 150., genMuPairInvMass, theWeight);
      
      // -- Invariant mass of the 4 leptons mu+ mu- e nu -- // 
      double genFourLepInvMass = (genMu1.p4() + genMu2.p4() + genEl.p4() + genNu.p4()).M();
      //double genFourLepInvMass = sqrt( pow(mu1.p4().Energy() + mu2.p4().Energy() + el.p4().Energy() + nu.p4().Energy(),2) - pow(mu1.p4().Px() + mu2.p4().Px() + el.p4().Px() + nu.p4().Px(),2) - pow(mu1.p4().Py() + mu2.p4().Py() + el.p4().Py() + nu.p4().Py(),2) - pow(mu1.p4().Pz() + mu2.p4().Pz() + el.p4().Pz() + nu.p4().Pz(),2) );
      theHistograms->fill("GEN_four_leptons_invariant_mass_mode2", "Gen Leptons mu+mu-e nu invariant mass (mode 2)", 150, 0, 300., genFourLepInvMass, theWeight);      
      
      
      // -- Transverse mass of the 4 leptons mu+ mu- e nu -- //
      double genFourLepTransverseMass = (genMu1.p4() + genMu2.p4() + genEl.p4() + genNu.p4()).Mt();;
      theHistograms->fill("GEN_four_leptons_transverse_mass_mode2", "Gen Four Leptons mu+mu-e nu transverse mass (mode 2)", 150, 0., 300., genFourLepTransverseMass, theWeight);
      
      // -- Invariant mass of the muon pair mu+mu- vs Invariant mass of the 4 leptons mu+mu-e nu -- //
      theHistograms->fill("GEN_mu_pair_four_lep_inv_mass_mode2", "Gen Muon pair vs Four Leptons invariant mass (mode 2)", 150, 0., 300., 150, 0., 300., genFourLepInvMass, genMuPairInvMass, theWeight);
      
      // -- Gen Neutrino pt vs Invariant mass of the 4 leptons mu+mu-e nu -- //
      theHistograms->fill("GEN_nu_pt_four_lep_inv_mass_mode2", "Gen Neutrino pt vs Four Leptons invariant mass (mode 2)", 150, 0., 300., 150, 0., 300., genFourLepInvMass, genNu.pt(), theWeight);
      
      // -- Gen Electron pt vs Invariant mass of the 4 leptons mu+mu-e nu -- //
      theHistograms->fill("GEN_el_pt_four_lep_inv_mass_mode2", "Gen Electron pt vs Four Leptons invariant mass (mode 2)", 250, 0., 500., 150, 0., 300., genFourLepInvMass, genEl.pt(), theWeight);
      
      // -- Transverse mass of the Gen Lepton pair el nu vs Invariant mass of the 4 Gen Leptons mu+mu-el nu -- //
      double genLepPairTransverseMass = (genEl.p4() + genNu.p4()).Mt();
      theHistograms->fill("GEN_lep_pair_transv_mass_four_lep_inv_mass_mode2", "Gen Lepton pair el nu transverse mass vs Four Leptons invariant mass (mode 2)", 250, 0., 500., 75, 0., 150., genFourLepInvMass, genLepPairTransverseMass, theWeight);
      
      // -- Transverse vs Invariant mass of the 4 gen leptons mu+mu-e nu -- // 
      theHistograms->fill("GEN_four_lep_transv_mass_four_lep_inv_mass_mode2", "Gen Four Lepton transverse vs invariant mass (mode 2); Inv mass; Transv mass", 100, 0., 500., 100, 0., 500., genFourLepInvMass, genFourLepTransverseMass, theWeight);
      
      
      // -- RECONSTRUCTED LEPTONS mu+ mu- e -- //
      if(muons->size()==2 && electrons->size()==1){
        phys::Particle recMu1 = muons->at(0);
        phys::Particle recMu2 = muons->at(1);
        phys::Particle recEl = electrons->at(0);
        phys::Particle recNu = *met;
      
        // -- Gen Leptons & Rec Leptons match -- //
        foreach(const phys::Particle &genEl, *genElectrons_){
          theHistograms->fill("REC_el_compatibility_mode2", "Reconstructed electron e compatibility", 75, 0, 150.,  genEl.pt(), theWeight);
        }
        foreach(const phys::Particle &genMu, *genMuons_){
          theHistograms->fill("REC_mu_compatibility_mode2", "Reconstructed muon pair mu+mu- compatibility", 75, 0, 150.,  genMu.pt(), theWeight);
        }
        
        // -- Gen Leptons & Rec Leptons pt match (deltaR < 0.01) -- //
        reconstructionLepCompatibility(genElectrons_, electrons, "REC_el_pt_compatibility_mode2", "Reconstructed electron e compatibility (deltaR < 0.01)");
        reconstructionLepCompatibility(genMuons_, muons, "REC_mu_pt_compatibility_mode2", "Reconstructed muon pair mu+mu- compatibility (deltaR < 0.01)");
        
        // -- Gen Electron & Rec Electron charge compatibility -- //
        if(genElectrons_->at(0).id() == electrons->at(0).id()){
          foreach(const phys::Particle &genEl, *genElectrons_){
            theHistograms->fill("REC_el_charge_compatibility_mode2", "Reconstructed electron e charge compatibility", 75, 0., 150., genEl.pt(), theWeight);
          }
        }
        
        // -- Transverse mass of the 4 REC leptons mu+mu-e nu -- //
        double recFourLepTransverseMass = (recMu1.p4() + recMu2.p4() + recEl.p4() + recNu.p4()).Mt();
        theHistograms->fill("REC_four_lep_transv_mass_mode2", "Rec Four Leptons mu+mu-e nu transverse mass (mode 2)", 250, 0., 500., recFourLepTransverseMass, theWeight);
        
        // -- Transverse mass of the 4 REC leptons mu+mu-e nu vs Transverse mass of the 4 GEN leptons mu+mu-e nu -- //
        theHistograms->fill("REC_four_lep_transv_mass_GEN_four_lep_transv_mass_mode2", "Rec Four Leptons transverse mass vs Gen Four Leptons transverse mass mu+mu-e nu (mode 2)", 250, 0., 500., 250, 0., 500., genFourLepTransverseMass, recFourLepTransverseMass, theWeight);
        
        // -- Transverse mass of the 4 REC leptons e+e-mu nu vs Invariant mass of the 4 GEN leptons mu+mu-e nu -- //
        theHistograms->fill("REC_four_lep_transv_mass_GEN_four_lep_inv_mass_mode2", "Rec Four Leptons transverse mass vs Gen Four Leptons invariant mass mu+mu-e nu (mode 2)", 250, 0., 500., 250, 0., 500., genFourLepInvMass, recFourLepTransverseMass, theWeight);
        
        // -- REC nu (met) pt vs Invariant mass of the 4 GEN leptons mu+mu-e nu -- //
        theHistograms->fill("REC_nu_pt_GEN_four_lep_inv_mass_mode2", "Rec nu (met) pt vs Gen Four Leptons invariant mass mu+mu-e nu (mode 2)", 150, 0., 300., 250, 0., 500., recNu.pt(), genFourLepInvMass, theWeight);
        
        // -- REC nu (met) pt vs Transverse mass of the 4 GEN leptons mu+mu-e nu -- //
        theHistograms->fill("REC_nu_pt_GEN_four_lep_transv_mass_mode2", "Rec nu (met) pt vs Gen Four Leptons transverse mass mu+mu-e nu (mode 2)", 150, 0., 300., 250, 0., 500., recNu.pt(), genFourLepTransverseMass, theWeight);
        
        // -- REC nu (met) pt vs Transverse mass of the 4 REC leptons mu+mu-e nu -- //
        theHistograms->fill("REC_nu_pt_REC_four_lep_transv_mass_mode2", "Rec nu (met) pt vs Rec Four Leptons transverse mass mu+mu-e nu (mode 2)", 150, 0., 300., 250, 0., 500., recNu.pt(), recFourLepTransverseMass, theWeight);
        
        
        
        
      
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
    
}


void WlllnuAnalyzer::genEventSetup(){
  
  genChLeptons_->clear();
  genElectrons_->clear();
  genMuons_->clear();
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
  
  std::sort(genChLeptons_  ->begin(), genChLeptons_  ->end(), [](const Particle& a, const Particle& b){ return a.pt() < b.pt(); });
  std::sort(genElectrons_  ->begin(), genElectrons_  ->end(), [](const Particle& a, const Particle& b){ return a.pt() < b.pt(); });
  std::sort(genMuons_      ->begin(), genMuons_      ->end(), [](const Particle& a, const Particle& b){ return a.pt() < b.pt(); });
  std::sort(genNeutrinos_  ->begin(), genNeutrinos_  ->end(), [](const Particle& a, const Particle& b){ return a.pt() < b.pt(); });
  
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

/*
void histogramFill(std::vector<phys::Particle>* genLep, string histName, string histTitle, int nBin, double start, double end, string var, double weight){  
  
  foreach(const phys::Particle &genL, *genLep){
    theHistograms->fill(histName, histTitle, nBin, start, end, genL.var, weight);
  }
  
  
}  
*/  
  
