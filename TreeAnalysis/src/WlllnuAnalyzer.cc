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
  
  
  
  // -- W boson 3 leptons decay possible configurations -- //
  if(genChLeptons_->size()==3){
    
    if(genElectrons_->size()==2 && genMuons_->size()==1 && (genElectrons_->at(0).id() + genElectrons_->at(1).id())==0){  // MODE 1 e+ e- mu nu
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
    
      // -- Invariant mass of the gen e+e- pair -- //
      double genElPairInvMass = (genEl1.p4() + genEl2.p4()).M();
      theHistograms->fill("GEN_el_pair_invariant_mass_mode1", "Gen Electron pair e+e- invariant mass (mode 1)", 75, 0, 150., genElPairInvMass, theWeight);
      
      // -- Invariant mass of the 4 gen leptons e+ e- mu nu -- // 
      double genFourLepInvMass = (genEl1.p4() + genEl2.p4() + genMu.p4() + genNu.p4()).M();
      theHistograms->fill("GEN_four_leptons_invariant_mass_mode1", "Gen Leptons e+e-mu nu invariant mass (mode 1)", 150, 0, 300., genFourLepInvMass, theWeight);
      
      // -- Transverse mass of the 4 gen leptons e+ e- mu nu -- //
      double genFourLepTransverseMass = (genEl1.p4() + genEl2.p4() + genMu.p4() + genNu.p4()).Mt();;
      theHistograms->fill("GEN_four_leptons_transverse_mass_mode1", "Gen Four Leptons e+e-mu nu transverse mass (mode 1)", 100, 0., 500., genFourLepTransverseMass, theWeight);
      
      // -- Invariant mass of the gen electron pair e+e- vs Invariant mass of the 4 gen leptons e+e-mu nu -- //
      theHistograms->fill("GEN_el_pair_four_lep_inv_mass_mode1", "Gen Electron pair vs GEN Four Leptons invariant mass (mode 1);4l mass;e+e- mass", 100, 0., 500., 60, 0., 300., genFourLepInvMass, genElPairInvMass, theWeight);
      
      // -- Gen Neutrino pt vs Invariant mass of the 4 gen leptons e+e-mu nu -- //
      theHistograms->fill("GEN_nu_pt_four_lep_inv_mass_mode1", "Gen Neutrino pt vs GEN Four Leptons invariant mass (mode 1);4l mass;nu pt", 100, 0., 500., 60, 0., 300., genFourLepInvMass, genNu.pt(), theWeight);
      
      // -- Gen Muon pt vs Invariant mass of the 4 gen leptons e+e-mu nu -- //
      theHistograms->fill("GEN_mu_pt_four_lep_inv_mass_mode1", "Gen Muon pt vs GEN Four Leptons invariant mass (mode 1);4l mass;mu pt", 100, 0., 500., 60, 0., 300., genFourLepInvMass, genMu.pt(), theWeight);
      
      // -- Transverse mass of the gen lepton pair mu nu vs Invariant mass of the 4 gen leptons e+e-mu nu -- //
      double genLepPairTransverseMass = (genMu.p4() + genNu.p4()).Mt();
      theHistograms->fill("GEN_lep_pair_transv_mass_mode1", "Gen Lepton pair mu nu transverse mass (mode 1)", 75, 0., 150., genLepPairTransverseMass, theWeight);
      theHistograms->fill("GEN_lep_pair_transv_mass_four_lep_inv_mass_mode1", "Gen Lepton pair mu nu transverse mass vs Four Leptons invariant mass (mode 1);4l mass;mu nu mass", 250, 0., 500., 30, 0., 150., genFourLepInvMass, genLepPairTransverseMass, theWeight);
      
      // -- Transverse vs Invariant mass of the 4 gen leptons e+e-mu nu -- // 
      theHistograms->fill("GEN_four_lep_transv_mass_four_lep_inv_mass_mode1", "Gen Four Lepton transverse vs invariant mass (mode 1); Inv mass; Transv mass", 100, 0., 500., 100, 0., 500., genFourLepInvMass, genFourLepTransverseMass, theWeight);
      
      
      
      // ------------------------------ GEN LEVEL SIGNAL DEFINITION ----------------------------- //
      if(isGen_mode1(genFourLepInvMass, genElPairInvMass, genLepPairTransverseMass, genMu.pt(), genNu.pt()) == true){
        theHistograms->fill("GEN_signal_mode1", "Gen signal (mode 1)", 3, 0., 3., 1, theWeight);
      }
      else{
        theHistograms->fill("not_GEN_signal_mode1", "NOT Gen signal (mode 1)", 3, 0., 3., 1, theWeight);
      }
      
      
      
        
        
        
        
        
      
      
      
      
      // -- RECONSTRUCTED LEPTONS e+ e- mu nu -- //
      if(electrons->size()==2 && muons->size()==1 && (electrons->at(0).id() + electrons->at(1).id())==0){
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
        
        
        // -- Transverse mass of the 4 REC leptons e+e-mu nu -- //
        double recFourLepTransverseMass = (recEl1.p4() + recEl2.p4() + recMu.p4() + recNu.p4()).Mt();
        theHistograms->fill("REC_four_lep_transv_mass_mode1", "Rec Four Leptons e+e-mu nu transverse mass (mode 1)", 250, 0., 500., recFourLepTransverseMass, theWeight);
        
        // -- Transverse mass of the 4 REC leptons e+e-mu nu vs Transverse mass of the 4 GEN leptons e+e-mu nu -- //
        theHistograms->fill("REC_four_lep_transv_mass_GEN_four_lep_transv_mass_mode1", "Rec Four Leptons transverse mass vs Gen Four Leptons transverse mass e+e-mu nu (mode 1); Gen 4l;Rec 4l", 100, 0., 500., 100, 0., 500., genFourLepTransverseMass, recFourLepTransverseMass, theWeight);
        
        
        // -- Transverse mass of the 4 REC leptons e+e-mu nu vs Invariant mass of the 4 GEN leptons e+e-mu nu -- //
        theHistograms->fill("REC_four_lep_transv_mass_GEN_four_lep_inv_mass_mode1", "Rec Four Leptons transverse mass vs Gen Four Leptons invariant mass e+e-mu nu (mode 1);Gen 4l;Rec 4l", 100, 0., 500., 100, 0., 500., genFourLepInvMass, recFourLepTransverseMass, theWeight);
        
        // -- Invariant mass of the rec electron pair e+e- vs Transverse mass of the 4 rec leptons e+e-mu nu -- //
        double recElPairInvMass = (recEl1.p4() + recEl2.p4()).M();
        theHistograms->fill("REC_el_pair_invariant_mass_mode1", "Rec Electron pair e+e- invariant mass (mode 1)", 75, 0, 150., recElPairInvMass, theWeight);
        theHistograms->fill("REC_el_pair_inv_mass_REC_four_lep_transv_mass_mode1", "Rec Electron pair invariant mass vs Rec Four Leptons transverse mass (mode 1);4l mass;e+e- mass", 100, 0., 500., 60, 0., 300., recFourLepTransverseMass, recElPairInvMass, theWeight);
        
        // -- Transverse mass of the rec lepton pair mu nu vs Invariant mass of the 4 rec leptons e+e-mu nu -- //
        double recLepPairTransverseMass = (recMu.p4() + recNu.p4()).Mt();
        theHistograms->fill("REC_lep_pair_transv_mass_mode1", "Rec Lepton pair mu nu transverse mass (mode 1)", 75, 0., 150., recLepPairTransverseMass, theWeight);
        theHistograms->fill("REC_lep_pair_transv_mass_REC_four_lep_transv_mass_mode1", "Rec Lepton pair mu nu transverse mass vs Rec Four Leptons transverse mass (mode 1);4l mass;mu nu mass", 100, 0., 500., 60, 0., 300., recFourLepTransverseMass, recLepPairTransverseMass, theWeight);
        
        
        
        
        // -- REC nu (met) pt vs Transverse mass of the 4 REC leptons e+e-mu nu -- //
        theHistograms->fill("REC_nu_pt_REC_four_lep_transv_mass_mode1", "Rec nu (met) pt vs Rec Four Leptons transverse mass e+e-mu nu (mode 1);4l mass;nu pt", 100, 0., 500., 60, 0., 300., recFourLepTransverseMass, recNu.pt(), theWeight);
        
        // -- REC mu pt vs Transverse mass of the 4 REC leptons e+e-mu nu -- //
        theHistograms->fill("REC_mu_pt_REC_four_lep_transv_mass_mode1", "Rec mu pt vs Rec Four Leptons transverse mass e+e-mu nu (mode 1); 4l mass;mu pt", 100, 0., 500., 30, 0., 150., recFourLepTransverseMass, recMu.pt(), theWeight);
        
        
        // ------------------------------ REC LEVEL SIGNAL DEFINITION ----------------------------- //
        if(isGen_mode1(genFourLepInvMass, genElPairInvMass, genLepPairTransverseMass, genMu.pt(), genNu.pt()) == true && isRec_mode1(recFourLepTransverseMass, recElPairInvMass, recLepPairTransverseMass, recMu.pt(), recNu.pt()) == true){
          // -- REC Z0 & W bosons pt -- //
          theHistograms->fill("REC_Z_pt_mode1", "Rec Z boson pt (mode 1);pt", 250, 0, 500., ZW->first().pt(), theWeight);
          theHistograms->fill("REC_W_pt_mode1", "Rec W boson pt (mode 1);pt", 250, 0, 500., ZW->second().pt(), theWeight);
          // deltaR between REC Z0 & W bosons -- //
          theHistograms->fill("REC_ZW_bosons_deltaR_mode1", "Rec ZW Bosons deltaR (mode 1)", 20, 0., 5., deltaR(ZW->first(),ZW->second()), theWeight);
          if(deltaR(ZW->first(),ZW->second())<0.1){
            theHistograms->fill("REC_ZW_bosons_deltaR_zoom_mode1", "Rec ZW Bosons deltaR zoom (mode 1)", 10, 0., 0.1, deltaR(ZW->first(),ZW->second()), theWeight);
          }
          
          // -- Signal Efficiency -- //       
          theHistograms->fill("GEN_REC_signal_mode1", "Signal Efficiency (mode 1)", 3, 0., 3., 1, theWeight);
        }
        else if(isGen_mode1(genFourLepInvMass, genElPairInvMass, genLepPairTransverseMass, genMu.pt(), genNu.pt()) == false && isRec_mode1(recFourLepTransverseMass, recElPairInvMass, recLepPairTransverseMass, recMu.pt(), recNu.pt()) == true){
          
          // -- Background Efficiency -- //
          theHistograms->fill("not_GEN_REC_signal_mode1", "Background Efficiency (mode 1)", 3, 0., 3., 1, theWeight);
        }
        
        
        }
        
        
        
        
      }
      
      
      
      
      
  
      
      
      
      
      
      
      
      
    }
    
    else if(genMuons_->size()==2 && genElectrons_->size()==1 && (genMuons_->at(0).id() + genMuons_->at(1).id())==0){    // MODE 2 mu+ mu- e nu
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
      theHistograms->fill("GEN_mu_pair_invariant_mass_mode2", "Gen Muon pair mu+mu- invariant mass (mode 2)", 75, 0, 150., genMuPairInvMass, theWeight);
      
      // -- Invariant mass of the 4 leptons mu+ mu- e nu -- // 
      double genFourLepInvMass = (genMu1.p4() + genMu2.p4() + genEl.p4() + genNu.p4()).M();
      theHistograms->fill("GEN_four_leptons_invariant_mass_mode2", "Gen Leptons mu+mu-e nu invariant mass (mode 2)", 150, 0, 300., genFourLepInvMass, theWeight);      
      
      // -- Transverse mass of the 4 leptons mu+ mu- e nu -- //
      double genFourLepTransverseMass = (genMu1.p4() + genMu2.p4() + genEl.p4() + genNu.p4()).Mt();;
      theHistograms->fill("GEN_four_leptons_transverse_mass_mode2", "Gen Four Leptons mu+mu-e nu transverse mass (mode 2)", 150, 0., 300., genFourLepTransverseMass, theWeight);
      
      // -- Invariant mass of the muon pair mu+mu- vs Invariant mass of the 4 leptons mu+mu-e nu -- //
      theHistograms->fill("GEN_mu_pair_four_lep_inv_mass_mode2", "Gen Muon pair vs Four Leptons invariant mass (mode 2);4l mass;mu+mu- mass", 100, 0., 500., 60, 0., 300., genFourLepInvMass, genMuPairInvMass, theWeight);
      
      // -- Gen Neutrino pt vs Invariant mass of the 4 leptons mu+mu-e nu -- //
      theHistograms->fill("GEN_nu_pt_four_lep_inv_mass_mode2", "Gen Neutrino pt vs Four Leptons invariant mass (mode 2);4l mass;nu pt", 100, 0., 500., 60, 0., 300., genFourLepInvMass, genNu.pt(), theWeight);
      
      // -- Gen Electron pt vs Invariant mass of the 4 leptons mu+mu-e nu -- //
      theHistograms->fill("GEN_el_pt_four_lep_inv_mass_mode2", "Gen Electron pt vs Four Leptons invariant mass (mode 2);4l mass;nu pt", 100, 0., 500., 60, 0., 300., genFourLepInvMass, genEl.pt(), theWeight);
      
      // -- Transverse mass of the Gen Lepton pair el nu vs Invariant mass of the 4 Gen Leptons mu+mu-el nu -- //
      double genLepPairTransverseMass = (genEl.p4() + genNu.p4()).Mt();
      theHistograms->fill("GEN_lep_pair_transv_mass_mode2", "Gen Lepton pair e nu transverse mass (mode 2)", 75, 0., 150., genLepPairTransverseMass, theWeight);
      theHistograms->fill("GEN_lep_pair_transv_mass_four_lep_inv_mass_mode2", "Gen Lepton pair e nu transverse mass vs Four Leptons invariant mass (mode 2);4l mass;e nu mass", 100, 0., 500., 75, 0., 150., genFourLepInvMass, genLepPairTransverseMass, theWeight);
      
      // -- Transverse vs Invariant mass of the 4 gen leptons mu+mu-e nu -- // 
      theHistograms->fill("GEN_four_lep_transv_mass_four_lep_inv_mass_mode2", "Gen Four Lepton transverse vs invariant mass (mode 2); Inv mass; Transv mass", 100, 0., 500., 100, 0., 500., genFourLepInvMass, genFourLepTransverseMass, theWeight);
      
      // ------------------------------ GEN LEVEL SIGNAL DEFINITION ----------------------------- //
      if(isGen_mode2(genFourLepInvMass, genMuPairInvMass, genLepPairTransverseMass, genEl.pt(), genNu.pt()) == true){
        theHistograms->fill("GEN_signal_mode2", "Gen signal (mode 2)", 3, 0., 3., 1, theWeight);
      }
      else{
        theHistograms->fill("not_GEN_signal_mode2", "NOT Gen signal (mode 2)", 3, 0., 3., 1, theWeight);
      }
      
      
      
      
      
      
      
      
      // -- RECONSTRUCTED LEPTONS mu+ mu- e nu -- //
      if(muons->size()==2 && electrons->size()==1 && (muons->at(0).id() + muons->at(1).id())==0){
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
        
        
        // -- Transverse mass of the 4 REC leptons mu+mu-e nu -- //
        double recFourLepTransverseMass = (recMu1.p4() + recMu2.p4() + recEl.p4() + recNu.p4()).Mt();
        theHistograms->fill("REC_four_lep_transv_mass_mode2", "Rec Four Leptons mu+mu-e nu transverse mass (mode 2)", 250, 0., 500., recFourLepTransverseMass, theWeight);
        
        // -- Transverse mass of the 4 REC leptons mu+mu-e nu vs Transverse mass of the 4 GEN leptons mu+mu-e nu -- //
        theHistograms->fill("REC_four_lep_transv_mass_GEN_four_lep_transv_mass_mode2", "Rec Four Leptons transverse mass vs Gen Four Leptons transverse mass mu+mu-e nu (mode 2);Gen 4l;Rec 4l", 100, 0., 500., 100, 0., 500., genFourLepTransverseMass, recFourLepTransverseMass, theWeight);
        
        // -- Transverse mass of the 4 REC leptons e+e-mu nu vs Invariant mass of the 4 GEN leptons mu+mu-e nu -- //
        theHistograms->fill("REC_four_lep_transv_mass_GEN_four_lep_inv_mass_mode2", "Rec Four Leptons transverse mass vs Gen Four Leptons invariant mass mu+mu-e nu (mode 2);Gen 4l;Rec 4l", 100, 0., 500., 100, 0., 500., genFourLepInvMass, recFourLepTransverseMass, theWeight);
        
        // -- Invariant mass of the rec muon pair mu+mu- vs Transverse mass of the 4 rec leptons mu+mu-e nu -- //
        double recMuPairInvMass = (recMu1.p4() + recMu2.p4()).M();
        theHistograms->fill("REC_mu_pair_invariant_mass_mode2", "Rec Muon pair mu+mu- invariant mass (mode 2)", 75, 0, 150., recMuPairInvMass, theWeight);
        theHistograms->fill("REC_mu_pair_inv_mass_REC_four_lep_transv_mass_mode2", "Rec Muon pair invariant mass vs Rec Four Leptons transverse mass (mode 2);4l mass;mu+mu- mass", 100, 0., 500., 60, 0., 300., recFourLepTransverseMass, recMuPairInvMass, theWeight);
        
        // -- Transverse mass of the rec lepton pair e nu vs Transverse mass of the 4 rec leptons mu+mu-e nu -- //
        double recLepPairTransverseMass = (recEl.p4() + recNu.p4()).Mt();
        theHistograms->fill("REC_lep_pair_transv_mass_mode2", "Rec Lepton pair e nu transverse mass (mode 2)", 75, 0., 150., recLepPairTransverseMass, theWeight);
        theHistograms->fill("REC_lep_pair_transv_mass_REC_four_lep_transv_mass_mode2", "Rec Lepton pair e nu transverse mass vs Rec Four Leptons transverse mass (mode 2);4l mass;e nu mass", 100, 0., 500., 60, 0., 300., recFourLepTransverseMass, recLepPairTransverseMass, theWeight);
        
        
        
        // -- REC nu (met) pt vs Transverse mass of the 4 REC leptons mu+mu-e nu -- //
        theHistograms->fill("REC_nu_pt_REC_four_lep_transv_mass_mode2", "Rec nu (met) pt vs Rec Four Leptons transverse mass mu+mu-e nu (mode 2);4l mass;nu pt", 100, 0., 500., 60, 0., 300., recFourLepTransverseMass, recNu.pt(), theWeight);
        
        // -- REC el pt vs Transverse mass of the 4 REC leptons mu+mu-e nu -- //
        theHistograms->fill("REC_el_pt_REC_four_lep_transv_mass_mode2", "Rec el pt vs Rec Four Leptons transverse mass mu+mu-e nu (mode 2);4l mass;e pt", 100, 0., 500., 60, 0., 300., recFourLepTransverseMass, recEl.pt(), theWeight);
        
        
        
        // ------------------------------ REC LEVEL SIGNAL DEFINITION ----------------------------- //
        if(isGen_mode2(genFourLepInvMass, genMuPairInvMass, genLepPairTransverseMass, genEl.pt(), genNu.pt()) == true && isRec_mode2(recFourLepTransverseMass, recMuPairInvMass, recLepPairTransverseMass, recEl.pt(), recNu.pt()) == true){
          // -- REC Z0 & W bosons pt -- // 
          theHistograms->fill("REC_Z_pt_mode2", "Rec Z boson pt (mode 2);pt", 250, 0, 500., ZW->first().pt(), theWeight);
          theHistograms->fill("REC_W_pt_mode2", "Rec W boson pt (mode 2);pt", 250, 0, 500., ZW->second().pt(), theWeight);
          // -- deltaR between REC Z0 & W bosons -- // 
          theHistograms->fill("REC_ZW_bosons_deltaR_mode2", "Rec ZW Bosons deltaR (mode 2)", 20, 0., 5., deltaR(ZW->first(),ZW->second()), theWeight);
          if(deltaR(ZW->first(),ZW->second())<0.1){
            theHistograms->fill("REC_ZW_bosons_deltaR_zoom_mode2", "Rec ZW Bosons deltaR zoom (mode 2)", 10, 0., 0.1, deltaR(ZW->first(),ZW->second()), theWeight);
          }
          
          // -- Signal Efficiency -- //
          theHistograms->fill("GEN_REC_signal_mode2", "Signal Efficiency (mode 2)", 3, 0., 3., 1, theWeight);
        }
        else if(isGen_mode2(genFourLepInvMass, genMuPairInvMass, genLepPairTransverseMass, genEl.pt(), genNu.pt()) == false && isRec_mode2(recFourLepTransverseMass, recMuPairInvMass, recLepPairTransverseMass, recEl.pt(), recNu.pt()) == true){
          
          // -- Background Efficiency -- //
          theHistograms->fill("not_GEN_REC_signal_mode2", "Background Efficiency (mode 2)", 3, 0., 3., 1, theWeight);
        }
        
        
        }
      
      } 
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
    }
    
    else if(genElectrons_->size()==3 && genMuons_->size()==0){      // mode 3 e+ e- e nu
    
      // --------------------- GEN LEVEL --------------------- //
      
      // Gen W --> e nu
      for(auto l : *genElectrons_){
        for(auto v : *genNeutrinos_){
	  if( abs(l.id() + v.id()) == 1 ){
	    Boson<Particle> Wcand(l,v);
	    if(GenWtoLNuDefinition(Wcand)){
	      genWlepCandidates_->push_back(Wcand);
            }
	  }
        }
      }
      // Gen Z --> e+e-    
      for(size_t i=0; i<genElectrons_->size(); i++){
        Particle& l1 = genElectrons_->at(i);
        for(size_t j = i+1; j < genElectrons_->size(); j++){
	  Particle& l2 = genElectrons_->at(j);
	  if( l1.id() + l2.id() == 0 ){
	    Boson<Particle> Zcand(l1,l2);
	    if(ZBosonDefinition(Zcand)){
	      genZlepCandidates_->push_back(Zcand);
            }
	  }
        }
      }
      // Gen ZW --> e+e-e nu
      if(genZlepCandidates_->size() >= 1 && genWlepCandidates_->size() >= 1){	
        std::sort(genZlepCandidates_->begin(), genZlepCandidates_->end(), MassComparator(phys::ZMASS));
        Boson<Particle>& genZ0 = genZlepCandidates_->front();

        std::sort(genWlepCandidates_->begin(), genWlepCandidates_->end(), MassComparator(phys::WMASS));
        auto genW = std::find_if(genWlepCandidates_->begin(), genWlepCandidates_->end(),
			     [&genZ0](auto W){ return !haveCommonDaughter(genZ0, W); }
			     );
        if(genW != genWlepCandidates_->end()){
          genZW_ = DiBoson<Particle, Particle>(genZ0, *genW);
        }
      }
      
      // -- Gen Leptons definition -- //
      phys::Particle genEl1 = genZW_.first().daughter(0);
      phys::Particle genEl2 = genZW_.first().daughter(1);
      phys::Particle genEl3 = genZW_.second().daughter(0);
      phys::Particle genNu = genZW_.second().daughter(1);
      
      
      // -- Invariant mass of the gen e+e- pair -- //
      double genElPairInvMass = (genEl1.p4() + genEl2.p4()).M();
      theHistograms->fill("GEN_el_pair_invariant_mass_mode3", "Gen Electron pair e+e- invariant mass (mode 3)", 75, 0, 150., genElPairInvMass, theWeight);
      
      // -- Invariant mass of the 4 gen leptons e+ e- e nu -- // 
      double genFourLepInvMass = (genEl1.p4() + genEl2.p4() + genEl3.p4() + genNu.p4()).M();
      theHistograms->fill("GEN_four_leptons_invariant_mass_mode3", "Gen Leptons e+e-e nu invariant mass (mode 3)", 150, 0, 300., genFourLepInvMass, theWeight);
      
      // -- Transverse mass of the 4 gen leptons e+ e- e nu -- //
      double genFourLepTransverseMass = (genEl1.p4() + genEl2.p4() + genEl3.p4() + genNu.p4()).Mt();;
      theHistograms->fill("GEN_four_leptons_transverse_mass_mode3", "Gen Four Leptons e+e-e nu transverse mass (mode 3)", 100, 0., 500., genFourLepTransverseMass, theWeight);
      
      // -- Invariant mass of the gen electron pair e+e- vs Invariant mass of the 4 gen leptons e+e-e nu -- //
      theHistograms->fill("GEN_el_pair_four_lep_inv_mass_mode3", "Gen Electron pair vs GEN Four Leptons invariant mass (mode 3);4l mass;e+e- mass", 100, 0., 500., 60, 0., 300., genFourLepInvMass, genElPairInvMass, theWeight);
      
      // -- Gen Neutrino pt vs Invariant mass of the 4 gen leptons e+e-e nu -- //
      theHistograms->fill("GEN_nu_pt_four_lep_inv_mass_mode3", "Gen Neutrino pt vs GEN Four Leptons invariant mass (mode 3);4l mass;nu pt", 100, 0., 500., 60, 0., 300., genFourLepInvMass, genNu.pt(), theWeight);
      
      // -- Gen Electron pt vs Invariant mass of the 4 gen leptons e+e-e nu -- //
      theHistograms->fill("GEN_el_pt_four_lep_inv_mass_mode3", "Gen Electron pt vs GEN Four Leptons invariant mass (mode 3);4l mass;mu pt", 100, 0., 500., 60, 0., 300., genFourLepInvMass, genEl3.pt(), theWeight);
      
      // -- Transverse mass of the gen lepton pair e nu vs Invariant mass of the 4 gen leptons e+e-e nu -- //
      double genLepPairTransverseMass = (genEl3.p4() + genNu.p4()).Mt();
      theHistograms->fill("GEN_lep_pair_transv_mass_mode3", "Gen Lepton pair e nu transverse mass (mode 3)", 75, 0., 150., genLepPairTransverseMass, theWeight);
      theHistograms->fill("GEN_lep_pair_transv_mass_four_lep_inv_mass_mode3", "Gen Lepton pair mu nu transverse mass vs Four Leptons invariant mass (mode 3);4l mass;e nu mass", 250, 0., 500., 30, 0., 150., genFourLepInvMass, genLepPairTransverseMass, theWeight);
      
      // -- Transverse vs Invariant mass of the 4 gen leptons e+e-e nu -- // 
      theHistograms->fill("GEN_four_lep_transv_mass_four_lep_inv_mass_mode3", "Gen Four Lepton transverse vs invariant mass (mode 3); Inv mass; Transv mass", 100, 0., 500., 100, 0., 500., genFourLepInvMass, genFourLepTransverseMass, theWeight);
    
    
    
      
      
      
      
      
      
      
      
      // --------------------- REC LEVEL --------------------- //
      if(abs(ZW->first().daughter(0).id())==11 && abs(ZW->first().daughter(1).id())==11 && abs(ZW->second().daughter(0).id())==11){      // electrons->size()==3 && muons->size()==0
        
        // deltaR(Z0,W)
        theHistograms->fill("REC_ZW_bosons_deltaR_mode3", "Rec ZW Bosons deltaR (mode 3)", 20, 0., 5., deltaR(ZW->first(),ZW->second()), theWeight);
        if(deltaR(ZW->first(),ZW->second())<0.1){
          theHistograms->fill("REC_ZW_bosons_deltaR_zoom_mode3", "Rec ZW Bosons deltaR zoom (mode 3)", 10, 0., 0.1, deltaR(ZW->first(),ZW->second()), theWeight);
        }      
        // Z0 & W bosons pt
        theHistograms->fill("REC_Z_pt_mode3", "Rec Z boson pt (mode 3);pt", 150, 0, 300., ZW->first().pt(), theWeight);
        theHistograms->fill("REC_W_pt_mode3", "Rec W boson pt (mode 3);pt", 150, 0, 300., ZW->second().pt(), theWeight);
        // REC bosons & REC leptons definition -- //
        phys::Particle recEl1 = ZW->first().daughter(0);
        phys::Particle recEl2 = ZW->first().daughter(1);
        phys::Particle recEl3 = ZW->second().daughter(0);
        phys::Particle recNu = ZW->second().daughter(1);
        
        
        if(genEl3.id()==recEl3.id()){
        // -- Transverse mass of the 4 REC leptons e+e-e nu -- //
        double recFourLepTransverseMass = (recEl1.p4() + recEl2.p4() + recEl3.p4() + recNu.p4()).Mt();
        theHistograms->fill("REC_four_lep_transv_mass_mode3", "Rec Four Leptons mu+mu-e nu transverse mass (mode 3)", 250, 0., 500., recFourLepTransverseMass, theWeight);
        
        // -- Transverse mass of the 4 REC leptons e+e-e nu vs Transverse mass of the 4 GEN leptons e+e-e nu -- //
        theHistograms->fill("REC_four_lep_transv_mass_GEN_four_lep_transv_mass_mode3", "Rec Four Leptons transverse mass vs Gen Four Leptons transverse mass mu+mu-e nu (mode 2);Gen 4l;Rec 4l", 100, 0., 500., 100, 0., 500., genFourLepTransverseMass, recFourLepTransverseMass, theWeight);
        
        // -- Transverse mass of the 4 REC leptons e+e-mu nu vs Invariant mass of the 4 GEN leptons e+e-e nu -- //
        theHistograms->fill("REC_four_lep_transv_mass_GEN_four_lep_inv_mass_mode3", "Rec Four Leptons transverse mass vs Gen Four Leptons invariant mass mu+mu-e nu (mode 2);Gen 4l;Rec 4l", 100, 0., 500., 100, 0., 500., genFourLepInvMass, recFourLepTransverseMass, theWeight);
        
        // -- Invariant mass of the rec muon pair e+e- vs Transverse mass of the 4 rec leptons e+e-e nu -- //
        double recElPairInvMass = (recEl1.p4() + recEl2.p4()).M();
        theHistograms->fill("REC_el_pair_invariant_mass_mode3", "Rec Electron pair e+e- invariant mass (mode 3)", 75, 0, 150., recElPairInvMass, theWeight);
        theHistograms->fill("REC_el_pair_inv_mass_REC_four_lep_transv_mass_mode3", "Rec Electron pair invariant mass vs Rec Four Leptons transverse mass (mode 3);4l mass;e+e- mass", 100, 0., 500., 60, 0., 300., recFourLepTransverseMass, recElPairInvMass, theWeight);
        
        // -- Transverse mass of the rec lepton pair e nu vs Transverse mass of the 4 rec leptons e+e-e nu -- //
        double recLepPairTransverseMass = (recEl3.p4() + recNu.p4()).Mt();
        theHistograms->fill("REC_lep_pair_transv_mass_mode3", "Rec Lepton pair e nu transverse mass (mode 3)", 75, 0., 150., recLepPairTransverseMass, theWeight);
        theHistograms->fill("REC_lep_pair_transv_mass_REC_four_lep_transv_mass_mode3", "Rec Lepton pair e nu transverse mass vs Rec Four Leptons transverse mass (mode 3);4l mass;e nu mass", 100, 0., 500., 60, 0., 300., recFourLepTransverseMass, recLepPairTransverseMass, theWeight);       
        
        // -- REC nu (met) pt vs Transverse mass of the 4 REC leptons e+e-e nu -- //
        theHistograms->fill("REC_nu_pt_REC_four_lep_transv_mass_mode3", "Rec nu (met) pt vs Rec Four Leptons transverse mass e+e-e nu (mode 3);4l mass;nu pt", 100, 0., 500., 60, 0., 300., recFourLepTransverseMass, recNu.pt(), theWeight);
        
        // -- REC el (W daughter) pt vs Transverse mass of the 4 REC leptons e+e-e nu -- //
        theHistograms->fill("REC_el_pt_REC_four_lep_transv_mass_mode3", "Rec el (W daughter) pt vs Rec Four Leptons transverse mass mu+mu-e nu (mode 3);4l mass;e pt", 100, 0., 500., 60, 0., 300., recFourLepTransverseMass, recEl3.pt(), theWeight);
        
        
        
        
        
        
        
        
        
        }
      }
      
      
      
      
      
      
      
      
      
      
      
      
    }
    
    else if(genMuons_->size()==3 && genElectrons_->size()==0){          // mode 4 mu+ mu- mu nu 
      
      
      // --------------------- GEN LEVEL --------------------- //
      // Gen W --> mu nu
      for(auto l : *genMuons_){
        for(auto v : *genNeutrinos_){
	  if( abs(l.id() + v.id()) == 1 ){
	    Boson<Particle> Wcand(l,v);
	    if(GenWtoLNuDefinition(Wcand)){
	      genWlepCandidates_->push_back(Wcand);
            }
	  }
        }
      }
      // Gen Z --> e+e-    
      for(size_t i=0; i<genMuons_->size(); i++){
        Particle& l1 = genMuons_->at(i);
        for(size_t j = i+1; j < genMuons_->size(); j++){
	  Particle& l2 = genMuons_->at(j);
	  if( l1.id() + l2.id() == 0 ){
	    Boson<Particle> Zcand(l1,l2);
	    if(ZBosonDefinition(Zcand)){
	      genZlepCandidates_->push_back(Zcand);
            }
	  }
        }
      }
      // Gen ZW --> mu+mu-mu nu
      if(genZlepCandidates_->size() >= 1 && genWlepCandidates_->size() >= 1){	
        std::sort(genZlepCandidates_->begin(), genZlepCandidates_->end(), MassComparator(phys::ZMASS));
        Boson<Particle>& genZ0 = genZlepCandidates_->front();

        std::sort(genWlepCandidates_->begin(), genWlepCandidates_->end(), MassComparator(phys::WMASS));
        auto genW = std::find_if(genWlepCandidates_->begin(), genWlepCandidates_->end(),
			     [&genZ0](auto W){ return !haveCommonDaughter(genZ0, W); }
			     );
        if(genW != genWlepCandidates_->end()){
          genZW_ = DiBoson<Particle, Particle>(genZ0, *genW);
        }
      }
      
      // -- Gen Leptons definition -- //
      phys::Particle genMu1 = genZW_.first().daughter(0);
      phys::Particle genMu2 = genZW_.first().daughter(1);
      phys::Particle genMu3 = genZW_.second().daughter(0);
      phys::Particle genNu = genZW_.second().daughter(1);
      
      
      // -- Invariant mass of the mu+mu- pair -- //
      double genMuPairInvMass = (genMu1.p4() + genMu2.p4()).M();
      theHistograms->fill("GEN_mu_pair_invariant_mass_mode4", "Gen Muon pair mu+mu- invariant mass (mode 4)", 75, 0, 150., genMuPairInvMass, theWeight);
      
      // -- Invariant mass of the 4 leptons mu+ mu- mu nu -- // 
      double genFourLepInvMass = (genMu1.p4() + genMu2.p4() + genMu3.p4() + genNu.p4()).M();
      theHistograms->fill("GEN_four_leptons_invariant_mass_mode4", "Gen Leptons mu+mu-mu nu invariant mass (mode 4)", 150, 0, 300., genFourLepInvMass, theWeight);      
      
      // -- Transverse mass of the 4 leptons mu+ mu- mu nu -- //
      double genFourLepTransverseMass = (genMu1.p4() + genMu2.p4() + genMu3.p4() + genNu.p4()).Mt();;
      theHistograms->fill("GEN_four_leptons_transverse_mass_mode4", "Gen Four Leptons mu+mu-mu nu transverse mass (mode 4)", 150, 0., 300., genFourLepTransverseMass, theWeight);
      
      // -- Invariant mass of the muon pair mu+mu- vs Invariant mass of the 4 leptons mu+mu-mu nu -- //
      theHistograms->fill("GEN_mu_pair_four_lep_inv_mass_mode4", "Gen Muon pair vs Four Leptons invariant mass (mode 4);4l mass;mu+mu- mass", 100, 0., 500., 60, 0., 300., genFourLepInvMass, genMuPairInvMass, theWeight);
      
      // -- Gen Neutrino pt vs Invariant mass of the 4 leptons mu+mu-mu nu -- //
      theHistograms->fill("GEN_nu_pt_four_lep_inv_mass_mode4", "Gen Neutrino pt vs Four Leptons invariant mass (mode 4);4l mass;nu pt", 100, 0., 500., 60, 0., 300., genFourLepInvMass, genNu.pt(), theWeight);
      
      // -- Gen Muon pt vs Invariant mass of the 4 leptons mu+mu-mu nu -- //
      theHistograms->fill("GEN_mu_pt_four_lep_inv_mass_mode4", "Gen Muon pt vs Four Leptons invariant mass (mode 4);4l mass;mu pt", 100, 0., 500., 60, 0., 300., genFourLepInvMass, genMu3.pt(), theWeight);
      
      // -- Transverse mass of the Gen Lepton pair mu nu vs Invariant mass of the 4 Gen Leptons mu+mu-mu nu -- //
      double genLepPairTransverseMass = (genMu3.p4() + genNu.p4()).Mt();
      theHistograms->fill("GEN_lep_pair_transv_mass_mode4", "Gen Lepton pair mu nu transverse mass (mode 4)", 75, 0., 150., genLepPairTransverseMass, theWeight);
      theHistograms->fill("GEN_lep_pair_transv_mass_four_lep_inv_mass_mode4", "Gen Lepton pair mu nu transverse mass vs Four Leptons invariant mass (mode 4);4l mass;mu nu mass", 100, 0., 500., 75, 0., 150., genFourLepInvMass, genLepPairTransverseMass, theWeight);
      
      // -- Transverse vs Invariant mass of the 4 gen leptons mu+mu-mu nu -- // 
      theHistograms->fill("GEN_four_lep_transv_mass_four_lep_inv_mass_mode4", "Gen Four Lepton transverse vs invariant mass (mode 4); Inv mass; Transv mass", 100, 0., 500., 100, 0., 500., genFourLepInvMass, genFourLepTransverseMass, theWeight);
      
      
      
      
      
      
      
      
      
      
      
      // --------------------- REC LEVEL --------------------- //
      if(abs(ZW->first().daughter(0).id())==13 && abs(ZW->first().daughter(1).id())==13 && abs(ZW->second().daughter(0).id())==13){     // muons->size()==3 && electrons->size()==0 
        
        // deltaR(Z0,W)
        theHistograms->fill("REC_ZW_bosons_deltaR_mode4", "Rec ZW Bosons deltaR (mode 4)", 20, 0., 5., deltaR(ZW->first(),ZW->second()), theWeight);
        if(deltaR(ZW->first(),ZW->second())<0.1){
          theHistograms->fill("REC_ZW_bosons_deltaR_zoom_mode4", "Rec ZW Bosons deltaR zoom (mode 4)", 10, 0., 0.1, deltaR(ZW->first(),ZW->second()), theWeight);
        }      
        // Z0 & W bosons pt
        theHistograms->fill("REC_Z_pt_mode4", "Rec Z boson pt (mode 4);pt", 150, 0, 300., ZW->first().pt(), theWeight);
        theHistograms->fill("REC_W_pt_mode4", "Rec W boson pt (mode 4);pt", 150, 0, 300., ZW->second().pt(), theWeight);
        // REC bosons & REC leptons definition -- //
        Boson<phys::Lepton> recZ0 = ZW->first();
        Boson<phys::Lepton> recW = ZW->second();
        phys::Particle recMu1 = recZ0.daughter(0);
        phys::Particle recMu2 = recZ0.daughter(1);
        phys::Particle recMu3 = recW.daughter(0);
        phys::Particle recNu = recW.daughter(1);
        
        
        if(genMu3.id()==recMu3.id()){
        // -- Transverse mass of the 4 REC leptons mu+mu-mu nu -- //
        double recFourLepTransverseMass = (recMu1.p4() + recMu2.p4() + recMu3.p4() + recNu.p4()).Mt();
        theHistograms->fill("REC_four_lep_transv_mass_mode4", "Rec Four Leptons mu+mu-mu nu transverse mass (mode 4)", 250, 0., 500., recFourLepTransverseMass, theWeight);
        
        // -- Transverse mass of the 4 REC leptons mu+mu-mu nu vs Transverse mass of the 4 GEN leptons mu+mu-mu nu -- //
        theHistograms->fill("REC_four_lep_transv_mass_GEN_four_lep_transv_mass_mode4", "Rec Four Leptons transverse mass vs Gen Four Leptons transverse mass mu+mu-mu nu (mode 4);Gen 4l;Rec 4l", 100, 0., 500., 100, 0., 500., genFourLepTransverseMass, recFourLepTransverseMass, theWeight);
        
        // -- Transverse mass of the 4 REC leptons mu+mu-mu nu vs Invariant mass of the 4 GEN leptons mu+mu-mu nu -- //
        theHistograms->fill("REC_four_lep_transv_mass_GEN_four_lep_inv_mass_mode4", "Rec Four Leptons transverse mass vs Gen Four Leptons invariant mass mu+mu-mu nu (mode 4);Gen 4l;Rec 4l", 100, 0., 500., 100, 0., 500., genFourLepInvMass, recFourLepTransverseMass, theWeight);
        
        // -- Invariant mass of the rec muon pair mu+mu- vs Transverse mass of the 4 rec leptons mu+mu-mu nu -- //
        double recMuPairInvMass = (recMu1.p4() + recMu2.p4()).M();
        theHistograms->fill("REC_mu_pair_invariant_mass_mode4", "Rec Muon pair mu+mu- invariant mass (mode 4)", 75, 0, 150., recMuPairInvMass, theWeight);
        theHistograms->fill("REC_mu_pair_inv_mass_REC_four_lep_transv_mass_mode4", "Rec Muon pair invariant mass vs Rec Four Leptons transverse mass (mode 4);4l mass;mu+mu- mass", 100, 0., 500., 60, 0., 300., recFourLepTransverseMass, recMuPairInvMass, theWeight);
        
        // -- Transverse mass of the rec lepton pair mu nu vs Transverse mass of the 4 rec leptons mu+mu-mu nu -- //
        double recLepPairTransverseMass = (recMu3.p4() + recNu.p4()).Mt();
        theHistograms->fill("REC_lep_pair_transv_mass_mode4", "Rec Lepton pair mu nu transverse mass (mode 2)", 75, 0., 150., recLepPairTransverseMass, theWeight);
        theHistograms->fill("REC_lep_pair_transv_mass_REC_four_lep_transv_mass_mode4", "Rec Lepton pair mu nu transverse mass vs Rec Four Leptons transverse mass (mode 4);4l mass;mu nu mass", 100, 0., 500., 60, 0., 300., recFourLepTransverseMass, recLepPairTransverseMass, theWeight);
        
        
        
        // -- REC nu (met) pt vs Transverse mass of the 4 REC leptons mu+mu-mu nu -- //
        theHistograms->fill("REC_nu_pt_REC_four_lep_transv_mass_mode4", "Rec nu (met) pt vs Rec Four Leptons transverse mass mu+mu-mu nu (mode 4);4l mass;nu pt", 100, 0., 500., 60, 0., 300., recFourLepTransverseMass, recNu.pt(), theWeight);
        
        // -- REC mu pt vs Transverse mass of the 4 REC leptons mu+mu-mu nu -- //
        theHistograms->fill("REC_mu_pt_REC_four_lep_transv_mass_mode4", "Rec mu pt vs Rec Four Leptons transverse mass mu+mu-mu nu (mode 4);4l mass;mu pt", 100, 0., 500., 60, 0., 300., recFourLepTransverseMass, recMu3.pt(), theWeight);
        
        
        
        
        
        
        
        }
      }
      

      
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
  genZlepCandidates_->clear();
  genWlepCandidates_->clear();

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


bool WlllnuAnalyzer::isGen_mode1(double var1, double var2, double var3, double var4, double var5){
  if(var1>170. &&                         // genFourLepInvMass
     80.<var2<100. &&                     // genElPairInvMass           (var2>80.)
     80.<var3<90. &&                      // genLepPairTransverseMass   (var3>80.)
     20.<var4<50. &&                      // genMu.pt()                 (var4>20.)
     20.<var5<50. )                       // genNu.pt()                 (var5>20.)
  {
    return true;
  }
  else{
    return false;
  }
  
}


bool WlllnuAnalyzer::isRec_mode1(double var1, double var2, double var3, double var4, double var5){
  if( var1>150. &&                         // recFourLepTransverseMass
      var2>80. &&                          // recElPairInvMass          (80.<var2<100.)
      var3>70. &&                          // recLepPairTransverseMass  (70.<var3<150.)
      var4>20. &&                          // recMu.pt()                (20.<var4<50. )
      var5>20. )                           // recNu.pt()                (20.<var5<70. )
  {
    return true;
  }
  else{
    return false;
  }
  
}


bool WlllnuAnalyzer::isGen_mode2(double var1, double var2, double var3, double var4, double var5){
  if(var1>170. &&                          // genFourLepInvMass 
     80.<var2<100. &&                      // genMuPairInvMass          (var2>80.)
     80.<var3<90. &&                       // genLepPairTransverseMass  (var3>80.)
     20.<var4<50. &&                       // genEl.pt()                (var4>20.)
     20.<var5<50. )                        // genNu.pt()                (var5>20.)
  {
    return true;
  }
  else{
    return false;
  }
  
}


bool WlllnuAnalyzer::isRec_mode2(double var1, double var2, double var3, double var4, double var5){
  if( var1>150. &&                          // recFourLepTransverseMass
      var2>80. &&                           // recMuPairInvMass           (80.<var2<100.)
      var3>70. &&                           // recLepPairTransverseMass   (70.<var3<130.)
      var4>20. &&                           // recEl.pt()                 (20.<var4<50. )
      var5>20. )                            // recNu.pt()                 (20.<var5<70. )
  {
    return true;
  }
  else{
    return false;
  }
  
}



/*
void histogramFill(std::vector<phys::Particle>* genLep, string histName, string histTitle, int nBin, double start, double end, string var, double weight){  
  
  foreach(const phys::Particle &genL, *genLep){
    theHistograms->fill(histName, histTitle, nBin, start, end, genL.var, weight);
  }
  
  
}  
*/  
  
