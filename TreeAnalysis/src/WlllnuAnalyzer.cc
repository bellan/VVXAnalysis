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

using namespace std;

using namespace phys;
using namespace physmath;


void WlllnuAnalyzer::begin(){
  genElectrons_ = new std::vector<Particle>;
  genMuons_ = new std::vector<Particle>;
  leptons = new std::vector<Lepton>;
  
  
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
  
  
  
  // -- W boson to 3 leptons decay four possible configurations -- //
  int event_mode = eventMode(genElectrons_, genMuons_, electrons, muons);
  bool gen_cut_m3lnu = false;
  
  bool rec_cut_3l = false;
  bool rec_cut_mT = false;
  bool rec_cut_llmI = false;
  bool rec_cut_lpt = false;
  bool rec_cut_pt3l = false;

  if(event_mode == 1){         // MODE 1 e+ e- mu nu
    
    // ------------ GENERATED LEPTONS e+e-mu nu ------------ //
    if( genElectrons_->size()==2 && genMuons_->size()==1 
        && (genElectrons_->at(0).id() + genElectrons_->at(1).id())==0
        && genElectrons_->at(0).pt()>10 && genElectrons_->at(1).pt()>10 && genMuons_->at(0).pt()>10
        && -2.5<genElectrons_->at(0).eta() && genElectrons_->at(0).eta()<2.5
        && -2.5<genElectrons_->at(1).eta() && genElectrons_->at(1).eta()<2.5 
        && -2.5<genMuons_->at(0).eta() && genMuons_->at(0).eta()<2.5 ){  
      //cout << "Provaaaaaaa1" << endl;
      genEvents = true;  
      const phys::Particle& genEl1 = genElectrons_->at(0);
      const phys::Particle& genEl2 = genElectrons_->at(1);
      const phys::Particle& genMu = genMuons_->at(0);
      const phys::Particle& genNu = genNeutrinos_->at(0);

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

      // -- GEN LEVEL SIGNAL DEFINITION -- //
      gen_cut_m3lnu = isGen_mode1(genFourLepInvMass);
    }
    else{
      genEvents = false;
    }

    
    // ------------------ RECONSTRUCTED LEPTONS e+ e- mu nu -------------------------- //
    if(electrons->size()==2 && muons->size()==1
       && (electrons->at(0).id() + electrons->at(1).id())==0 ){
      phys::Particle recEl1 = electrons->at(0);
      phys::Particle recEl2 = electrons->at(1);
      phys::Particle recMu = muons->at(0);
      phys::Particle recNu = *met;
      //cout << "Prova1" << endl; 
              
      // -- 3l e+e-mu events recostruction efficiency -- //
      theHistograms->fill("GEN_REC_3l_events_mode1", "3l events recostruction efficiency (mode 1)", 3, 0., 3., 1, theWeight);


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
      
      // -- REC mu pt -- //
      theHistograms->fill("REC_mu_pt_mode1", "Rec Mu pt (mode 1)", 150, 0., 300., recMu.pt(), theWeight);
      theHistograms->fill("REC_mu_pt_REC_four_lep_transv_mass_mode1", "Rec mu pt vs Rec Four Leptons transverse mass e+e-mu nu (mode 1); 4l mass;mu pt", 100, 0., 500., 30, 0., 150., recFourLepTransverseMass, recMu.pt(), theWeight);
      
      // -- REC e+e- pair Invariant mass vs REC mu pt -- //
      theHistograms->fill("REC_el_pair_inv_mass_REC_mu_pt_mode1", "Rec e+e- pair invariant vs Rec mu pt (mode 1); mu pt;e+e- inv mass", 100, 0., 500., 30, 0., 150., recMu.pt(), recElPairInvMass, theWeight);
      
      // -- REC charged leptons e+e-mu pt sum -- //
      theHistograms->fill("REC_ch_lep_pt_sum_mode1", "Rec Charged Lepton e+e-mu pt sum (mode 1)", 150, 0., 300., recEl1.pt()+recEl2.pt()+recMu.pt(), theWeight);
      theHistograms->fill("REC_ch_lep_pt_sum_REC_four_lep_transv_mass_mode1", "Rec Charged Lepton e+e-mu pt sum vs Rec Four Leptons transverse mass e+e-mu nu (mode 1);4l transv mass;e+e-mu pt", 100, 0., 500., 60, 0., 300., recFourLepTransverseMass, recEl1.pt()+recEl2.pt()+recMu.pt(), theWeight);
      theHistograms->fill("REC_el_pair_inv_mass_REC_ch_lep_pt_sum_mode1", "Rec e+e- pair invariant vs Rec Charged Lepton e+e-mu pt sum (mode 1);e+e-mu pt;e+e- inv mass", 100, 0., 500., 60, 0., 300., recEl1.pt()+recEl2.pt()+recMu.pt(), recElPairInvMass, theWeight);
      theHistograms->fill("REC_mu_pt_REC_ch_lep_pt_sum_mode1", "Rec mu pt vs Rec Charged Lepton e+e-mu pt sum (mode 1);e+e-mu pt;mu pt", 100, 0., 500., 60, 0., 300., recEl1.pt()+recEl2.pt()+recMu.pt(), recMu.pt(), theWeight);  
      
      // ------------------------------ REC LEVEL SIGNAL REGION DEFINITION ----------------------------- //
      //rec_cut_3l = true;
      rec_cut_mT = isRec_mode1(recFourLepTransverseMass);
      rec_cut_llmI = isRec_mode11(recElPairInvMass);
      rec_cut_lpt = isRec_mode111(recMu.pt());
      rec_cut_pt3l = isRec_mode1111(recEl1.pt()+recEl2.pt()+recMu.pt());
      
    }
      
  }    
      
// ------------------------------------------------------------------------------------------------------------------------------------------------------------------------ //
  
  else if(event_mode == 2){   // MODE 2 mu+ mu- e nu
    
    // ---------- GENERATED LEPTONS mu+mu-e nu ---------- //
    if( genMuons_->size()==2 && genElectrons_->size()==1
        && (genMuons_->at(0).id() + genMuons_->at(1).id())==0
        && genMuons_->at(0).pt()>10 && genMuons_->at(1).pt()>10 && genElectrons_->at(0).pt()>10
        && -2.5<genMuons_->at(0).eta() && genMuons_->at(0).eta()<2.5
        && -2.5<genMuons_->at(1).eta() && genMuons_->at(1).eta()<2.5 
        && -2.5<genElectrons_->at(0).eta() && genElectrons_->at(0).eta()<2.5 ){
      //cout << "Provaaaaaaa2" << endl;
      genEvents = true;    
      const phys::Particle& genMu1 = genMuons_->at(0);
      const phys::Particle& genMu2 = genMuons_->at(1);
      const phys::Particle& genEl = genElectrons_->at(0);
      const phys::Particle& genNu = genNeutrinos_->at(0);
      
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
      
      // -- GEN LEVEL SIGNAL DEFINITION -- //
      gen_cut_m3lnu = isGen_mode2(genFourLepInvMass);
    }  
    else{
      genEvents = false;
    }  
      

    // ---------- RECONSTRUCTED LEPTONS mu+ mu- e nu -------------- //
    if(muons->size()==2 && electrons->size()==1 
       && (muons->at(0).id() + muons->at(1).id())==0 /*&& genElectrons_->at(0).id()==electrons->at(0).id()*/
       /*&& 60.<(muons->at(0).p4()+muons->at(1).p4()+electrons->at(0).p4()+met->p4()).Mt() && (muons->at(0).p4()+muons->at(1).p4()+electrons->at(0).p4()+met->p4()).Mt()<120.*/
       /*&& 4.<(muons->at(0).p4() + muons->at(1).p4()).M() && (muons->at(0).p4() + muons->at(1).p4()).M()<12.*/
       /*&& jets->size()==0*/){
      phys::Particle recMu1 = muons->at(0);
      phys::Particle recMu2 = muons->at(1);
      phys::Particle recEl = electrons->at(0);
      phys::Particle recNu = *met;
      //cout << "Prova2" << endl;
        
      //// -- 3l mu+mu-e events recostruction efficiency -- //
      //theHistograms->fill("GEN_REC_3l_events_mode2", "3l events recostruction efficiency (mode 2)", 3, 0., 3., 1, theWeight);
                
                
      // -- Transverse mass of the 4 REC leptons mu+mu-e nu -- //
      double recFourLepTransverseMass = (recMu1.p4() + recMu2.p4() + recEl.p4() + recNu.p4()).Mt();
      theHistograms->fill("REC_four_lep_transv_mass_mode2", "Rec Four Leptons mu+mu-e nu transverse mass (mode 2)", 250, 0., 500., recFourLepTransverseMass, theWeight);
        
      // -- Transverse mass of the 4 REC leptons mu+mu-e nu vs Transverse mass of the 4 GEN leptons mu+mu-e nu -- //
      theHistograms->fill("REC_four_lep_transv_mass_GEN_four_lep_transv_mass_mode2", "Rec Four Leptons transverse mass vs Gen Four Leptons transverse mass mu+mu-e nu (mode 2);Gen 3lnu;Rec 3lnu", 100, 0., 500., 100, 0., 500., genFourLepTransverseMass, recFourLepTransverseMass, theWeight);
        
      // -- Transverse mass of the 4 REC leptons e+e-mu nu vs Invariant mass of the 4 GEN leptons mu+mu-e nu -- //
      theHistograms->fill("REC_four_lep_transv_mass_GEN_four_lep_inv_mass_mode2", "Rec Four Leptons transverse mass vs Gen Four Leptons invariant mass mu+mu-e nu (mode 2);Gen 3lnu;Rec 3lnu", 100, 0., 500., 100, 0., 500., genFourLepInvMass, recFourLepTransverseMass, theWeight);
        
      // -- Invariant mass of the rec muon pair mu+mu- vs Transverse mass of the 4 rec leptons mu+mu-e nu -- //
      double recMuPairInvMass = (recMu1.p4() + recMu2.p4()).M();
      theHistograms->fill("REC_mu_pair_invariant_mass_mode2", "Rec Muon pair mu+mu- invariant mass (mode 2)", 75, 0, 150., recMuPairInvMass, theWeight);
      theHistograms->fill("REC_mu_pair_inv_mass_REC_four_lep_transv_mass_mode2", "Rec Muon pair invariant mass vs Rec Four Leptons transverse mass (mode 2);3lnu mass;mu+mu- mass", 100, 0., 500., 60, 0., 300., recFourLepTransverseMass, recMuPairInvMass, theWeight);
        
      // -- Transverse mass of the rec lepton pair e nu vs Transverse mass of the 4 rec leptons mu+mu-e nu -- //
      double recLepPairTransverseMass = (recEl.p4() + recNu.p4()).Mt();
      theHistograms->fill("REC_lep_pair_transv_mass_mode2", "Rec Lepton pair e nu transverse mass (mode 2)", 75, 0., 150., recLepPairTransverseMass, theWeight);
      theHistograms->fill("REC_lep_pair_transv_mass_REC_four_lep_transv_mass_mode2", "Rec Lepton pair e nu transverse mass vs Rec Four Leptons transverse mass (mode 2);4l mass;e nu mass", 100, 0., 500., 60, 0., 300., recFourLepTransverseMass, recLepPairTransverseMass, theWeight);
        
      // -- REC nu (met) pt vs Transverse mass of the 4 REC leptons mu+mu-e nu -- //
      theHistograms->fill("REC_nu_pt_REC_four_lep_transv_mass_mode2", "Rec nu (met) pt vs Rec Four Leptons transverse mass mu+mu-e nu (mode 2);3lnu mass;nu pt", 100, 0., 500., 60, 0., 300., recFourLepTransverseMass, recNu.pt(), theWeight);
        
      // -- REC el pt -- //
      theHistograms->fill("REC_el_pt_mode2", "Rec El pt (mode 2)", 150, 0., 300., recEl.pt(), theWeight);
      theHistograms->fill("REC_el_pt_REC_four_lep_transv_mass_mode2", "Rec el pt vs Rec Four Leptons transverse mass mu+mu-e nu (mode 2);3lnu mass;e pt", 100, 0., 500., 60, 0., 300., recFourLepTransverseMass, recEl.pt(), theWeight);
        
      // -- REC mu+mu- pair Invariant mass vs REC e pt -- //
      theHistograms->fill("REC_mu_pair_inv_mass_REC_el_pt_mode2", "Rec mu+mu- pair invariant vs Rec e pt (mode 2); e pt;mu+mu- inv mass", 100, 0., 500., 30, 0., 150., recEl.pt(), recMuPairInvMass, theWeight);
        
      // -- REC charged leptons e+e-mu pt sum -- //
      theHistograms->fill("REC_ch_lep_pt_sum_mode2", "Rec Charged Lepton mu+mu-e pt sum (mode 2)", 150, 0., 300., recMu1.pt()+recMu2.pt()+recEl.pt(), theWeight);
      theHistograms->fill("REC_ch_lep_pt_sum_REC_four_lep_transv_mass_mode2", "Rec Charged Lepton mu+mu-e pt sum vs Rec Four Leptons transverse mass mu+mu-e nu (mode 2);4l transv mass;mu+mu-e pt", 100, 0., 500., 60, 0., 300., recFourLepTransverseMass, recMu1.pt()+recMu2.pt()+recEl.pt(), theWeight);
      theHistograms->fill("REC_mu_pair_inv_mass_REC_ch_lep_pt_sum_mode2", "Rec mu+mu- pair invariant vs Rec Charged Lepton mu+mu-e pt sum (mode 2);mu+mu-e pt;mu+mu- inv mass", 100, 0., 500., 60, 0., 300., recMu1.pt()+recMu2.pt()+recEl.pt(), recMuPairInvMass, theWeight);
      theHistograms->fill("REC_el_pt_REC_ch_lep_pt_sum_mode2", "Rec e pt vs Rec Charged Lepton mu+mu-e pt sum (mode 2);mu+mu-e pt;e pt", 100, 0., 500., 60, 0., 300., recMu1.pt()+recMu2.pt()+recEl.pt(), recEl.pt(), theWeight);
        
        
      // ------------------------------ REC LEVEL SIGNAL DEFINITION ----------------------------- //
      //rec_cut_3l = true;
      rec_cut_mT = isRec_mode2(recFourLepTransverseMass);
      rec_cut_llmI = isRec_mode22(recMuPairInvMass);
      rec_cut_lpt = isRec_mode222(recEl.pt());
      rec_cut_pt3l = isRec_mode2222(recMu1.pt()+recMu2.pt()+recEl.pt());
      
    }
        
        

  
  }    
      
// ------------------------------------------------------------------------------------------------------------------------------------------------------------------------ //      
  
  else if(eventMode(genElectrons_, genMuons_, electrons, muons)==3){    // MODE 3 e+e-e nu
  
    // ---------- GENERATED LEPTONS e+e-e nu ---------- //
    if( genElectrons_->size()==3 && genMuons_->size()==0
        && checkLeptonsCharge(genElectrons_->at(0), genElectrons_->at(1), genElectrons_->at(2))==true
        && genElectrons_->at(0).pt()>10 && genElectrons_->at(1).pt()>10 && genElectrons_->at(2).pt()>10
        && -2.5<genElectrons_->at(0).eta() && genElectrons_->at(0).eta()<2.5
        && -2.5<genElectrons_->at(1).eta() && genElectrons_->at(1).eta()<2.5 
        && -2.5<genElectrons_->at(2).eta() && genElectrons_->at(2).eta()<2.5 ){
        
      genEvents = true;
      
      
      // -- Invariant mass of the 4 gen leptons e+ e- mu nu -- // 
      double genFourLepInvMass = (genElectrons_->at(0).p4() + genElectrons_->at(1).p4() + genElectrons_->at(2).p4() + genNeutrinos_->at(0).p4()).M();
      theHistograms->fill("GEN_four_leptons_invariant_mass_mode3", "Gen Leptons e+e-e nu invariant mass (mode 3)", 150, 0., 300., genFourLepInvMass, theWeight);
      
      
      
      // -- GEN LEVEL SIGNAL DEFINITION -- //
      gen_cut_m3lnu = isGen_mode3(genFourLepInvMass);
      
      
      // -- Min & Max Invariant mass of the 2 possible e+e- pair -- //
      std::vector<Boson<Particle> > genChElPairs = possibleLepPair(genElectrons_);
      if(gen_cut_m3lnu){
        Boson<Particle> genZ0Min_ = minInvMassChLepPair(genChElPairs);
        theHistograms->fill("GEN_signal_min_el_pair_inv_mass_mode3", "GEN signal minimum electron pair inv mass", 150, 0., 300., genZ0Min_.p4().M(), theWeight);
      }
      else{
        Boson<Particle> genZ0Min_ = minInvMassChLepPair(genChElPairs);
        theHistograms->fill("GEN_not_signal_min_el_pair_inv_mass_mode3", "GEN not signal minimum electron pair inv mass", 150, 0., 300., genZ0Min_.p4().M(), theWeight);  
      }
      
      if(gen_cut_m3lnu){
        Boson<Particle> genZ0Max_ = maxInvMassChLepPair(genChElPairs);
        theHistograms->fill("GEN_signal_max_el_pair_inv_mass_mode3", "GEN signal maximum electron pair inv mass", 150, 0., 300., genZ0Max_.p4().M(), theWeight);
      }
      else{
        Boson<Particle> genZ0Max_ = maxInvMassChLepPair(genChElPairs);
        theHistograms->fill("GEN_not_signal_max_el_pair_inv_mass_mode3", "GEN not signal maximum electron pair inv mass", 150, 0., 300., genZ0Max_.p4().M(), theWeight);  
      }
      
      
      
    }
    else{
      genEvents = false;
    }
    
    
    /*
    // ------------------- RECOSTRUCTED LEPTONS e+ e- e nu ----------------- //
    if(electrons->size()==3 && muons->size()==0 
       && checkLeptonsCharge(electrons->at(0), electrons->at(1), electrons->at(2))==true){
    		
    		
    		
    
    		
    		
    } 
     */ 
  
    
    
  }
  
// ------------------------------------------------------------------------------------------------------------------------------------------------------------------------ //
  
  
  else if(eventMode(genElectrons_, genMuons_, electrons, muons)==4){      // MODE 4 mu+mu-mu nu
    
    // ---------- GENERATED LEPTONS e+e-e nu ---------- //
    if( genMuons_->size()==3 && genElectrons_->size()==0
        && checkLeptonsCharge(genMuons_->at(0), genMuons_->at(1), genMuons_->at(2))==true
        && genMuons_->at(0).pt()>10 && genMuons_->at(1).pt()>10 && genMuons_->at(2).pt()>10
        && -2.5<genMuons_->at(0).eta() && genMuons_->at(0).eta()<2.5
        && -2.5<genMuons_->at(1).eta() && genMuons_->at(1).eta()<2.5 
        && -2.5<genMuons_->at(2).eta() && genMuons_->at(2).eta()<2.5 ){
      
      genEvents = true;
      
      
      // -- Invariant mass of the 4 gen leptons e+ e- mu nu -- // 
      double genFourLepInvMass = (genMuons_->at(0).p4() + genMuons_->at(1).p4() + genMuons_->at(2).p4() + genNeutrinos_->at(0).p4()).M();
      theHistograms->fill("GEN_four_leptons_invariant_mass_mode4", "Gen Leptons mu+mu-mu nu invariant mass (mode 4)", 150, 0., 300., genFourLepInvMass, theWeight);
      
      
      
      // -- GEN LEVEL SIGNAL DEFINITION -- //
      gen_cut_m3lnu = isGen_mode4(genFourLepInvMass);
      
      
      // -- Min & Max Invariant mass of the 2 possible mu+mu- pair -- //
      std::vector<Boson<Particle> > genChMuPairs = possibleLepPair(genMuons_);
      if(gen_cut_m3lnu){
        Boson<Particle> genZ0Min_ = minInvMassChLepPair(genChMuPairs);
        theHistograms->fill("GEN_signal_min_mu_pair_inv_mass_mode4", "GEN signal minimum muon pair inv mass", 150, 0., 300., genZ0Min_.p4().M(), theWeight);
      }
      else{
        Boson<Particle> genZ0Min_ = minInvMassChLepPair(genChMuPairs);
        theHistograms->fill("GEN_not_signal_min_mu_pair_inv_mass_mode4", "GEN not signal minimum muon pair inv mass", 150, 0., 300., genZ0Min_.p4().M(), theWeight);  
      }
      
      if(gen_cut_m3lnu){
        Boson<Particle> genZ0Max_ = maxInvMassChLepPair(genChMuPairs);
        theHistograms->fill("GEN_signal_max_mu_pair_inv_mass_mode4", "GEN signal maximum muon pair inv mass", 150, 0., 300., genZ0Max_.p4().M(), theWeight);
      }
      else{
        Boson<Particle> genZ0Max_ = maxInvMassChLepPair(genChMuPairs);
        theHistograms->fill("GEN_not_signal_max_mu_pair_inv_mass_mode4", "GEN not signal maximum muon pair inv mass", 150, 0., 300., genZ0Max_.p4().M(), theWeight);  
      }
      
    }
    else{
      genEvents = false;
    }  
                   
    
    
    
    
    
    
    
    
    /*  
    // ------------------- RECOSTRUCTED LEPTONS mu+ mu- mu nu ----------------- //
    if(muons->size()==3 && electrons->size()==0 
       && checkLeptonsCharge(muons->at(0), muons->at(1), muons->at(2))==true){
    		
    }
    */
      
    
    
    
  }  
  
  
// ------------------------------------------------------------------------------------------------------------------------------------------------------------------------ //  
    
    
  if(event_mode != 0){  
    bool gen_cut = genEvents && gen_cut_m3lnu;
    bool rec_cut = rec_cut_mT && rec_cut_llmI && rec_cut_lpt && rec_cut_pt3l; //rec_cut_3l
    
    if(gen_cut){      
      theHistograms->fill(Form("GEN_signal_mode%d", event_mode), Form("Gen signal (mode %d)", event_mode), 3, 0., 3., 1, theWeight);
      if(rec_cut){
	// -- Signal Efficiency -- //       
	theHistograms->fill(Form("GEN_REC_signal_mode%d", event_mode), Form("Rec Signal Efficiency (mode %d)", event_mode), 3, 0., 3., 1, theWeight);
      }
      else{
        // -- Gen Signal but Not Rec Events -- //
        theHistograms->fill(Form("Leptons_num_GEN_not_REC_signal_mode%d", event_mode), Form("Leptons number in Gen Not Rec Signal events (mode %d)", event_mode), 21, -0.5, 20.5, leptons->size(), theWeight);
        theHistograms->fill(Form("GEN_not_REC_signal_mode%d", event_mode), Form("Gen Not Rec Signal events (mode %d)", event_mode), 3, 0., 3., 1, theWeight);
        /*for(int i;i<leptons->size();i++){
        
        }*/
      }
    
      
      // -------------------- CUTFLOW --------------------- //
      // Bin 0 //
      theHistograms->fill(Form("GEN_signal_CUTFLOW_mode%d", event_mode), Form("Gen signal CUTFLOW (mode %d)", event_mode), 10, -0.5, 9.5, 0, theWeight);
      if(leptons->size()==3){
        // Bin 1 //
        theHistograms->fill(Form("GEN_signal_CUTFLOW_mode%d", event_mode), Form("Gen signal CUTFLOW (mode %d)", event_mode), 10, -0.5, 9.5, 1, theWeight);
        if(event_mode==1){
          if(electrons->size()==2 && muons->size()==1){
            // Bin 2 (mode 1) //
            theHistograms->fill(Form("GEN_signal_CUTFLOW_mode%d", event_mode), Form("Gen signal CUTFLOW (mode %d)", event_mode), 10, -0.5, 9.5, 2, theWeight);
            if((electrons->at(0).id()+electrons->at(1).id())==0){
              // Bin 3 (mode 1)
              theHistograms->fill(Form("GEN_signal_CUTFLOW_mode%d", event_mode), Form("Gen signal CUTFLOW (mode %d)", event_mode), 10, -0.5, 9.5, 3, theWeight);
              if( isRec_mode1( ( electrons->at(0).p4()+electrons->at(1).p4()+muons->at(0).p4()+met->p4() ).Mt() ) ){
                // Bin 4 (mode 1) mT
                theHistograms->fill(Form("GEN_signal_CUTFLOW_mode%d", event_mode), Form("Gen signal CUTFLOW (mode %d)", event_mode), 10, -0.5, 9.5, 4, theWeight);
                if(isRec_mode11( ( electrons->at(0).p4()+electrons->at(1).p4() ).M() )){
                  // Bin 5 (mode 1) mll
                  theHistograms->fill(Form("GEN_signal_CUTFLOW_mode%d", event_mode), Form("Gen signal CUTFLOW (mode %d)", event_mode), 10, -0.5, 9.5, 5, theWeight);
                  if(isRec_mode111( muons->at(0).pt() )){
                    // Bin 6 (mode1) lpt
                    theHistograms->fill(Form("GEN_signal_CUTFLOW_mode%d", event_mode), Form("Gen signal CUTFLOW (mode %d)", event_mode), 10, -0.5, 9.5, 6, theWeight);
                    if(isRec_mode1111( electrons->at(0).pt()+electrons->at(1).pt()+muons->at(0).pt() )){
                      // Bin 7 (mode 1) 3lpt
                      theHistograms->fill(Form("GEN_signal_CUTFLOW_mode%d", event_mode), Form("Gen signal CUTFLOW (mode %d)", event_mode), 10, -0.5, 9.5, 7, theWeight);
                    }
                  }
                }
              }
            }
          }
        }
        else if(event_mode==2){
          if(electrons->size()==1 && muons->size()==2){
            // Bin 2 (mode 2) //
            theHistograms->fill(Form("GEN_signal_CUTFLOW_mode%d", event_mode), Form("Gen signal CUTFLOW (mode %d)", event_mode), 10, -0.5, 9.5, 2, theWeight);
            if((muons->at(0).id()+muons->at(1).id())==0){
              // Bin 3 (mode 2) //
              theHistograms->fill(Form("GEN_signal_CUTFLOW_mode%d", event_mode), Form("Gen signal CUTFLOW (mode %d)", event_mode), 10, -0.5, 9.5, 3, theWeight);
              if( isRec_mode1( ( muons->at(0).p4()+muons->at(1).p4()+electrons->at(0).p4()+met->p4() ).Mt() ) ){
                // Bin 4 (mode 2) mT
                theHistograms->fill(Form("GEN_signal_CUTFLOW_mode%d", event_mode), Form("Gen signal CUTFLOW (mode %d)", event_mode), 10, -0.5, 9.5, 4, theWeight);
                if(isRec_mode11( ( muons->at(0).p4()+muons->at(1).p4() ).M() )){
                  // Bin 5 (mode 2) mll
                  theHistograms->fill(Form("GEN_signal_CUTFLOW_mode%d", event_mode), Form("Gen signal CUTFLOW (mode %d)", event_mode), 10, -0.5, 9.5, 5, theWeight);
                  if(isRec_mode111( electrons->at(0).pt() )){
                    // Bin 6 (mode2) lpt
                    theHistograms->fill(Form("GEN_signal_CUTFLOW_mode%d", event_mode), Form("Gen signal CUTFLOW (mode %d)", event_mode), 10, -0.5, 9.5, 6, theWeight);
                    if(isRec_mode1111( muons->at(0).pt()+muons->at(1).pt()+electrons->at(0).pt() )){
                      // Bin 7 (mode 2) 3lpt
                      theHistograms->fill(Form("GEN_signal_CUTFLOW_mode%d", event_mode), Form("Gen signal CUTFLOW (mode %d)", event_mode), 10, -0.5, 9.5, 7, theWeight);
                    }
                  }
                }
              }
            }
          }
        }
      }
      // ------------------------------------------------- //
      
    }
    
    
    else{
      theHistograms->fill(Form("not_GEN_signal_mode%d", event_mode), Form("Background Efficiency (mode %d)", event_mode), 3, 0., 3., 1, theWeight);
      if(rec_cut){
	// -- Background Efficiency -- //            
	theHistograms->fill(Form("not_GEN_REC_signal_mode%d", event_mode), Form("Background Efficiency (mode %d)", event_mode), 3, 0., 3., 1, theWeight);
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
  leptons->clear();

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
  
  // -- Sort leptons -- //
  for(auto p : *electrons){
    leptons->push_back(p);
  }
  for(auto p : *muons){
    leptons->push_back(p);
  }
  
  std::sort(leptons ->begin(), leptons  ->end(), [](const Particle& a, const Particle& b){ return a.pt() < b.pt(); });
  
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


int WlllnuAnalyzer::eventMode(std::vector<phys::Particle>* genEl_, std::vector<phys::Particle>* genMu_, std::vector<phys::Lepton>* el, std::vector<phys::Lepton>* mu){
  if( (genEl_->size()==2 && genMu_->size()==1) || (el->size()==2 && mu->size()==1) ){
    return 1;
  }
  else if( (genEl_->size()==1 && genMu_->size()==2) || (el->size()==1 && mu->size()==2) ){
    return 2;
  }
  else if( (genEl_->size()==3 && genMu_->size()==0) || (el->size()==3 && mu->size()==0) ){ 
    return 3;
  }
  else if( (genEl_->size()==0 && genMu_->size()==3) || (el->size()==0 && mu->size()==3) ){
    return 4;
  }
  else{
    return 0;
  }
}


bool WlllnuAnalyzer::isGen_mode1(double var1){
  return 70.<var1 && var1<90;
}

bool WlllnuAnalyzer::isRec_mode1(double var1){   // double var2, double var3, double var4
  return 60.<var1 && var1<120.;                  // && 4.<var2 && var2<12. && 10.<var3 && var3<30. && 30<var4 && var4<70
}
bool WlllnuAnalyzer::isRec_mode11(double var2){
  return 4.<var2 && var2<12.;  
}
bool WlllnuAnalyzer::isRec_mode111(double var3){
  return 10.<var3 && var3<30.;  
}
bool WlllnuAnalyzer::isRec_mode1111(double var4){
  return 30.<var4 && var4<70.;  
}


bool WlllnuAnalyzer::isGen_mode2(double var1){
  return 70.<var1 && var1<90.;
}

bool WlllnuAnalyzer::isRec_mode2(double var1){   // double var2, double var3, double var4
  return 60.<var1 && var1<120.;                  // && 4.<var2 && var2<12. && 10.<var3 && var3<30. && 30.<var4 && var4<70.
}
bool WlllnuAnalyzer::isRec_mode22(double var2){   
  return 4.<var2 && var2<12.;
}
bool WlllnuAnalyzer::isRec_mode222(double var3){
  return 10.<var3 && var3<30.;  
}
bool WlllnuAnalyzer::isRec_mode2222(double var4){
  return 30.<var4 && var4<70.;  
}


bool WlllnuAnalyzer::isGen_mode3(double var1){
  return 70.<var1 && var1<90;
}



bool WlllnuAnalyzer::isGen_mode4(double var1){
  return 70.<var1 && var1<90;
}







std::vector<Boson<Particle> > WlllnuAnalyzer::possibleLepPair(std::vector<phys::Particle>* genChLeptons){
  std::vector<Boson<Particle> > genZlepCandidates;
  for(int i=0;i<genChLeptons->size();i++){
    Particle& l1 = genChLeptons->at(i);
    for(int j=i;j<genChLeptons->size();j++){
      Particle& l2 = genChLeptons->at(j);
      if((l1.id()+l2.id())==0){
        Boson<Particle> Zcand(l1,l2);
        genZlepCandidates.push_back(Zcand);
      }
    }
  }
  std::sort(genZlepCandidates.begin(), genZlepCandidates.end(), [](const Boson<Particle>& a, const Boson<Particle>& b){ return (a.daughter(0).p4()+a.daughter(1).p4()).M() < (b.daughter(0).p4()+b.daughter(1).p4()).M(); });
  return genZlepCandidates;
}  
  
Boson<Particle> WlllnuAnalyzer::minInvMassChLepPair(std::vector<Boson<Particle> > genChLepPairs){
  Boson<Particle>& Z0 = genChLepPairs.front();
  return Z0;
}

Boson<Particle> WlllnuAnalyzer::maxInvMassChLepPair(std::vector<Boson<Particle> > genChLepPairs){
  Boson<Particle>& Z0 = genChLepPairs.back();
  return Z0;
}


bool WlllnuAnalyzer::checkLeptonsCharge(Particle lep1, Particle lep2, Particle lep3){
	if(abs(lep1.charge()+lep2.charge()+lep3.charge())==1){
		return true;
	}
	else{
		return false;
	}

}



