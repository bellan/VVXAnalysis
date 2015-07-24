#include "VVXAnalysis/TreeAnalysis/interface/VVXAnalyzer.h"
#include "VVXAnalysis/Commons/interface/SignalDefinitions.h"
#include "VVXAnalysis/Commons/interface/Utilities.h"

#include <boost/foreach.hpp>
#define foreach BOOST_FOREACH

#include <boost/assign/std/vector.hpp> 
using namespace boost::assign;

using std::cout;
using std::endl;


using namespace phys;

Int_t VVXAnalyzer::cut() {
  return 1;
}


void VVXAnalyzer::ZZplots(int id){

  if(ZZ->id() != id && id != -1) return; // -1 here means generic 4l final state

  std::string decay  = "4l";
  std::string decay1 = "l";
  std::string decay2 = "l";
  if      (id == 52) {decay = "4m"  ; decay1 = "m"; decay1 = "m";}
  else if (id == 48) {decay = "2e2m"; decay1 = "e"; decay2 = "m";}
  else if (id == 44) {decay = "4e"  ; decay1 = "e"; decay2 = "e";}
  theHistograms.fill(std::string("ZZTo")+decay+std::string("_mZ1To2")+decay1, std::string("Invariant mass of Z_{1}#rightarrow 2")+decay1,  15, 60,  120, ZZ->first().mass() , theWeight); 
  theHistograms.fill(std::string("ZZTo")+decay+std::string("_mZ2To2")+decay2, std::string("Invariant mass of Z_{2}#rightarrow 2")+decay2,  15, 60,  120, ZZ->second().mass(), theWeight); 
  theHistograms.fill(std::string("ZZTo")+decay+std::string("_mZZTo") +decay , std::string("Invariant mass of ZZ#rightarrow ")    +decay ,  40,  0, 1000, ZZ->mass()         , theWeight); 

  theHistograms.fill(std::string("ZZTo")+decay+"_nJets"       , "Number of jets (|#eta|<4.7 and p_T > 30 GeV)"        , 10, 0, 10, jets->size(), theWeight); 
  theHistograms.fill(std::string("ZZTo")+decay+"_nCentralJets", "Number of central jets (|#eta|<2.5 and p_T > 30 GeV)", 10, 0, 10, centralJets->size(), theWeight); 

  if(jets->size() >= 2)
    theHistograms.fill(std::string("ZZTo")+decay+"_deltaEtaJJ", "#Delta #eta(j,j) between the two most energetic jets",  10, 0, 8, fabs(jets->at(0).eta() - jets->at(1).eta()), theWeight); 
  
  
  if(centralJets->size() >= 2){
    theHistograms.fill(std::string("ZZTo")+decay+"_deltaEtaJJcentral", "#Delta #eta(j,j) between the two most energetyc central jets",  10, 0, 8, fabs(centralJets->at(0).eta() - centralJets->at(1).eta()), theWeight); 
    theHistograms.fill(std::string("ZZTo")+decay+"_mJJ", "m_{jj}",  20, 0, 1000, (centralJets->at(0).p4() + centralJets->at(1).p4()).M(), theWeight); 
  }

  theHistograms.fill(std::string("ZZTo")+decay+"_nExtraMuons"    , "Number of extra muons in the event"    , 10, 0, 10, muons->size(), theWeight); 
  theHistograms.fill(std::string("ZZTo")+decay+"_nExtraElectrons", "Number of extra electrons in the event", 10, 0, 10, electrons->size(), theWeight); 
  theHistograms.fill(std::string("ZZTo")+decay+"_nExtraLeptons"  , "Number of extra leptons in the event"  , 10, 0, 10, muons->size()+electrons->size(), theWeight); 

}


int VVXAnalyzer::pairing(const phys::Particle &p, const phys::Boson<phys::Particle> &vb1, const phys::Boson<phys::Particle> &vb2){
  
  double dR1 = 99999.;
  if      (p.id() == vb1.daughter(0).id()) dR1 = physmath::deltaR(p,vb1.daughter(0));
  else if (p.id() == vb1.daughter(1).id()) dR1 = physmath::deltaR(p,vb1.daughter(1));
  
  double dR2 = 99999.;
  if      (p.id() == vb2.daughter(0).id()) dR2 = physmath::deltaR(p,vb2.daughter(0));
  else if (p.id() == vb2.daughter(1).id()) dR2 = physmath::deltaR(p,vb2.daughter(1));

  if(dR1 > 0.5 && dR2 > 0.5) return 0;
  if(dR1 < dR2)              return 1;
  else                       return 2;

  
}


void VVXAnalyzer::analyze(){

  //if(event != 481932 
     // && event != 1742299
     // && event != 2655836
     // && event != 3851846
     // && event != 3852047
     // && event != 3850333
     // && event != 3853736
     // && event != 4219998
  //   ) return; 

  //cout << "------------------------------------------------------------------"<<endl;
  //cout << "Run: " << run << " event: " << event << endl;
  
  if(genVBParticles->size()<2) return;
  phys::Boson<phys::Particle> genZ1;
  phys::Boson<phys::Particle> genZ2;

  foreach(const phys::BosonParticle &gen, *genVBParticles)
    if(gen.id() == 23 && (abs(gen.daughter(0).id()) == 11 || abs(gen.daughter(0).id()) == 13)){
      if(!genZ1.isValid()) genZ1 = gen;
      else if(!genZ2.isValid()) genZ2 = gen;
      else cout << "3 Z boson in the event" << endl;
    }
  
  phys::DiBoson<phys::Particle,phys::Particle> genZZ(genZ1,genZ2);

 

  theHistograms.fill("Z1MassResolution","Resolution on Z1 mass", 100, -50, 50, genZ1.mass()-ZZ->first().mass());
  theHistograms.fill("Z2MassResolution","Resolution on Z2 mass", 100, -50, 50, genZ2.mass()-ZZ->second().mass());
  theHistograms.fill("ZZMassResolution","Resolution on ZZ mass", 100, -50, 50, genZZ.mass()-ZZ->mass());

  if(topology.test(0)){
      theHistograms.fill("Z1MassResolution_Signal","Resolution on Z1 mass", 100, -50, 50, genZ1.mass()-ZZ->first().mass());
      theHistograms.fill("Z2MassResolution_Signal","Resolution on Z2 mass", 100, -50, 50, genZ2.mass()-ZZ->second().mass());
      theHistograms.fill("ZZMassResolution_Signal","Resolution on ZZ mass", 100, -50, 50, genZZ.mass()-ZZ->mass());
    }
  if(!topology.test(0)){
      theHistograms.fill("Z1MassResolution_NoSignal","Resolution on Z1 mass", 100, -50, 50, genZ1.mass()-ZZ->first().mass());
      theHistograms.fill("Z2MassResolution_NoSignal","Resolution on Z2 mass", 100, -50, 50, genZ2.mass()-ZZ->second().mass());
      theHistograms.fill("ZZMassResolution_NoSignal","Resolution on ZZ mass", 100, -50, 50, genZZ.mass()-ZZ->mass());
    }

  int pairingCode_Z1_d0 = pairing(ZZ->first().daughter(0),genZ1,genZ2);
  int pairingCode_Z1_d1 = pairing(ZZ->first().daughter(1),genZ1,genZ2);
 
  int pairingCode_Z2_d0 = pairing(ZZ->second().daughter(0),genZ1,genZ2);
  int pairingCode_Z2_d1 = pairing(ZZ->second().daughter(1),genZ1,genZ2);
  
  if(pairingCode_Z1_d0 == pairingCode_Z1_d1 && pairingCode_Z2_d0 == pairingCode_Z2_d1 && pairingCode_Z1_d0*pairingCode_Z2_d0 != 0) {
    //cout << "Good! Right pairing combination! " << pairingCode_Z1_d0 << " " << pairingCode_Z1_d1 << " " << pairingCode_Z2_d0 << " " << pairingCode_Z2_d1 << endl;
    theHistograms.fill("Pairing","Pairing", 3, 0, 3, 2);
  }
  else{
    //cout << "Naaa! Wrong pairing combination! " << pairingCode_Z1_d0 << " " << pairingCode_Z1_d1 << " " << pairingCode_Z2_d0 << " " << pairingCode_Z2_d1 << endl;
    if(pairingCode_Z1_d0*pairingCode_Z2_d0 != 0)
      theHistograms.fill("Pairing","Pairing", 3, 0, 3, 1);
    else
      theHistograms.fill("Pairing","Pairing", 3, 0, 3, 0);
  }


  // cout << "Gen ZZ: "  << genZZ << endl;
  // cout << "Gen Bosons:" << endl
  //      << "Z1: " << genZ1 << endl
  //      << genZ1.daughter(0) << endl
  //      << genZ1.daughter(1) << endl
  //      << "Z2: " << genZ2 << endl
  //      << genZ2.daughter(0) << endl
  //      << genZ2.daughter(1) << endl;


  // cout << "Reco ZZ: " << *ZZ << endl;
  // cout << "Z1: " << ZZ->first()   << endl
  //      << ZZ->first().daughter(0) << endl
  //      << ZZ->first().daughter(1) << endl
  //      << "Z2: " << ZZ->second()  << endl
  //      << ZZ->second().daughter(0) << endl
  //      << ZZ->second().daughter(1) << endl;

  
  // cout << "Gen masses, ZZ = " << genZZ.mass() << " Z1 = " << genZ1.mass() << " Z2= " << genZ2.mass() << endl;
  // cout << "Reco masses, ZZ = " << ZZ->mass() << " Z1 = " << ZZ->first().mass() << " Z2= " << ZZ->second().mass() << endl;
  // cout << endl;

  // Some basic plots on ZZ
  //ZZplots();   // ZZ --> 4l
  //ZZplots(52); // ZZ --> 4m
  //ZZplots(48); // ZZ --> 2e2m
  //ZZplots(44); // ZZ --> 4e




}
  
