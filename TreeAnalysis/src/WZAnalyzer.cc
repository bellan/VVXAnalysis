#include "VVXAnalysis/TreeAnalysis/interface/WZAnalyzer.h"
#include "VVXAnalysis/Commons/interface/AriEle.h"
#include "VVXAnalysis/Commons/interface/Comparators.h"
#include "VVXAnalysis/Commons/interface/Constants.h"
#include "VVXAnalysis/Commons/interface/SignalDefinitions.h"
#include "VVXAnalysis/Commons/interface/Utilities.h"

#include <boost/foreach.hpp>
#define foreach BOOST_FOREACH

#include <boost/assign/std/vector.hpp> 

using namespace boost::assign;
using namespace phys;
using namespace std;

void WZAnalyzer::begin() {
  zahl = 0;
  eventcounter = 0;
  begintime = ((float)clock())/CLOCKS_PER_SEC;
}

Int_t WZAnalyzer::cut() {
  return 1;
}


void WZAnalyzer::analyze(){
  zahl++;

  vector<Particle> electron;
  vector<Particle> muon;
  
  cout << "\n--------------------------------------------------------------------------"<<endl;
  cout << "Run: " << run << " event: " << event << " number: " << zahl << endl;
  
  foreach(const Particle &gen, *genParticles){
    if((abs(gen.id()) != 11 && abs(gen.id()) != 13) || (!(gen.genStatusFlags().test(phys::GenStatusBit::isPrompt)) || !(gen.genStatusFlags().test(phys::GenStatusBit::fromHardProcess)))) continue;
    cout << "id: " << gen.id() << " pt: " << gen.pt() << "\t eta: " << gen.eta() << endl;
    theHistograms.fill("ptAllGenParticle",  "p_{t}", 100, 0, 100, gen.pt());
    theHistograms.fill("massAllGenParticle","mass",  100, 0, 100, gen.mass());
    theHistograms.fill("etaAllGenParticle", "#eta",  100, 0, 100, gen.eta());
    theHistograms.fill("YAllGenParticle",   "Y",     100, 0, 100, gen.eta());
       
    if(gen.id() == 11)       electron.insert(electron.begin(),gen);
    else if(gen.id() == -11) electron.push_back(gen);
    else if(gen.id() == 13)  muon.insert(muon.begin(),gen);
    else if(gen.id() == -13) muon.push_back(gen);
  }

 
  // ~~~~~~ tests on ZZ Analysis ~~~~~~
  // /*
  if(electron.size()+muon.size()!=4) {
    cout << "There are not enough or too many final leptons in this event." << endl;
    return;
  }

  // Reconstruction of the two Zs
  
  vector<Ztype> Zet;
  vector<Ztype> possibleZ;

  if(electron.size()==2 && muon.size()==2){
    possibleZ.push_back(Ztype(electron[0], electron[1], 23));
    possibleZ.push_back(Ztype(muon[0], muon[1], 23));
  }
  else if(electron.size()==4){
    possibleZ.push_back(Ztype(electron[0], electron[2], 23));
    possibleZ.push_back(Ztype(electron[0], electron[3], 23));
    possibleZ.push_back(Ztype(electron[1], electron[2], 23));
    possibleZ.push_back(Ztype(electron[1], electron[3], 23));
  }
  else if(muon.size()==4){      
    possibleZ.push_back(Ztype(muon[0], muon[2], 23));
    possibleZ.push_back(Ztype(muon[0], muon[3], 23));
    possibleZ.push_back(Ztype(muon[1], muon[2], 23));
    possibleZ.push_back(Ztype(muon[1], muon[3], 23));
  }
      
  stable_sort(possibleZ.begin(), possibleZ.end(), MassComparator(ZMASS));
  Zet.push_back(possibleZ[0]);
  
  foreach(const Boson<Particle> candidates, possibleZ){
    if(candidates.daughter(0).pt() != Zet[0].daughter(0).pt() && candidates.daughter(1).pt() != Zet[0].daughter(1).pt()) Zet.push_back(candidates);
  } 

  theHistograms.fill("massZ", "Z mass", 1000, 50, 150, Zet[0].mass());
  theHistograms.fill("massZ", "Z mass", 1000, 50, 150, Zet[1].mass());
  
  cout << "Reconstructed Zs: Z0 " << Zet[0].mass() << " daughter's ID: " << Zet[0].daughter(0).id() << " and " << Zet[0].daughter(1).id() << "\n\t\t  Z1 " << Zet[1].mass() << " daughter's ID: " << Zet[1].daughter(0).id() << " and " << Zet[1].daughter(1).id() << endl;


  //Reconstruction of the ZZ
  
  ZZtype doppelZ(Zet[0], Zet[1]);

  theHistograms.fill("ptZZ",   "ZZ p_{t}", 100,  0, 100, doppelZ.pt());
  theHistograms.fill("etaZZ",  "ZZ #eta",  100,-10,  10, doppelZ.eta());
  theHistograms.fill("massZZ", "ZZ mass", 1000,  0, 400, doppelZ.mass());
  theHistograms.fill("YZZ",    "ZZ Y",     100,  0, 100, doppelZ.rapidity());
  
  cout << "Reconstructed ZZ: " << doppelZ << endl;

  //add transvers mass
  
  // */
  // ~~~~~~ end of tests on ZZ Analysis ~~~~~~
 

  // ~~~~~~ WZ Analysis ~~~~~~
  //
  /*
  if(electron.size()+muon.size()!=3)  {
    cout << "There are not enough or too many final leptons in this event." << endl;
    return;
  }
  
  // ------ ZL reconstructed ------
  
  Ztype Zet;
  vector<Particle> lepton;
  vector<Ztype> possibleZ;
  vector<Zltype> Zls;

  // first filter on pt and eta
  foreach(const Particle ele, electron){
    if(ele.pt() < 7 || abs(ele.eta()) > 2.5){
      cout << "Electrons: pt less than 7 GeV or eta's absolute value greater than 2.5" << endl;
      return;
    }
    lepton.push_back(ele);
  }
  foreach(const Particle mu, muon){
    if(mu.pt() < 5 || abs(mu.eta()) > 2.4){
      cout << "Muons: pt less than 5 GeV or eta's absolute value greater than 2.4" << endl;
      return;
    }
    lepton.push_back(mu);
  }
  stable_sort(electron.begin(), electron.end(), PtComparator());
  stable_sort(muon.begin(), muon.end(), PtComparator());
  stable_sort(lepton.begin(), lepton.end(), PtComparator());

  // second filter on pt and eta
  if(lepton[0].pt() < 20){
    cout << "First lepton pt less than 20 GeV" << endl;
    return;
  }
  if(abs(lepton[1].id() == 11) && lepton[1].pt() < 12){
    cout << "Second lepton is an electron and has pt less than 12 GeV" << endl;
    return;
  }
  if(abs(lepton[1].id() == 13) && lepton[1].pt() < 10){
    cout << "Second lepton is a muon and has pt less than 10 GeV" << endl;
    return;
  }
  
  // Zl recontructed
  if(electron.size()==3){
    possibleZ.push_back(Ztype(electron[0], electron[2], 23));
    Zls.push_back(Zltype(possibleZ[0], electron[1]));
    
    if(electron[0].id()==electron[1].id()){
      possibleZ.push_back(Ztype(electron[1], electron[2], 23));
      Zls.push_back(Zltype(possibleZ[1], electron[0]));
    }
    else{
      possibleZ.push_back(Ztype(electron[0], electron[1], 23));
      Zls.push_back(Zltype(possibleZ[1], electron[2]));
    }
  }
  
  else if(electron.size()==2 && muon.size()==1){
    Zet = Ztype(electron[0], electron[1], 23);
    Zls.push_back(Zltype(Zet, muon[0]));
  }
  
  else if(electron.size()==1 && muon.size()==2){
    Zet = Ztype(muon[0], muon[1], 23);
    Zls.push_back(Zltype(Zet, electron[0]));    
  }
  
  else if(muon.size()==3){
    possibleZ.push_back(Ztype(muon[0], muon[2], 23));
    Zls.push_back(Zltype(possibleZ[0], muon[1]));
    
    if(muon[0].id()==muon[1].id()){
      possibleZ.push_back(Ztype(muon[1], muon[2], 23));
      Zls.push_back(Zltype(possibleZ[1], muon[0]));
    }
    else{
      possibleZ.push_back(Ztype(muon[0], muon[1], 23));
      Zls.push_back(Zltype(possibleZ[1], muon[2]));
    }
  
  }

  cout << "\nZl candidates are: " << Zls.size() << endl;
  
  //
 */
   
}
  
void WZAnalyzer::end(TFile &){

  cout << "\n--------------------------------------------------------------------------"<<endl;
  
  cout << "\nNumber of event analyzed: " << eventcounter << endl;
  
  // execution time
  endtime = ((float)clock())/CLOCKS_PER_SEC;
  cout << "Execution time of the analysis: " << endtime - begintime << " seconds." << endl;
  //cout << "prova" << (int)((endtime - begintime)/3600) << " h " << (int)(((endtime - begintime) - (endtime - begintime)/3600)) << " s " << endl; to be rewritten with %
  cout << "\n--------------------------------------------------------------------------"<<endl;
}
