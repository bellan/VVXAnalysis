#include "VVXAnalysis/TreeAnalysis/interface/WZAnalyzer.h"
#include "VVXAnalysis/Commons/interface/AriEle.h"
#include "VVXAnalysis/Commons/interface/Colours.h"
#include "VVXAnalysis/Commons/interface/Comparators.h"
#include "VVXAnalysis/Commons/interface/Constants.h"
#include "VVXAnalysis/Commons/interface/SignalDefinitions.h"
#include "VVXAnalysis/Commons/interface/Utilities.h"

#include <boost/foreach.hpp>
#define foreach BOOST_FOREACH

#include <boost/assign/std/vector.hpp> 

using namespace boost::assign;
using namespace colour;
using namespace phys;
using namespace std;

void WZAnalyzer::begin() {

  nunumber = 0;
  totalevent = 0;
  WZevent = 0;
  zahl = 0;

  begintime = ((float)clock())/CLOCKS_PER_SEC;
}

Int_t WZAnalyzer::cut() {
  return 1;
}


void WZAnalyzer::analyze(){
  zahl++;

  vector<Particle> electron;
  vector<Particle> muon;
  vector<Particle> neutrino;
  
  cout << "\n--------------------------------------------------------------------------"<< endl;
  cout << "Run: " << run << " event: " << event << endl;
  
  foreach(const Particle &gen, *genParticles){
    if((abs(gen.id()) != 11 && abs(gen.id()) != 13 && abs(gen.id()) != 12 && abs(gen.id()) != 14) || (!(gen.genStatusFlags().test(phys::GenStatusBit::isPrompt)) || !(gen.genStatusFlags().test(phys::GenStatusBit::fromHardProcess)))) continue;
    cout << "id: " << gen.id() << " pt: " << gen.pt() << "\t eta: " << gen.eta() << endl;
    theHistograms.fill("ptAllGenParticle",  "p_{t} all particles", 100,  0  , 100  , gen.pt());
    theHistograms.fill("massAllGenParticle","mass all particles",  100,  0  , 100  , gen.mass());
    theHistograms.fill("etaAllGenParticle", "#eta all particles",  100,  0  , 100  , gen.eta());
    theHistograms.fill("YAllGenParticle",   "Y all particles",     100,  0  , 100  , gen.eta());
    theHistograms.fill("idAllGenParticle",  "ids all particles",     6, 10.5,  15.5, abs(gen.id()));
       
    if(gen.id() == 11)       electron.insert(electron.begin(), gen);
    else if(gen.id() == -11) electron.push_back(gen);
    else if(gen.id() == 13)  muon.insert(muon.begin(), gen);
    else if(gen.id() == -13) muon.push_back(gen);
    else if(abs(gen.id()) == 12 || abs(gen.id()) == 14)  neutrino.push_back(gen);
  }

  // ~~~~~~ tests on ZZ Analysis ~~~~~~
  //
  /*
  if(electron.size()+muon.size()!=4) {
    cout << "There are not enough or too many final leptons in this event." << endl;
    return;
  }

  totalevent++;

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
  
  //
  */
  // ~~~~~~ end of tests on ZZ Analysis ~~~~~~
 

  // ~~~~~~ WZ Analysis ~~~~~~
  //  /*
  if(electron.size()+muon.size()!=3)  {
    cout << "\nThere are not enough or too many final leptons in this event." << endl;
    return;
  }

  if(neutrino.size() != 1){
    cout << "\nThere are not enough or too many final neutrinos in this event." << endl;
    return;
  }

  nunumber++;
  
  Ztype Weh;
  Ztype Zet;
  vector<Particle> lepton;
  vector<Ztype> possibleW;
  vector<Ztype> possibleZ;
  vector<Zltype> Zls;
  
  // ------ filters on leptons ------
  // first filter on pt and eta
  foreach(const Particle ele, electron){
    if(ele.pt() < 7 || abs(ele.eta()) > 2.5){
      cout << "\nElectrons: pt less than 7 GeV or eta's absolute value greater than 2.5" << endl;
      return;
    }
    lepton.push_back(ele);
  }
  
  foreach(const Particle mu, muon){
    if(mu.pt() < 5 || abs(mu.eta()) > 2.4){
      cout << "\nMuons: pt less than 5 GeV or eta's absolute value greater than 2.4" << endl;
      return;
    }
    lepton.push_back(mu);
  }
  
  stable_sort(electron.begin(), electron.end(), PtComparator());
  stable_sort(lepton.begin(), lepton.end(), PtComparator());
  stable_sort(muon.begin(), muon.end(), PtComparator());

  // second filter on pt and eta -> unnecessary at the moment
  /*
  if(lepton[0].pt() < 20){
    cout << "\nFirst lepton pt less than 20 GeV" << endl;
    return;
  }
  if(abs(lepton[1].id()) == 11 && lepton[1].pt() < 12){
    cout << "\nSecond lepton is an electron and has pt less than 12 GeV" << endl;
    return;
  }
  if(abs(lepton[1].id()) == 13 && lepton[1].pt() < 10){
    cout << "\nSecond lepton is a muon and has pt less than 10 GeV" << endl;
    return;
  }
  */

  totalevent++;
  
  // Z and W must be on shell
  TLorentzVector Ptot = neutrino[0].p4();
  foreach(const Particle l, lepton)
    Ptot += l.p4();

  double masslllnu = Ptot.M();
  double trmasslllnu = Ptot.Mt();
  
  theHistograms.fill("allmasslllnu", "m 3 leptons and #nu", 1200, 0, 1200, masslllnu);
  theHistograms.fill("alltrmasslllnu", "m_{T} 3 leptons and #nu", 1200, 0, 1200, trmasslllnu);

  if(masslllnu < 165){
    cout << "\nTotal mass of the products insufficient for the WZ analysis." << endl;
    return;
  }

  // ------ Z & W ------
  WZevent++;
  
  theHistograms.fill("allmassWZ", "m 3 leptons and #nu", 1200, 160, 1360, masslllnu);
  
  if(electron.size()==2 && muon.size()==1){
    Zet = Ztype(electron[0], electron[1], 23);
    Zls.push_back(Zltype(Zet, muon[0]));
    
    if(muon[0].charge() > 0)
      Weh = Ztype(muon[0], neutrino[0], 24);
    else if(muon[0].charge() < 0)
      Weh = Ztype(muon[0], neutrino[0], -24);
  }
  
  else if(electron.size()==1 && muon.size()==2){
    Zet = Ztype(muon[0], muon[1], 23);
    Zls.push_back(Zltype(Zet, electron[0]));
    
    if(electron[0].charge() > 0)
      Weh = Ztype(electron[0], neutrino[0], 24);
    else if(electron[0].charge() < 0)
      Weh = Ztype(electron[0], neutrino[0], -24);
  }
    
  else if(electron.size()==3){
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

    stable_sort(possibleZ.begin(), possibleZ.end(), MassComparator(ZMASS));
    Zet = possibleZ[0];

    bool isSortOk = abs(possibleZ[0].mass() - ZMASS) < abs(possibleZ[1].mass() - ZMASS);
    theHistograms.fill("SortOK", "Is Sort OK", 2, -0.5, 1.5, isSortOk);
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

    stable_sort(possibleZ.begin(), possibleZ.end(), MassComparator(ZMASS));
    Zet = possibleZ[0];  

    bool isSortOk = abs(possibleZ[0].mass() - ZMASS) < abs(possibleZ[1].mass() - ZMASS);
    theHistograms.fill("SortOK", "Is Sort OK", 2, -0.5, 1.5, isSortOk);
  }
  
  cout << "\nZl candidates are: " << Zls.size() << endl;
  foreach(const Zltype zl, Zls){
    cout << " 1. Z " << zl.first << "\n 2. l " << zl.second << endl << endl;
  }
    
  if( ( isTheSame(Zet.daughter(0), Zls[0].first.daughter(0)) && isTheSame(Zet.daughter(1), Zls[0].first.daughter(1)) ) || ( isTheSame(Zet.daughter(0), Zls[0].first.daughter(1)) && isTheSame(Zet.daughter(1), Zls[0].first.daughter(0)) ) ){
    if(Zls[0].second.id() < 0)
      Weh = Ztype(Zls[0].second, neutrino[0], -24);
    if(Zls[0].second.id() > 0)
      Weh = Ztype(Zls[0].second, neutrino[0], 24);
    cout << " The best Z is in the first Zl couple." << endl;
  }
  
  if( ( isTheSame(Zet.daughter(0), Zls[1].first.daughter(0)) && isTheSame(Zet.daughter(1), Zls[1].first.daughter(1)) ) || ( isTheSame(Zet.daughter(0), Zls[1].first.daughter(1)) && isTheSame(Zet.daughter(1), Zls[1].first.daughter(0)) ) ){
    if(Zls[1].second.id() < 0)
      Weh = Ztype(Zls[1].second, neutrino[0], -24);
    if(Zls[1].second.id() > 0)
      Weh = Ztype(Zls[1].second, neutrino[0], 24);
    cout << " The best Z is in the second Zl couple." << endl;
  }
  
  cout << "Z is: " << Zet << endl;
  cout << "W is: " << Weh << "\n  her lepton daughter is: " << Weh.daughter(0) << endl;

  //To do: add histograms to study the distribution of Z and W (pt, eta, rapidity, transverse mass... ; their ids and their daughter's ids); make the prints out readable
  // */
   
}
  
void WZAnalyzer::end(TFile &){

  cout << "\n--------------------------------------------------------------------------"<< endl;
  
  cout << "\nNumber of events of the sample:                         " << setw(7) << zahl << endl;
  cout << "Number of events ending with 3 leptons and 1 neutrino:  " << setw(7) << nunumber << endl;
  cout << "Number of events useful for Wlllnu and WZ analysis:     " << setw(7) << totalevent << endl;
  cout << "Number of events useful for WZ analysis:                " << setw(7) << WZevent << endl;
  
  // execution time
  endtime = ((float)clock())/CLOCKS_PER_SEC;
  cout << "\nExecution time: " << (int)((endtime - begintime)/3600) << " h " << (((int)(endtime - begintime)%3600)/60) << " m " << endtime - begintime - (int)((endtime - begintime)/3600)*3600 - (((int)(endtime - begintime)%3600)/60)*60 << " s." << endl;
  cout << "\n--------------------------------------------------------------------------"<<endl;
}
