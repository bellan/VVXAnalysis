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

  threemuonsplus = 0;
  threemuonsminus = 0;
  threeelesplus = 0;
  threeelesminus = 0;

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
    //cout << "id: " << gen.id() << " pt: " << gen.pt() << "\t eta: " << gen.eta() << endl;
    theHistograms.fill("AllGenParticleid",  "ids all particles",     4, 10.5,  14.5, abs(gen.id()));
    theHistograms.fill("AllGenParticlept",  "p_{t} all particles", 200,  0  , 200  , gen.pt());
    theHistograms.fill("AllGenParticleY",   "Y all particles",     100,-10  ,  10  , gen.rapidity());
    theHistograms.fill("AllGenParticleeta", "#eta all particles",  100,-10  ,  10  , gen.eta());
       
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
  
  vector<Vtype> Zet;
  vector<Vtype> possibleZ;

  if(electron.size()==2 && muon.size()==2){
    possibleZ.push_back(Vtype(electron[0], electron[1], 23));
    possibleZ.push_back(Vtype(muon[0], muon[1], 23));
  }
  else if(electron.size()==4){
    possibleZ.push_back(Vtype(electron[0], electron[2], 23));
    possibleZ.push_back(Vtype(electron[0], electron[3], 23));
    possibleZ.push_back(Vtype(electron[1], electron[2], 23));
    possibleZ.push_back(Vtype(electron[1], electron[3], 23));
  }
  else if(muon.size()==4){      
    possibleZ.push_back(Vtype(muon[0], muon[2], 23));
    possibleZ.push_back(Vtype(muon[0], muon[3], 23));
    possibleZ.push_back(Vtype(muon[1], muon[2], 23));
    possibleZ.push_back(Vtype(muon[1], muon[3], 23));
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
    cout << Red("\nThere are not enough or too many final leptons in this event.") << endl;
    return;
  }

  if(neutrino.size() != 1){
    cout << Red("\nThere are not enough or too many final neutrinos in this event.") << endl;
    return;
  }

  if(electron.size() == 3){
    if(electron[0].charge() == electron[1].charge() && electron[1].charge() == electron[2].charge() && electron[0].charge() > 0){
      cout << Red("\nThere are three positrons.") << endl;
      threeelesplus++;
      return;
    }
    else {if(electron[0].charge() == electron[1].charge() && electron[1].charge() == electron[2].charge() && electron[0].charge() < 0){
	cout << Red("\nThere are three electrons.") << endl;
	threeelesminus++;
	return;
      }
    }
  }

  if(muon.size() == 3){
    if(muon[0].charge() == muon[1].charge() && muon[1].charge() == muon[2].charge() && muon[0].charge() > 0){
      cout << Red("\nThere are three antimuons.") << endl;
      threemuonsplus++;
      return;
    }
    else {if(muon[0].charge() == muon[1].charge() && muon[1].charge() == muon[2].charge() && muon[0].charge() < 0){
	cout << Red("\nThere are three muons.") << endl;
	threemuonsminus++;
	return;
      }
    }
  }
  
  nunumber++;
  
  Vtype Weh;
  Vtype Zet;
  vector<Particle> lepton;
  vector<Vtype> possibleW;
  vector<Vtype> possibleZ;
  vector<Zltype> Zls;
  
  // ------ filters on leptons ------
  // first filter on pt and eta
  foreach(const Particle ele, electron){
    if(ele.pt() < 7 || abs(ele.eta()) > 2.5){
      cout << Violet("\nElectrons: pt less than 7 GeV or eta's absolute value greater than 2.5") << endl;
      return;
    }
    lepton.push_back(ele);
    theHistograms.fill("Lcharge", "leptons charge", 2, -2, 2, ele.charge());
  }
  
  foreach(const Particle mu, muon){
    if(mu.pt() < 5 || abs(mu.eta()) > 2.4){
      cout << Violet("\nMuons: pt less than 5 GeV or eta's absolute value greater than 2.4") << endl;
      return;
    }
    lepton.push_back(mu);
    theHistograms.fill("Lcharge", "leptons charge", 2, -2, 2, mu.charge());
  }
  
  stable_sort(electron.begin(), electron.end(), PtComparator());
  stable_sort(lepton.begin(), lepton.end(), PtComparator());
  stable_sort(muon.begin(), muon.end(), PtComparator());

  // second filter on pt and eta
  // /*
  if(lepton[0].pt() < 20){
    cout << Violet("\nFirst lepton pt less than 20 GeV") << endl;
    return;
  }
  if(abs(lepton[1].id()) == 11 && lepton[1].pt() < 12){
    cout << Violet("\nSecond lepton is an electron and has pt less than 12 GeV") << endl;
    return;
  }
  if(abs(lepton[1].id()) == 13 && lepton[1].pt() < 10){
    cout << Violet("\nSecond lepton is a muon and has pt less than 10 GeV") << endl;
    return;
  }
  // */

  totalevent++;
  
  // Z and W must be on shell
  TLorentzVector Ptot = neutrino[0].p4();
  foreach(const Particle l, lepton)
    Ptot += l.p4();

  double masslllnu = Ptot.M();
  double trmasslllnu = Ptot.Mt();

  theHistograms.fill("allmasslllnu", "m 3 leptons and #nu", 1200, 0, 1200, masslllnu); //what happens between 90 and 160 GeV?
  theHistograms.fill("alltrmasslllnu", "m_{T} 3 leptons and #nu", 1200, 0, 1200, trmasslllnu);

  if(masslllnu < 165){
    cout << Yellow("\nTotal mass of the products insufficient for the WZ analysis.") << endl;
    return;
  }

  foreach(const Particle &gen, *genParticles){
    if((abs(gen.id()) != 11 && abs(gen.id()) != 13 && abs(gen.id()) != 12 && abs(gen.id()) != 14) || (!(gen.genStatusFlags().test(phys::GenStatusBit::isPrompt)) || !(gen.genStatusFlags().test(phys::GenStatusBit::fromHardProcess)))) continue;
    cout << "id: " << gen.id() << " pt: " << gen.pt() << "\t eta: " << gen.eta() << endl;
  }
  
  // ------ Z & W ------
  WZevent++;

  // Construction of the two possible Zs
  if(electron.size()==2 && muon.size()==1 && electron[0].charge() != electron[1].charge()){
    Zet = Vtype(electron[0], electron[1], 23);
    Zls.push_back(Zltype(Zet, muon[0]));
  }
  
  else if(electron.size()==1 && muon.size()==2 && muon[0].charge() != muon[1].charge()){
    Zet = Vtype(muon[0], muon[1], 23);
    Zls.push_back(Zltype(Zet, electron[0]));
  }
    
  else if(electron.size()==3){
    if(electron[0].id() != electron[2].id()){
      possibleZ.push_back(Vtype(electron[0], electron[2], 23));
      Zls.push_back(Zltype(possibleZ[0], electron[1]));
    }
      
    if(electron[0].id() == electron[1].id() && electron[1].id() != electron[2].id()){
      possibleZ.push_back(Vtype(electron[1], electron[2], 23));
      Zls.push_back(Zltype(possibleZ[1], electron[0]));
    }
    else{if(electron[0].id() != electron[1].id()){
	possibleZ.push_back(Vtype(electron[0], electron[1], 23));
	Zls.push_back(Zltype(possibleZ[1], electron[2]));
      }
    }
      
    // Z is made up of the couple which gives a better Zmass 
    stable_sort(possibleZ.begin(), possibleZ.end(), MassComparator(ZMASS));
    Zet = possibleZ[0];

    /*
      bool isSortOk = abs(possibleZ[0].mass() - ZMASS) < abs(possibleZ[1].mass() - ZMASS);
      theHistograms.fill("SortOK", "Is Sort OK", 2, -0.5, 1.5, isSortOk);
    */
  }
  
  else if(muon.size()==3){
    if(muon[0].id() != muon[2].id()){
      possibleZ.push_back(Vtype(muon[0], muon[2], 23));
      Zls.push_back(Zltype(possibleZ[0], muon[1]));
    }
      
    if(muon[0].id() == muon[1].id() && muon[1].id() != muon[2].id()){
      possibleZ.push_back(Vtype(muon[1], muon[2], 23));
      Zls.push_back(Zltype(possibleZ[1], muon[0]));
    }
    else{if(muon[0].id() != muon[1].id()){
	possibleZ.push_back(Vtype(muon[0], muon[1], 23));
	Zls.push_back(Zltype(possibleZ[1], muon[2]));
      }
    }
    
    stable_sort(possibleZ.begin(), possibleZ.end(), MassComparator(ZMASS));
    Zet = possibleZ[0];  

    /*    
      bool isSortOk = abs(possibleZ[0].mass() - ZMASS) < abs(possibleZ[1].mass() - ZMASS);
      theHistograms.fill("SortOK", "Is Sort OK", 2, -0.5, 1.5, isSortOk);
    */
  }
  
  cout << "\nZl candidates are: " << Zls.size() << endl;
  
  foreach(const Zltype zl, Zls){
    cout << "   Z " << zl.first << "\n   l " << zl.second << endl << endl;
    
    // W is made up of the remaining lepton and the neutrino
    if( ( isTheSame(Zet.daughter(0), zl.first.daughter(0)) && isTheSame(Zet.daughter(1), zl.first.daughter(1)) ) || ( isTheSame(Zet.daughter(0), zl.first.daughter(1)) && isTheSame(Zet.daughter(1), zl.first.daughter(0)) ) ){
      Weh = Vtype(zl.second, neutrino[0], copysign(24, zl.second.charge()));
    }
  }
  
  cout << "Z is: " << Zet << endl;
  cout << "W is: " << Weh << "\n  her lepton daughter is: " << Weh.daughter(0) << endl;
  
  //W histograms
  theHistograms.fill("Wcharge", "W's charge",   3,-1.5, 1.5, Weh.charge()); //why are there more W+ than W-?
  theHistograms.fill("Wmass",   "W's mass",   350,   0, 350, Weh.mass());
  theHistograms.fill("Wpt"  ,   "W's p_{t}",  300,   0, 600, Weh.pt());
  theHistograms.fill("WY",      "W's Y",       50,  -5,   5, Weh.rapidity());
  theHistograms.fill("Weta",    "W's #eta",    50,  -9,   9, Weh.eta());
  
  //Z histograms
  theHistograms.fill("Zcharge", "Z's charge", 5, -2.5, 2.5, Zet.charge()); //just to be sure...
  theHistograms.fill("Zmass",   "Z's mass",  350,   0, 350, Zet.mass());
  theHistograms.fill("Zpt",     "Z's p_{t}", 300,   0, 600, Zet.pt());
  theHistograms.fill("ZY",      "Z's Y",      50,  -4,   4, Zet.rapidity());
  theHistograms.fill("Zeta",    "Z's #eta",  100,  -7,   7, Zet.eta());

  //W&Z histograms
  theHistograms.fill("allmassWZ", "m 3 leptons and #nu", 1200, 160, 1360, masslllnu);
    
  bool areZWOnShell = Zet.mass() >= 86 && Zet.mass() <= 96 && Weh.mass() >= 75 && Weh.mass() <= 85;
  theHistograms.fill("ZWOnShell", "Are W and Z on shell?", 2, -0.5, 1.5, areZWOnShell);
  
  if(areZWOnShell == 0){
    int cases = 0;
    if(Zet.mass() < 86 && Weh.mass() < 75)                            // Z sotto  W sotto
      cases = 1;
    else if(Zet.mass() < 86 && Weh.mass() > 85)                       // Z sotto  W sopra
      cases = 2;
    else if(Zet.mass() > 96 && Weh.mass() < 75)                       // Z sopra  W sotto 
      cases = 3;
    else if(Zet.mass() > 96 && Weh.mass() > 85)                       // Z sopra  W sopra
      cases = 4;
    else if(Zet.mass() >= 86 && Zet.mass() <= 96 && Weh.mass() < 75)  // Z dentro W sotto
      cases = 5;
    else if(Zet.mass() >= 86 && Zet.mass() <= 96 && Weh.mass() > 85)  // Z dentro W sopra
      cases = 6;
    else if(Weh.mass() >= 75 && Weh.mass() <= 85 && Zet.mass() < 86)  // Z sotto  W dentro
      cases = 7;
    else if(Weh.mass() >= 75 && Weh.mass() <= 85 && Zet.mass() > 96)  // Z sopra  W dentro
      cases = 8;

    theHistograms.fill("ZWcases", "ZWcases", 9, -0.5, 8.5, cases);
  }
  
  // */
   
}
  
void WZAnalyzer::end(TFile &){

  cout << "\n--------------------------------------------------------------------------"<< endl;
  
  cout << "\nNumber of events of the sample:                         " << setw(7) << zahl << endl;
  cout << "Number of events ending with 3 leptons and 1 neutrino:  " << setw(7) << nunumber << endl;
  cout << "Number of events useful for Wlllnu and WZ analysis:     " << setw(7) << totalevent << endl;
  cout << "Number of events useful for WZ analysis:                " << setw(7) << WZevent << endl;

  /*
  cout << "\nNumber of events ending with 3 positrons                " << setw(7) << threeelesplus << endl;
  cout << "Number of events ending with 3 electrons                " << setw(7) << threeelesminus << endl;
  cout << "Number of events ending with 3 antimuons                " << setw(7) << threemuonsplus << endl;
  cout << "Number of events ending with 3 muons                    " << setw(7) << threemuonsminus << endl;
  */
  
  // execution time
  endtime = ((float)clock())/CLOCKS_PER_SEC;
  WZAnalyzer::printTime(begintime, endtime);
  cout << "\n--------------------------------------------------------------------------"<<endl;
}
