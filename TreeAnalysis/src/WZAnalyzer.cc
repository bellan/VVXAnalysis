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
  vector<Particle> jets;
  
  cout << "\n--------------------------------------------------------------------------"<< endl;
  cout << "Run: " << run << " event: " << event << endl;
  
  foreach(const Particle &gen, *genParticles){
    //if((abs(gen.id()) != 11 && abs(gen.id()) != 13 && abs(gen.id()) != 12 && abs(gen.id()) != 14) || (!(gen.genStatusFlags().test(phys::GenStatusBit::isPrompt)) || !(gen.genStatusFlags().test(phys::GenStatusBit::fromHardProcess)))) continue;
    if((!(gen.genStatusFlags().test(phys::GenStatusBit::isPrompt)) || !(gen.genStatusFlags().test(phys::GenStatusBit::fromHardProcess)))) continue;
    //cout << "id: " << gen.id() << " pt: " << gen.pt() << "\t eta: " << gen.eta() << endl;
    theHistograms.fill("AllGenParticleid",  "ids all particles",     4, 10.5,  14.5, abs(gen.id()));
    theHistograms.fill("AllGenParticlept",  "p_{t} all particles", 200,  0  , 200  , gen.pt());
    theHistograms.fill("AllGenParticleY",   "Y all particles",     100,-10  ,  10  , gen.rapidity());
    theHistograms.fill("AllGenParticleeta", "#eta all particles",  100,-10  ,  10  , gen.eta());
    
    if(abs(gen.id()) == 11)      electron.push_back(gen);
    else if(abs(gen.id()) == 13) muon.push_back(gen);
    else if(abs(gen.id()) == 12 || abs(gen.id()) == 14)  neutrino.push_back(gen);
  } 

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
  }
  
  foreach(const Particle mu, muon){
    if(mu.pt() < 5 || abs(mu.eta()) > 2.4){
      cout << Violet("\nMuons: pt less than 5 GeV or eta's absolute value greater than 2.4") << endl;
      return;
    }
    lepton.push_back(mu);
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

  foreach(const Particle ele, electron){
    theHistograms.fill("Lcharge", "leptons charge", 3, -1.5, 1.5, ele.charge());
  }

  foreach(const Particle mu, muon){
    theHistograms.fill("Lcharge", "leptons charge", 3, -1.5, 1.5, mu.charge());
  }
  
  totalevent++;
  
  // Z and W must be on shell
  TLorentzVector Ptot = neutrino[0].p4();
  foreach(const Particle l, lepton)
    Ptot += l.p4();

  double masslllnu = Ptot.M();
  double trmasslllnu = Ptot.Mt();

  theHistograms.fill("allmasslllnu", "m 3 leptons and #nu", 1500, 0, 1500, masslllnu); //what happens between 90 and 160 GeV?
  theHistograms.fill("alltrmasslllnu", "m_{T} 3 leptons and #nu", 1500, 0, 1500, trmasslllnu);

  if(masslllnu < 165){
    cout << Yellow("\nTotal mass of the products insufficient for the WZ analysis.") << endl;
    return;
  }

  foreach(const Particle &gen, *genParticles){
    if((abs(gen.id()) != 11 && abs(gen.id()) != 13 && abs(gen.id()) != 12 && abs(gen.id()) != 14) || (!(gen.genStatusFlags().test(phys::GenStatusBit::isPrompt)) || !(gen.genStatusFlags().test(phys::GenStatusBit::fromHardProcess)))) continue;
    cout << "id: " << gen.id() << " pt: " << setw(5) << gen.pt() << "\t eta: " << gen.eta() << endl;
  }
    
  // ------ Z & W ------
  WZevent++;

  // Construction of the two possible Zs
  
  if(electron.size()==2 && muon.size()==1 && electron[0].charge() != electron[1].charge()){
    //return;
    // /*
    Zet = Vtype(electron[0], electron[1], 23);
    Zls.push_back(Zltype(Zet, muon[0]));
    Weh = Vtype(Zls[0].second, neutrino[0], copysign(24, Zls[0].second.charge()) );
    // */
  }
  
  if(electron.size()==1 && muon.size()==2 && muon[0].charge() != muon[1].charge()){
    //return;
    // /*
    Zet = Vtype(muon[0], muon[1], 23);
    Zls.push_back(Zltype(Zet, electron[0]));
    Weh = Vtype(Zls[0].second, neutrino[0], copysign(24, Zls[0].second.charge()) );
    // */
  }
  
  if(electron.size()==3){
    //return;
    ///*
    if(electron[0].charge() != electron[2].charge()){
      possibleZ.push_back(Vtype(electron[0], electron[2], 23));
      Zls.push_back(Zltype(possibleZ.back(), electron[1]));
    }
    
    if(electron[1].charge() != electron[2].charge()){
      possibleZ.push_back(Vtype(electron[1], electron[2], 23));
      Zls.push_back(Zltype(possibleZ.back(), electron[0]));
    }

    if(electron[0].charge() != electron[1].charge()){
      possibleZ.push_back(Vtype(electron[0], electron[1], 23));
      Zls.push_back(Zltype(possibleZ.back(), electron[2]));
    }
    //*/
    /*
    int k = electron.size() - 1;
    for(int i = 0; i < (int)electron.size() -2; i++){
      for(int j = 1; j < (int)electron.size() -1; j++){
	if(electron[i].charge() != electron[j].charge()){
	  possibleZ.push_back(Vtype(electron[i], electron[j], 23));
	  Zls.push_back(Zltype (possibleZ.back(), electron[k]));
	}
	k--;
      }
    }
    */
    
    // Z is made up of the couple which gives a better Zmass 
    if(Zls.size() < 1){
      cout << Red("No Z formed.") << endl;
      threeelesplus++;
      return;
    }
    
    if(Zls.size() > 1){
      stable_sort(Zls.begin(), Zls.end(), ZlMassComparator(ZMASS));
    }

    Zet = Zls[0].first;

    // W is made up of the remaining lepton and the neutrino
    Weh = Vtype(Zls[0].second, neutrino[0], copysign(24, Zls[0].second.charge()));
    // */
    /*
      bool isSortOk = abs(possibleZ[0].mass() - ZMASS) < abs(possibleZ[1].mass() - ZMASS);
      theHistograms.fill("SortOK", "Is Sort OK", 2, -0.5, 1.5, isSortOk);
    */
  }
  
  else if(muon.size()==3){
    //return;
    ///*
    if(muon[0].charge() != muon[2].charge()){
      possibleZ.push_back(Vtype(muon[0], muon[2], 23));
      Zls.push_back(Zltype(possibleZ.back(), muon[1]));
    }
    
    if(muon[1].charge() != muon[2].charge()){
      possibleZ.push_back(Vtype(muon[1], muon[2], 23));
      Zls.push_back(Zltype(possibleZ.back(), muon[0]));
    }

    if(muon[0].charge() != muon[1].charge()){
      possibleZ.push_back(Vtype(muon[0], muon[1], 23));
      Zls.push_back(Zltype(possibleZ.back(), muon[2]));
    }
    //*/
    /*
    int k = muon.size() - 1;
    for(int i = 0; i < (int)muon.size() -2; i++){
      for(int j = 1; j < (int)muon.size() -1; j++){
	if(muon[i].charge() != muon[j].charge()){
	  possibleZ.push_back(Vtype(muon[i], muon[j], 23));
	  Zls.push_back(Zltype (possibleZ.back(), muon[k]));
	}
	k--;
      }
    }
    */

    if(Zls.size() < 1){
      cout << Red("No Z formed.") << endl;
      threemuonsplus++;
      return;
    }
    
    if(Zls.size() > 1){
      stable_sort(Zls.begin(), Zls.end(), ZlMassComparator(ZMASS));
    }

    Zet = Zls[0].first;
    Weh = Vtype(Zls[0].second, neutrino[0], copysign(24, Zls[0].second.charge()));
    // */
    /*    
      bool isSortOk = abs(possibleZ[0].mass() - ZMASS) < abs(possibleZ[1].mass() - ZMASS);
      theHistograms.fill("SortOK", "Is Sort OK", 2, -0.5, 1.5, isSortOk);
    */
  }
  
  cout << "\nZl candidates are: " << Zls.size() << endl;
  
  foreach(const Zltype zl, Zls){
    cout << "   Z " << zl.first << "\n   l " << zl.second << endl << endl;
  }
  
  cout << "Z is: " << Zet << endl;
  cout << "W is: " << Weh << "\n  her lepton daughter is: " << Weh.daughter(0) << endl;
  
  //W histograms
  theHistograms.fill("Wcharge", "W's charge",   5,-2.5, 2.5, Weh.charge());
  theHistograms.fill("Wmass",   "W's mass",   400,   0, 400, Weh.mass());
  theHistograms.fill("Wpt"  ,   "W's p_{t}",  650,   0, 650, Weh.pt());
  theHistograms.fill("WY",      "W's Y",       50,  -5,   5, Weh.rapidity());
  theHistograms.fill("Weta",    "W's #eta",    50,  -9,   9, Weh.eta());
  
  //Z histograms
  theHistograms.fill("Zcharge", "Z's charge", 5, -2.5, 2.5, Zet.charge()); //just to be sure...
  theHistograms.fill("Zmass",   "Z's mass",  400,   0, 400, Zet.mass());
  theHistograms.fill("Zpt",     "Z's p_{t}", 650,   0, 650, Zet.pt());
  theHistograms.fill("ZY",      "Z's Y",      50,  -4,   4, Zet.rapidity());
  theHistograms.fill("Zeta",    "Z's #eta",  100,  -7,   7, Zet.eta());

  //W&Z histograms
  theHistograms.fill("allmassWZ", "m 3 leptons and #nu", 1340, 160, 1500, masslllnu);
    
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

  // ------- Jet  -------

  cout << "\n~~~~~~~~~~~~~~~~ Jets ~~~~~~~~~~~~~~~~\n" << endl;

  foreach(const Particle jet, *pgenJets){
    jets.push_back(jet);

    cout << "ID: " << jet.id() << " pt: " << jet.pt() << "\t eta: " << jet.eta() << endl;

    theHistograms.fill("AllGenJetsnumber", "number of jets",  13,  -0.5,  12.5, jets.size());
    theHistograms.fill("AllGenJetpt",      "p_{t} all jets", 500,   0  , 500  , jet.pt());
    theHistograms.fill("AllGenJetY",       "Y all jets",      80,  -6  ,   6  , jet.rapidity());
    theHistograms.fill("AllGenJeteta",     "#eta all jets",   80,  -6  ,   6  , jet.eta());    
    theHistograms.fill("AllGenJetcharge",  "charge all jets",  5,  -2.5,   2.5, jet.charge());
    theHistograms.fill("AllGenJetmass",    "mass all jets",  100,   0  , 100  , jet.mass());
  }
  
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
