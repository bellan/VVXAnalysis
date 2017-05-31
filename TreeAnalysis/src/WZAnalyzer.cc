#include "VVXAnalysis/TreeAnalysis/interface/WZAnalyzer.h"
#include "VVXAnalysis/Commons/interface/SignalDefinitions.h"
#include "VVXAnalysis/Commons/interface/Utilities.h"
#include "VVXAnalysis/Commons/interface/Comparators.h"
#include "VVXAnalysis/Commons/interface/Constants.h"
#include "VVXAnalysis/Commons/interface/AriEle.h"

#include <boost/foreach.hpp>
#define foreach BOOST_FOREACH

#include <boost/assign/std/vector.hpp> 

using namespace boost::assign;
using namespace phys;
using namespace std;

Int_t WZAnalyzer::cut() {
  return 1;
}


void WZAnalyzer::analyze(){

  std::vector<Particle> electron;
  std::vector<Particle> muon;
  
  cout << "------------------------------------------------------------------"<<endl;
  cout << "Run: " << run << " event: " << event << endl;
  
  foreach(const Particle &gen, *genParticles){
    if((abs(gen.id()) != 11 && abs(gen.id()) != 13) || (!(gen.genStatusFlags().test(phys::GenStatusBit::isPrompt)) || !(gen.genStatusFlags().test(phys::GenStatusBit::fromHardProcess)))) continue;
    //if(abs(gen.id()) != 11 && abs(gen.id()) != 13) continue;
    cout << "id: " << gen.id() << " pt: " << gen.pt() << endl;
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
    cout << "There are not enough or too many leptons in this event. They are: " << electron.size()+muon.size() << endl;
    return;
  }

  // Reconstruction of the two Zs
  
  vector<Boson<Particle> > Zet;
  vector<Boson<Particle> > possibleZ;

  if(electron.size()==2 && muon.size()==2){
    possibleZ.push_back(Boson<Particle>(electron[0], electron[1], 23));
    possibleZ.push_back(Boson<Particle>(muon[0], muon[1], 23));
  }
  
  else if(electron.size()==4){
    possibleZ.push_back(Boson<Particle>(electron[0], electron[2], 23));
    possibleZ.push_back(Boson<Particle>(electron[0], electron[3], 23));
    possibleZ.push_back(Boson<Particle>(electron[1], electron[2], 23));
    possibleZ.push_back(Boson<Particle>(electron[1], electron[3], 23));
  }
  else if(muon.size()==4){      
    possibleZ.push_back(Boson<Particle>(muon[0], muon[2], 23));
    possibleZ.push_back(Boson<Particle>(muon[0], muon[3], 23));
    possibleZ.push_back(Boson<Particle>(muon[1], muon[2], 23));
    possibleZ.push_back(Boson<Particle>(muon[1], muon[3], 23));
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
  theHistograms.fill("massZZ", "ZZ mass", 1000,  0, 400, doppelZ.pt());
  theHistograms.fill("YZZ",    "ZZ Y",     100,  0, 100, doppelZ.rapidity());
  
  cout << "Reconstructed ZZ: " << doppelZ << endl;
  

  // */
  // ~~~~~~ end of tests on ZZ Analysis ~~~~~~
 

  // ~~~~~~ WZ Analysis ~~~~~~
   /*
  if(electron.size()+muon.size()!=3) return;
  
   */
 
  
}
  
