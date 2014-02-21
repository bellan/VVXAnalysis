#include "VVXAnalysis/TreeAnalysis/interface/ZZWAnalyzer.h"
                                                                                                                                                                                                                                                                                                                                                                             
#include <boost/foreach.hpp>
#define foreach BOOST_FOREACH

using namespace phys;
using namespace std;

Int_t ZZWAnalyzer::cut(){
  
  bool passZSize = (Zmm->size() + Zee->size()) >= 2;

  bool pass = true;
  
  bool passSize = passZSize && Wjj->size() >= 1;
  
  int numW = 0;
  foreach(const Boson<Jet>& w, *Wjj)
    if(w.daughter(0).pt() > 40 && w.daughter(1).pt() > 40) ++numW;
  
  pass = passSize && numW >=1;
  
  if(pass) ++theCutCounter;
  
  return pass ? 1 : -1;
  
}


void ZZWAnalyzer::analyze() {
  
  cout << "Event " << event << endl;
  
  const Particle* Z0;
  const Particle* Z1;
  Boson<Jet> W;
  
  std::vector<const Particle* > Zll;
  
  foreach(const Boson<Lepton>& z, *Zmm)
    Zll.push_back(&z);
  
  foreach(const Boson<Electron>& z, *Zee)
    Zll.push_back(&z);
  
  std::stable_sort(Zll.begin(),Zll.end(),MassComparator(ZMASS));
    
  Z0 = Zll.at(0);
  Z1 = Zll.at(1);
  
  std::stable_sort(Wjj->begin(),Wjj->end(),MassComparator(WMASS));
  
  W = Wjj->at(0);
  
  TLorentzVector p_Z0 = Z0->p4();
  TLorentzVector p_Z1 = Z1->p4();
  
  TLorentzVector p_j1 = W.daughter(0).p4();
  TLorentzVector p_j2 = W.daughter(1).p4();
  
  cout << "Z0_Mass= " << p_Z0.M() << endl;
  cout << "Z0_Mass= " << p_Z1.M() << endl;
  cout << "W_Mass= "  << W.p4().M() << endl;
 
  //================================Histograms=====================================
  
  //------------Mass--------------
  
  theHistograms.fill("Wjj_Mass", 200,0,200,W.p4().M());
  
  theHistograms.fill("Z0_Mass",200,0,200,p_Z0.M());
  theHistograms.fill("Z1_Mass",200,0,200,p_Z1.M());
  
  //------------Pt--------------
  
  theHistograms.fill("Z0_Pt", 100,0,100,Z0->pt());
  theHistograms.fill("Z1_Pt", 100,0,100,Z1->pt());
  theHistograms.fill("W_Pt", 100,0,100,W.pt());  
  
  //------------Mass 6f--------------
  
  TLorentzVector p_6f = p_Z0 + p_Z1 + p_j1 + p_j2;
  
  theHistograms.fill("6f_Mass",3000,0,3000,p_6f.M());
  
}
