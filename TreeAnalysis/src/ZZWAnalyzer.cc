#include "VVXAnalysis/TreeAnalysis/interface/ZZWAnalyzer.h"
                                                                                                                                                                                                                                                                                                                                                                             
#include <boost/foreach.hpp>
#define foreach BOOST_FOREACH

using namespace phys;
using namespace std;


//template< class PAR >
//bool ZZWAnalyzer::bosonDefinition(phys::Boson<PAR> vb) const{
//  return true;
//}

Int_t ZZWAnalyzer::cut() {
  
  theHistograms.fill<TH1I>("Number of events", "Number of events", 10, 0, 10, 0);  // Number of events without any extra cut

  float mZ = 90.19;

  bool passMass = true;
  bool passSize = Wjj->size() >= 1;  //// ???

  foreach(const Boson<Lepton>& zm, *Zmm) {
    foreach(const Boson<Electron>& ze, *Zee) {
      passMass = passMass && fabs(zm.p4().M() - mZ ) < 20. && fabs(ze.p4().M() - mZ ) < 20.;
    }
  }
 
  bool pass1 = passMass && passSize;

  if(pass1) {  
    theHistograms.fill("Number of events" , 10, 0, 10, 1);    //Number of events after the first cut: Zreco Mass range (mZ +/- 20 GeV) + At least 1 W  
  }                                                                   

  int numW = 0;
  foreach(const Boson<Jet>& w, *Wjj)
    if(w.daughter(0).pt() > 40 && w.daughter(1).pt() > 40) ++numW;       
  
  bool pass = pass1 && numW >= 1;

  if(pass) {
    ++theCutCounter;
    theHistograms.fill("Number of events" , 10, 0, 10, 2);   //Number of events after the second cut: Zreco Mass range (mZ +/- 20 GeV) + At least 1 W + 2jets Pt>40GeV 
  }

  return pass ? 1 : -1;
      
}


void ZZWAnalyzer::analyze() {
  
  //================================ genParticles ================================

  int numMu = 0;
  int numE  = 0;
  
  std::vector<const Particle* > Genq;
  std::vector<const Particle* > Genj;
  std::vector<const Particle* > Genl;    
  std::vector<const Particle* > Genlp;
  std::vector<const Particle* > Genlm;
  std::vector<const Particle* > GenZ;
  std::vector<const Particle* > GenW;
  
  foreach(const Particle& p, *genParticles) {
    
    int s_id = p.id();
    int id   = abs(s_id);
        
    if ( (id < 8 || id == 9)) {              // quarks and gluons
      Genj.push_back(&p);
    } 
    if ( id < 8) {                           // quarks
      Genq.push_back(&p);
    } else if ( id==23 ) {                   // Z
      GenZ.push_back(&p);
    } else if ( id==24 ) {                   // W
      GenW.push_back(&p);
    } else if ( id >= 11 && id <= 16 ) {     // leptons 
      numE  = id == 11 ? numE+1  : numE;
      numMu = id == 13 ? numMu+1 : numMu;
      Genl.push_back(&p);
      if(s_id>0) {Genlp.push_back(&p);
      } else {Genlm.push_back(&p);}
    }
    
  }

  const Particle* Z0gen = GenZ.at(0);
  const Particle* Z1gen = GenZ.at(1);
  const Particle* Wgen  = GenW.at(0);


  cout << "=============== genParticles: history information ===============" << endl; 

  cout << "Z0genMass= " << Z0gen->p4().M() << endl;
  cout << "Z1genMass= " << Z1gen->p4().M() << endl;
  cout << "WgenMass= "  << Wgen->p4().M() << endl;
  cout << "------------------------" << endl;
  cout << "Number of partons in the jets= " << Genj.size() << endl;
  cout << "Number of leptons= "             << Genl.size() << endl;


  /////-----------------Histograms-------------------
  
  //------------Mass--------------
  
  theHistograms.fill("Z0Gen_Mass", 200, 0, 200, Z0gen->p4().M(), theWeight);
  theHistograms.fill("Z1Gen_Mass", 200, 0, 200, Z1gen->p4().M(), theWeight);
  theHistograms.fill("WGen_Mass" , 200, 0, 200, Wgen->p4().M() , theWeight);
  
  //------------Pt--------------
  
  theHistograms.fill("Z0Gen_Pt", 300, 0, 300, Z0gen->pt(), theWeight);
  theHistograms.fill("Z1Gen_Pt", 300, 0, 300, Z1gen->pt(), theWeight);
  theHistograms.fill("WGen_Pt" , 300, 0, 300, Wgen->pt() , theWeight);
  
  //----------------------------

  theHistograms.fill<TH1I>("Number of leptons", "Number of leptons", 10, 0, 10, Genl.size());
  theHistograms.fill<TH1I>("Number of partons", "Number of partons", 10, 0, 10, Genj.size());


  //================================ reco Particles ================================
  
  const Particle* Z0;
  const Particle* Z1;
  Boson<Jet> W;
  
  std::vector<const Particle* > Zll;
  
  foreach(const Boson<Lepton>& z, *Zmm)
    Zll.push_back(&z);
  
  foreach(const Boson<Electron>& z, *Zee)
    Zll.push_back(&z);

  
  std::stable_sort(Zll.begin(),Zll.end(),MassComparator(ZMASS));
  std::stable_sort(Wjj->begin(),Wjj->end(),MassComparator(WMASS));

  Z0 = Zll.at(0);
  Z1 = Zll.at(1);
  W  = Wjj->at(0);
  
  TLorentzVector p_Z0 = Z0->p4();
  TLorentzVector p_Z1 = Z1->p4();
  
  TLorentzVector p_j1 = W.daughter(0).p4();
  TLorentzVector p_j2 = W.daughter(1).p4();
  
  
  /////-----------------Histograms-------------------
  
  //------------Mass--------------
  
  theHistograms.fill("Wjj_Mass", 200,0,200,W.p4().M(), theWeight);  
  theHistograms.fill("Z0_Mass",200,0,200,p_Z0.M()    , theWeight);
  theHistograms.fill("Z1_Mass",200,0,200,p_Z1.M()    , theWeight);
  
  //------------Pt--------------
  
  theHistograms.fill("Z0_Pt", 300,0,300,Z0->pt(), theWeight);
  theHistograms.fill("Z1_Pt", 300,0,300,Z1->pt(), theWeight);
  theHistograms.fill("W_Pt", 300,0,300,W.pt()   , theWeight);  
  
  //------------Mass 6f--------------
  
  TLorentzVector p_6f = p_Z0 + p_Z1 + p_j1 + p_j2;
  
  theHistograms.fill("6f_Mass",3000,0,3000,p_6f.M(), theWeight);
  
}
