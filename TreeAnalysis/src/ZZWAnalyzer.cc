#include "VVXAnalysis/TreeAnalysis/interface/ZZWAnalyzer.h"
#include "VVXAnalysis/TreeAnalysis/interface/Colours.h"
                                                                                                                                                                                                                                                                                                                                                                           
#include <boost/foreach.hpp>
#define foreach BOOST_FOREACH

using namespace phys;
using namespace std;
using namespace colour;

Int_t ZZWAnalyzer::cut() {
  
  bool category0 = genCategory == 0;

  bool passSize = ( Zmm->size() + Zee->size() ) >= 2;

  bool offset =  category0 && passSize; 

  if(offset) theHistograms.fill<TH1I>("Number of events", "Number of events", 10, 0, 10, 0);  // Number of total events --------------------------

 
  bool pass1 =  offset && Wjj->size() >= 1;

  if(pass1) theHistograms.fill("Number of events" , 10, 0, 10, 1);    //Number of events after the first cut: At least 1 W -----------------------                               


  int numW = 0;
  foreach(const Boson<Jet>& w, *Wjj)
    if(w.daughter(0).pt() > 40 && w.daughter(1).pt() > 40) ++numW;       
  
  bool pass = pass1;  //&& numW >= 1;; 

  if(pass) {
    ++theCutCounter;
    theHistograms.fill("Number of events" , 10, 0, 10, 2);   //Number of events after the second cut: 2jets Pt>40GeV -----------------------------
  }

  return pass ? 1 : -1;
      
}


void ZZWAnalyzer::analyze() {

  
  //================================================================//
  //                             genParticles                       //    
  //================================================================//
  

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
  cout << "Number of partons in the jets= " << Genj.size() << endl;
  cout << "Number of leptons= "             << Genl.size() << endl;
  cout << "Category= "                      << genCategory << endl;


  /////-----------------Histograms-------------------
  
  //------------Mass--------------
  
  theHistograms.fill("Z0Gen_Mass", "Z0Gen_Mass", 200, 0, 200, Z0gen->p4().M(), theWeight);
  theHistograms.fill("Z1Gen_Mass", "Z1Gen_Mass", 200, 0, 200, Z1gen->p4().M(), theWeight);
  theHistograms.fill("WGen_Mass" , "WGen_Mass" , 200, 0, 200, Wgen->p4().M() , theWeight);
  
  //------------Pt--------------
  
  theHistograms.fill("Z0Gen_Pt"  , "Z0Gen_Pt"  , 300, 0, 300, Z0gen->pt()    , theWeight);
  theHistograms.fill("Z1Gen_Pt"  , "Z1Gen_Pt"  , 300, 0, 300, Z1gen->pt()    , theWeight);
  theHistograms.fill("WGen_Pt"   , "WGen_Pt"   , 300, 0, 300, Wgen->pt()     , theWeight);
   
  //----------------------------

  theHistograms.fill<TH1I>("Number of leptons", "Number of leptons", 10, 0, 10, Genl.size());
  theHistograms.fill<TH1I>("Number of partons", "Number of partons", 10, 0, 10, Genj.size());



  //================================================================//
  //                          reco Particles                        //    
  //================================================================//

 
  //%%%%%%%% Comparison genParticles - recoParticles %%%%%%%%//

  std::vector< std::pair<const Particle*, const Particle* > > Vgen_reco;

  std::vector<const Particle* > Zll;

  foreach(const Boson<Lepton>& z, *Zmm) {
    Vgen_reco.push_back(make_pair(Z0gen, & z));
    Vgen_reco.push_back(make_pair(Z1gen, & z));
    Zll.push_back(& z);      
  }

  foreach(const Boson<Electron>& z, *Zee) {
    Vgen_reco.push_back(make_pair(Z0gen, & z));
    Vgen_reco.push_back(make_pair(Z1gen, & z));
    Zll.push_back(& z);
  }
  
  std::stable_sort(Vgen_reco.begin(), Vgen_reco.end(), deltaRComparator());
 
  const Particle* Z0 = Vgen_reco.at(0).second;
  const Particle* Z1 = Vgen_reco.at(1).second;


  //%%%%%%%%%%%% Definition of the signal %%%%%%%%%%%%//  
 
  std::stable_sort(Zll.begin() ,Zll.end() ,MassComparator(ZMASS));
  std::stable_sort(Wjj->begin(),Wjj->end(),MassComparator(WMASS));
  
  const Particle* myZ0  = Zll.at(0);
  const Particle* myZ1  = Zll.at(1);
  Boson<Jet> myW        = Wjj->at(0);
  
  TLorentzVector p_myZ0 = myZ0->p4();
  TLorentzVector p_myZ1 = myZ1->p4();
 
  TLorentzVector p_myj1 = myW.daughter(0).p4();
  TLorentzVector p_myj2 = myW.daughter(1).p4();
  
  
  cout <<  "--------------- MASSES COMPARISON: Z gen  ||  Z reco matched with gen  ||  Z reco ---------------" << endl;
  cout << "Z0gen= " << Z0gen->p4().M() << "\tZ0matched = " << Z0->p4().M() << "\tZ0reco = " << Green(p_myZ0.M()) <<endl;
  cout << "Z1gen= " << Z1gen->p4().M() << "\tZ1matched = " << Z1->p4().M() << "\tZ1reco = " << Green(p_myZ1.M()) <<endl;
  cout << "    " << endl;
  
  
  if ( (p_myZ0 == Z0->p4() && p_myZ1 == Z1->p4()) || (p_myZ0 == Z1->p4() && p_myZ1 == Z0->p4())) {
    theHistograms.fill<TH1I>("Efficiency of signal definition", "Efficiency of signal definition", 3, 0, 3, 1);
    
    
    /////-----------------Histograms-------------------
    
    //------------Mass--------------
    
    theHistograms.fill("Z0_Mass" , "Z0_Mass" , 200, 0, 200, p_myZ0.M()  , theWeight);
    theHistograms.fill("Z1_Mass" , "Z1_Mass" , 200, 0, 200, p_myZ1.M()  , theWeight);
    theHistograms.fill("Wjj_Mass", "Wjj_Mass", 200, 0, 200, myW.p4().M(), theWeight);  
    
    //------------Pt--------------
    
    theHistograms.fill("Z0_Pt"   , "Z0_Pt"   , 300, 0, 300, myZ0->pt()  , theWeight);
    theHistograms.fill("Z1_Pt"   , "Z1_Pt"   , 300, 0, 300, myZ1->pt()  , theWeight);
    theHistograms.fill("W_Pt"    , "W_Pt"    , 300, 0, 300, myW.pt()    , theWeight);  
    
    //------------Mass 6f--------------
    
    TLorentzVector p_6f = p_myZ0 + p_myZ1 + p_myj1 + p_myj2;
    
    theHistograms.fill("6f_Mass" , "6f_Mass" , 3000, 0, 3000, p_6f.M(), theWeight);
    
  }  
}


