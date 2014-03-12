#include "VVXAnalysis/TreeAnalysis/interface/ZZWAnalyzer.h"
#include "VVXAnalysis/TreeAnalysis/interface/Colours.h"
                                                                                                                                                                                                                                                                                                                                                                           
#include <boost/foreach.hpp>
#define foreach BOOST_FOREACH

using namespace phys;
using namespace std;
using namespace colour;



Int_t ZZWAnalyzer::cut() {


  bool background = genCategory != 0;

  if (background) theHistograms.fill("Number of background events", "Number of background events", 10, 0, 10, 0, theWeight);  // Number of total background events -----------------------------------------


  bool signal = genCategory == 0;
  
  if(signal)
    theHistograms.fill("Number of signal events", "Number of signal events", 10, 0, 10, 0, theWeight);           // Number of total signal events ----------------------------------------------------------


  bool passZsize = ( Zmm->size() + Zee->size() ) >= 2;
 
  if(signal && passZsize)     theHistograms.fill("Number of signal events"    , 10, 0, 10, 1, theWeight);               //Number of events with 2 well-defined Z bosons  ------------------------------------                              
  if(background && passZsize) theHistograms.fill("Number of background events", 10, 0, 10, 1, theWeight);


  bool passWsize =  Wjj->size() >= 1;

 
  if(background && passZsize && passWsize) theHistograms.fill("Number of background events", 10, 0, 10, 2, theWeight);  //Number of events with 2 well-defined Z bosons and 1 well-defined W boson -----------
  
  bool pass = signal && passZsize && passWsize;
  
  if(pass) {
    ++theCutCounter;
    theHistograms.fill("Number of signal events" , 10, 0, 10, 2, theWeight);              //Number of events with 2 well-defined Z bosons and 1 well-defined W boson -----------------------------------------
  }

  return pass ? 1 : -1;
      
}


void ZZWAnalyzer::analyze() {

  
  //================================================================//
  //                             genParticles                       //    
  //================================================================//
  

  std::vector<const Boson<Particle>* > GenZ;
  std::vector<const Boson<Particle>* > GenW;
  
  foreach(const Boson<Particle>& b, *genVBParticles) {
    
    int s_id = b.id();
    int id   = abs(s_id);

    if (id == 23)  GenZ.push_back(&b);      // Z
    else if (id == 24) GenW.push_back(&b);  // W
 
  }

  const Boson<Particle>* Z0gen = GenZ.at(0);
  const Boson<Particle>* Z1gen = GenZ.at(1);
  const Boson<Particle>* Wgen  = GenW.at(0);


  cout << "\n=============== genParticles: history information ===============" << endl; 
  cout << "Number of generated Z= "  << GenZ.size() << endl;
  cout << "Number of generated W= "  << GenW.size() << endl;
  cout << "Category= "               << genCategory << endl;


  /////-----------------Histograms-------------------
  
  //------------Mass--------------
  
  theHistograms.fill("Z0Gen_Mass", "Z0Gen_Mass", 200, 0, 200, Z0gen->p4().M(), theWeight);
  theHistograms.fill("Z1Gen_Mass", "Z1Gen_Mass", 200, 0, 200, Z1gen->p4().M(), theWeight);
  theHistograms.fill("WGen_Mass" , "WGen_Mass" , 200, 0, 200, Wgen->p4().M() , theWeight);
  
  //------------Pt--------------
  
  theHistograms.fill("Z0Gen_Pt"  , "Z0Gen_Pt"  , 300, 0, 300, Z0gen->pt()    , theWeight);
  theHistograms.fill("Z1Gen_Pt"  , "Z1Gen_Pt"  , 300, 0, 300, Z1gen->pt()    , theWeight);
  theHistograms.fill("WGen_Pt"   , "WGen_Pt"   , 300, 0, 300, Wgen->pt()     , theWeight);
   


  //================================================================//
  //                          reco Particles                        //    
  //================================================================//

 
  //%%%%%%%% Comparison genParticles - recoParticles %%%%%%%%//

  std::vector< std::pair<const Particle*, const Particle* > > ZcomparatorVector;
  std::vector< std::pair<const Particle*, const Particle* > > WcomparatorVector;
  std::vector<const Particle* > Zll;
  
  Boson<Lepton>& Zmm1 = Zmm->at(0);
  Boson<Lepton>& Zmm2 = Zmm->at(0);
  
  
  foreach(const Boson<Lepton>& z, *Zmm) {
    ZcomparatorVector.push_back(make_pair(Z0gen, & z));
    ZcomparatorVector.push_back(make_pair(Z1gen, & z));
    Zll.push_back(& z);
  }
  
  
  foreach(const Boson<Electron>& z, *Zee) {   
    ZcomparatorVector.push_back(make_pair(Z0gen, & z));
    ZcomparatorVector.push_back(make_pair(Z1gen, & z));
    Zll.push_back(& z);
  }
  
  foreach(const Boson<Jet>& w, *Wjj) {
    WcomparatorVector.push_back(make_pair(Wgen, & w));
  }
  

  std::stable_sort(ZcomparatorVector.begin(), ZcomparatorVector.end(), deltaRComparator());
  std::stable_sort(WcomparatorVector.begin(), WcomparatorVector.end(), deltaRComparator());
  
  const Particle* Z0 = ZcomparatorVector.at(0).second;         // Definition of correctly matched bosons
  const Particle* Z1 = ZcomparatorVector.at(1).second;         //
  const Particle* W  = WcomparatorVector.at(0).first;    
  



  //%%%%%%%%%%%% Definition of the signal %%%%%%%%%%%%//  
 
  std::stable_sort(Zll.begin() ,Zll.end() ,MassComparator(ZMASS));
  std::stable_sort(Wjj->begin(),Wjj->end(),MassComparator(WMASS));
  //std::stable_sort(Wjj->begin(),Wjj->end(),WPtComparator());
  
  const Particle* myZ0  = Zll.at(0);


//   bool passGhost = true;
//     for (int i=0; i<=1; ++i){
//       for (int j=0; j<=1; ++j) {
// 	double DR = deltaR(Zmm1.daughter(i).p4().Rapidity(), Zmm1.daughter(i).p4().Phi(), z.daughter(j).p4().Rapidity(), z.daughter(j).p4().Phi());
// 	if (DR < 0.02) {
// 	  passGhost=false;
// 	  break;
// 	}
//       }
//     }
    
//     if (passGhost)  { 
//       Z2=z;
//       break;
//     }


//     if (passGhost)  {    
//       
//     } else Zllwrong.push_back(& z);
//  }
  
  // Se tutte le Z falliscono il taglio, allora Zmm1==Zmm2
  // Altrimenti la rtua seconda Z e' Zmm2



  const Particle* myZ1  = Zll.at(1);
  Boson<Jet> myW        = Wjj->at(0);
  
  TLorentzVector p_myZ0 = myZ0->p4();
  TLorentzVector p_myZ1 = myZ1->p4();
  TLorentzVector p_myW  = myW.p4();
 
  TLorentzVector p_myj1 = myW.daughter(0).p4();
  TLorentzVector p_myj2 = myW.daughter(1).p4();
  
  
  cout <<  "--------------- MASSES COMPARISON: Z gen  ||  Z reco matched with gen  ||  Z reco ---------------"   << endl;
  cout << "Z0gen= " << Z0gen->p4().M() << "\tZ0matched = " << Z0->p4().M() << "\tZ0reco = " << Green(p_myZ0.M()) <<endl;
  cout << "Z1gen= " << Z1gen->p4().M() << "\tZ1matched = " << Z1->p4().M() << "\tZ1reco = " << Green(p_myZ1.M()) <<endl;
  cout << "Wgen= "  << Wgen->p4().M()  << "\tWmatched = "  << W->p4().M()  << "\tWreco = "  << Green(p_myW.M())  <<endl;
  cout << "WcomparatorVector size = "  << WcomparatorVector.size() << endl;
  
  bool ZcorrectMatch = (p_myZ0 == Z0->p4() && p_myZ1 == Z1->p4()) || (p_myZ0 == Z1->p4() && p_myZ1 == Z0->p4());
  bool WcorrectMatch = p_myW == W->p4();

  if ( ZcorrectMatch ) theHistograms.fill("Efficiency of Z definition", "Efficiency of Z definition", 3, 0, 3, 1, theWeight);
  
  if ( WcorrectMatch ) theHistograms.fill("Efficiency of W definition", "Efficiency of W definition", 3, 0, 3, 1, theWeight);

  if ( ZcorrectMatch && WcorrectMatch ) theHistograms.fill<TH1I>("Efficiency of signal definition", "Efficiency of signal definition", 3, 0, 3, 1);

  
  /////-----------------Histograms-------------------
  
  //------------Mass-------------
  
  theHistograms.fill("Z0_Mass" , "Z0_Mass" , 200, 0, 200, p_myZ0.M() , theWeight);
  theHistograms.fill("Z1_Mass" , "Z1_Mass" , 200, 0, 200, p_myZ1.M() , theWeight);
  theHistograms.fill("Wjj_Mass", "Wjj_Mass", 200, 0, 200, p_myW.M()  , theWeight);  
  
  //------------Pt--------------
  
  theHistograms.fill("Z0_Pt"   , "Z0_Pt"   , 300, 0, 300, myZ0->pt() , theWeight);
  theHistograms.fill("Z1_Pt"   , "Z1_Pt"   , 300, 0, 300, myZ1->pt() , theWeight);
  theHistograms.fill("W_Pt"    , "W_Pt"    , 300, 0, 300, myW.pt()   , theWeight);  
  
  //------------Mass 6f----------
  
  TLorentzVector p_6f = p_myZ0 + p_myZ1 + p_myj1 + p_myj2;
  
  theHistograms.fill("6f_Mass" , "6f_Mass" , 3000, 0, 3000, p_6f.M(), theWeight);
   
}


