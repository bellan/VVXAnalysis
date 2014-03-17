#include "VVXAnalysis/TreeAnalysis/interface/ZZWAnalyzer.h"
#include "VVXAnalysis/TreeAnalysis/interface/Colours.h"
                                                                                                                                                                                                                                                                                                                                                                           
#include <boost/foreach.hpp>
#define foreach BOOST_FOREACH

using namespace phys;
using namespace std;
using namespace colour;



Int_t ZZWAnalyzer::cut() {


  //================================================================//
  //                          reco Particles                        //    
  //================================================================//
    
 
  theHistograms.fill("Number of events", "Number of events", 10, 0, 10, 0, theWeight);  // 0: Number of total events --------------------------------------------------------------------------------------
  
  theHistograms.fill("MET", "MET", 300, 0, 300, met->p4().E(), theWeight);


  bool passZsize = ( Zmm->size() + Zee->size() ) >= 2;                //  1: Events with 2 well-defined Z bosons  -----------------------------------------------------------------------------------------
  
  if(passZsize)  theHistograms.fill("Number of events", 10, 0, 10, 1, theWeight);    
  
  
  bool passWsize =  Wjj->size() >= 1;
  
  bool pass = passZsize && passWsize;                     //  2: Events with 2 well-defined Z bosons and 1 well-defined W boson ---------------------------------------------------------------------------
  
  if(pass) {
    
    theHistograms.fill("Number of events", 10, 0, 10, 2, theWeight);  
       
    std::vector<Boson<Lepton> > Zll; 
    
    foreach(const Boson<Lepton>& z, *Zmm) {
      Zll.push_back(z.clone<Lepton>()); 
    }
  
    foreach(const Boson<Electron>& z, *Zee) { 
      Zll.push_back(z.clone<Lepton>());
    }
  
    std::stable_sort(Zll.begin() ,Zll.end() ,MassComparator(ZMASS));
    std::stable_sort(Wjj->begin(),Wjj->end(),MassComparator(WMASS));
    //std::stable_sort(Wjj->begin(),Wjj->end(),WPtComparator());
  
    myZ0 = Zll.at(0);
    myZ1 = Boson<Lepton>();
    myW  = Wjj->at(0);
  

    bool passGhost = true;                                 // 3: Events with 2 well-defined Z bosons, 1 well-defined W boson, no wrong leptons pairing------------------------------------------------------
  
    for(vector<Boson<Lepton> >::iterator b = Zll.begin()+1; b != Zll.end(); ++b) { 
    
      double DR00 = deltaR(myZ0.daughter(0).p4().Rapidity(), myZ0.daughter(0).p4().Phi(), b->daughter(0).p4().Rapidity(), b->daughter(0).p4().Phi());
      double DR01 = deltaR(myZ0.daughter(0).p4().Rapidity(), myZ0.daughter(0).p4().Phi(), b->daughter(1).p4().Rapidity(), b->daughter(1).p4().Phi());
      double DR10 = deltaR(myZ0.daughter(1).p4().Rapidity(), myZ0.daughter(1).p4().Phi(), b->daughter(0).p4().Rapidity(), b->daughter(0).p4().Phi());
      double DR11 = deltaR(myZ0.daughter(1).p4().Rapidity(), myZ0.daughter(1).p4().Phi(), b->daughter(1).p4().Rapidity(), b->daughter(1).p4().Phi());

      if (DR00 < 0.02 || DR01 < 0.02 || DR10 < 0.02 || DR11 < 0.02) {
	passGhost=false;
      }    
      if(passGhost) {  
	myZ1 = *b;
	break;
      }   
    }  

    if(passGhost == false) return -1;

  
    theHistograms.fill("Number of events", 10, 0, 10, 3, theWeight);   
  
  
    bool passllLowMass = true;                       // 4: Events with 2 well-defined Z bosons, 1 well-defined W boson, no wrong leptons pairing, massll > 4 GeV---------------------------------------------
    
    for(int i = 0; i <=1; ++i) {
      if ( myZ0.daughter(0).charge() != myZ1.daughter(i).charge() ) {
    	if ( (myZ0.daughter(0).p4() + myZ1.daughter(i).p4()).M() < 4 || (myZ0.daughter(1).p4() + myZ1.daughter((i+1)%2).p4()).M() < 4 ) passllLowMass = false;
      }
    }
    
    if (passllLowMass == false)  return -1;
    
    
    theHistograms.fill("Number of events", 10, 0, 10, 4, theWeight);   
    
    
    bool passLeptonsPt = false;                       // 5: Events with 2 well-defined Z bosons, 1 well-defined W boson, no wrong leptons pairing, massll > 4 GeV, 1lepton pt>10 and 1lepton pt >20------------
    
    if (myZ0.daughter(0).pt() > 10 || myZ0.daughter(1).pt() > 10 || myZ1.daughter(0).pt() > 10 || myZ1.daughter(1).pt() > 10 ) {
      if (myZ0.daughter(0).pt() > 20 || myZ0.daughter(1).pt() > 20 || myZ1.daughter(0).pt() > 20 || myZ1.daughter(1).pt() > 20 ) {
	passLeptonsPt = true;
      }
    }

    if(passLeptonsPt = false) return -1;



    theHistograms.fill("Number of events", 10, 0, 10, 5, theWeight);   


    bool passMET = true;                    // 6: Events with 2 well-defined Z bosons, 1 well-defined W boson, no wrong leptons pairing, massll > 4 GeV, 1lepton pt>10 and 1lepton pt >20, MET < 100--------------

    if( met->p4().E() > 100 ) passMET = false;

    if(passMET  = false) return -1;


    theHistograms.fill("Number of events", 10, 0, 10, 6, theWeight);  

    TLorentzVector p_myZ0 = myZ0.p4();
    TLorentzVector p_myZ1 = myZ1.p4();
    TLorentzVector p_myW  = myW.p4();
      
    TLorentzVector p_myj1 = myW.daughter(0).p4();
    TLorentzVector p_myj2 = myW.daughter(1).p4();

    cout << Red("\n============Event: ") << Red(event) << Red(" ==================")  << endl;
    cout << "-----------reco Particles Masses ------------" << endl;
    cout << "Z0reco = " << Green(myZ0.p4().M()) <<endl;
    cout << "Z1reco = " << Green(myZ1.p4().M()) <<endl;
    cout << "Wreco  = " << Green(myW.p4().M())  <<endl;
  

    /////-----------------Histograms-------------------
    
    //------------Mass-------------
    
    theHistograms.fill("Z0_Mass" , "Z0_Mass" , 200, 0, 200, p_myZ0.M() , theWeight);
    theHistograms.fill("Z1_Mass" , "Z1_Mass" , 200, 0, 200, p_myZ1.M() , theWeight);
    theHistograms.fill("Wjj_Mass", "Wjj_Mass", 200, 0, 200, p_myW.M()  , theWeight);  
      
    //------------Pt--------------
      
    theHistograms.fill("Z0_Pt"   , "Z0_Pt"   , 300, 0, 300, myZ0.pt() , theWeight);
    theHistograms.fill("Z1_Pt"   , "Z1_Pt"   , 300, 0, 300, myZ1.pt() , theWeight);
    theHistograms.fill("W_Pt"    , "W_Pt"    , 300, 0, 300, myW.pt()  , theWeight);  
      
    //------------Mass 6f----------
      
    TLorentzVector p_6f = p_myZ0 + p_myZ1 + p_myj1 + p_myj2;
      
    theHistograms.fill("6f_Mass" , "6f_Mass" , 3000, 0, 3000, p_6f.M(), theWeight);

    //-------------MET-------------

    theHistograms.fill("MET after cut", "MET after cut", 300, 0, 300, met->p4().M(), theWeight);

      
    theCutCounter += theWeight; 

  
    return pass ? 1 : -1;   

  }
    
  else return -1;

}


void ZZWAnalyzer::analyze() {

  if (genCategory == 999) genCategory = 10;
  
  cout << Yellow("\nCategory= ") << Yellow(genCategory) << endl;
  
  theHistograms.fill("Event category", "Event category", 11, 0, 11, genCategory, theWeight);  //// Event Category check //////
  
  
  //================================================================//
  //                             genParticles                       //    
  //================================================================//
  
  if (genVBParticles->size() >= 3) {
    
    std::vector<const Boson<Particle>* > GenZ;
    std::vector<const Boson<Particle>* > GenW;
    

    foreach(const Boson<Particle>& b, *genVBParticles) {
      
      int s_id = b.id();
      int id   = abs(s_id);
      
      if (id == 23)  GenZ.push_back(&b);      // Z
      else if (id == 24) GenW.push_back(&b);  // W
      
    }

    cout << "---------- genParticles: history information ----------" << endl; 
    cout << "Number of generated Z= "  << GenZ.size() << endl;
    cout << "Number of generated W= "  << GenW.size() << endl;
  
    if ( GenZ.size() >= 2 && GenW.size() >= 1 ) {
    
      const Boson<Particle>* Z0gen = GenZ.at(0);
      const Boson<Particle>* Z1gen = GenZ.at(1);
      const Boson<Particle>* Wgen  = GenW.at(0);
      
      
      /////-----------------Histograms-------------------
      
      //------------Mass--------------
  
      theHistograms.fill("Z0Gen_Mass", "Z0Gen_Mass", 200, 0, 200, Z0gen->p4().M(), theWeight);
      theHistograms.fill("Z1Gen_Mass", "Z1Gen_Mass", 200, 0, 200, Z1gen->p4().M(), theWeight);
      theHistograms.fill("WGen_Mass" , "WGen_Mass" , 200, 0, 200, Wgen->p4().M() , theWeight);
  
      //------------Pt--------------
  
      theHistograms.fill("Z0Gen_Pt"  , "Z0Gen_Pt"  , 300, 0, 300, Z0gen->pt()    , theWeight);
      theHistograms.fill("Z1Gen_Pt"  , "Z1Gen_Pt"  , 300, 0, 300, Z1gen->pt()    , theWeight);
      theHistograms.fill("WGen_Pt"   , "WGen_Pt"   , 300, 0, 300, Wgen->pt()     , theWeight);
   

 
      //%%%%%%%% Comparison genParticles - recoParticles %%%%%%%%//

      std::vector< std::pair<const Particle*, const Particle* > > ZcomparatorVector;
      std::vector< std::pair<const Particle*, const Particle* > > WcomparatorVector;  
  


      foreach(const Boson<Lepton>& z, *Zmm) {
	ZcomparatorVector.push_back(make_pair(Z0gen, & z));
	ZcomparatorVector.push_back(make_pair(Z1gen, & z));
      }
  
  
      foreach(const Boson<Electron>& z, *Zee) {   
	ZcomparatorVector.push_back(make_pair(Z0gen, & z));
	ZcomparatorVector.push_back(make_pair(Z1gen, & z));
      }
  
      foreach(const Boson<Jet>& w, *Wjj) {
	WcomparatorVector.push_back(make_pair(Wgen, & w));
      }
  
      std::stable_sort(ZcomparatorVector.begin(), ZcomparatorVector.end(), deltaRComparator());
      std::stable_sort(WcomparatorVector.begin(), WcomparatorVector.end(), deltaRComparator());

      const Particle* Z0 = ZcomparatorVector.at(0).second;         // Definition of correctly matched bosons
      const Particle* Z1 = ZcomparatorVector.at(1).second;         //
      const Particle* W  = WcomparatorVector.at(0).second;          //
  
      
      cout <<  "\n------- MASSES COMPARISON: Z gen  ||  Z reco matched with gen  ||  Z reco -------"   << endl;
      cout << "Z0gen= " << Z0gen->p4().M() << "\t\tZ0matched = " << Z0->p4().M() << "\t\tZ0reco = " << Green(myZ0.p4().M()) <<endl;
      cout << "Z1gen= " << Z1gen->p4().M() << "\t\tZ1matched = " << Z1->p4().M() << "\t\tZ1reco = " << Green(myZ1.p4().M()) <<endl;
      cout << "Wgen=  " << Wgen->p4().M()  << "\t\tWmatched =  " << W->p4().M()  << "\t\tWreco = "  << Green(myW.p4().M())  <<endl;

      
      bool ZcorrectMatch = (myZ0.p4() == Z0->p4() && myZ1.p4() == Z1->p4()) || (myZ0.p4() == Z1->p4() && myZ1.p4() == Z0->p4());
      bool WcorrectMatch = myW.p4() == W->p4();
      
      if ( ZcorrectMatch ) theHistograms.fill("Efficiency of Z definition", "Efficiency of Z definition", 3, 0, 3, 1, theWeight);
      
      if ( WcorrectMatch ) theHistograms.fill("Efficiency of W definition", "Efficiency of W definition", 3, 0, 3, 1, theWeight);
      
      if ( ZcorrectMatch && WcorrectMatch ) theHistograms.fill("Efficiency of signal definition", "Efficiency of signal definition", 3, 0, 3, 1, theWeight);
         
    }

    else {

      cout << "\nNOT A GENERATED SIGNAL EVENT" << endl;
      return;
    }

  }
  
  else {
    cout << "\nNOT A GENERATED SIGNAL EVENT" << endl;
    return;
  }
}

