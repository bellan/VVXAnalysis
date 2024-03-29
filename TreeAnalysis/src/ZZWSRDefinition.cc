#include "VVXAnalysis/TreeAnalysis/interface/ZZWSRDefinition.h"
#include "VVXAnalysis/Commons/interface/Colours.h"
#include "VVXAnalysis/Commons/interface/Comparators.h"

#include <boost/foreach.hpp>
#define foreach BOOST_FOREACH

using namespace phys;
using namespace std;
using namespace colour;



Int_t ZZWSRDefinition::cut() {



  //================================================================//
  //                          reco Particles                        //    
  //================================================================//
  
 
  theHistograms->fill("Number of events", "Number of events", 10, 0, 10, 0, theWeight);  // 0: Number of total events ----------------------------------------------------------------------------------------
  
  theHistograms->fill("MET_beforeCuts", "MET_beforeCuts", 300, 0, 300, met->pt(), theWeight);
  
  theHistograms->fill("Selection", "Selection", 3, 0, 3, 0, theWeight);


  foreach(const Lepton m, *muons) {
    theHistograms->fill("Leptons_P_{T}_beforeCuts", "Leptons_P_{T}_beforeCuts", 300, 0, 300, m.p4().Pt(), theWeight);
  }

  foreach(const Lepton e, *electrons) {
    theHistograms->fill("Leptons_P_{T}_beforeCuts", 300, 0, 300, e.p4().Pt(), theWeight );
  }

  foreach(const Jet j, *jets) {
    theHistograms->fill("Jets_P_{T}_beforeCuts", "Jets_P_{T}_beforeCuts", 300, 0, 300, j.p4().Pt(), theWeight );
  }
  
  //.......W request

  bool passWsize =  Vhad->size() >= 1;    
  
  if (passWsize) theHistograms->fill("Selection", "Selection", 3, 0, 3, 2, theWeight);
  
  
  //.......Z0,Z1 definition

  bool passZsize = 0; // FIXME
                   
  bool passGhost = true;   
  
  bool passllLowMass = true;
  
  if (passZsize) {
    
    std::vector<Boson<Lepton> > Zll; 
    
    //    foreach(const Boson<Lepton>& z, *Z) {
    //  Zll.push_back(z.clone<Lepton>()); 
    // }
    
    
    std::stable_sort(Zll.begin() ,Zll.end() ,MassComparator(ZMASS));
    
    myZ0 = Zll.at(0);  // Fixme, it will crash
    myZ1 = Boson<Lepton>();
        
    for(vector<Boson<Lepton> >::iterator b = Zll.begin()+1; b != Zll.end(); ++b) { 
    
      double DR00 = physmath::deltaR(myZ0.daughter(0), b->daughter(0));
      double DR01 = physmath::deltaR(myZ0.daughter(0), b->daughter(1));
      double DR10 = physmath::deltaR(myZ0.daughter(1), b->daughter(0));
      double DR11 = physmath::deltaR(myZ0.daughter(1), b->daughter(1));

      if (DR00 < 0.02 || DR01 < 0.02 || DR10 < 0.02 || DR11 < 0.02) {
	passGhost=false; 
	// BUG! If the first Z1 canidate fails, then the remaining candidate can never have passGhost == true  ?
	cout << "Check me! Most probably I am a bug!" << endl;
      } 
      
      if(passGhost) {  
	myZ1 = *b;
	break;
      }   
    } 
    
    if(passGhost == false) return -1;

    for(int i = 0; i <=1; ++i) {
      if ( myZ0.daughter(0).charge() != myZ1.daughter(i).charge() ) {
	if ( (myZ0.daughter(0).p4() + myZ1.daughter(i).p4()).M() < 4 || (myZ0.daughter(1).p4() + myZ1.daughter((i+1)%2).p4()).M() < 4 ) passllLowMass = false;
      }
    }
  
    if (passllLowMass == false)  return -1;
  
  }
  
  bool Zdef = (passZsize && passGhost && passllLowMass);
  if (Zdef) {
    theHistograms->fill("Number of events", 10, 0, 10, 1, theWeight);     //  1: Events with at lest 2 well-defined Z bosons  -----------------------------------------------------------------------------------
    theHistograms->fill("Selection", "Selection", 3, 0, 3, 1, theWeight);
  }  
 
  
 
  bool pass = Zdef && passWsize;                    
  
  if(pass) {
       
    std::stable_sort(Vhad->begin(),Vhad->end(),MassComparator(WMASS)); // FIXME: Z too!
    //std::stable_sort(Wjj->begin(),Wjj->end(),WJetPtComparator());
  
    myW  = Vhad->at(0);
    //  bool WmassRange = (fabs(myW.p4().M() - WMASS) < 40);
    
    cout << endl << Green("---------------Selected W---------------") << endl;
    cout << "Pt Jet1 myW = " << myW.daughter(0).pt() << endl;
    cout << "Pt Jet2 myW = " << myW.daughter(1).pt() << endl;
    cout << "----------------------------------------" << endl;
    
    theHistograms->fill("Number of events", 10, 0, 10, 2, theWeight);      //  2: Events with at lest 2 well-defined Z bosons and 1 well-defined W boson  -----------------------------------------------------------------------------------
   
    //   if(!WmassRange) return -1;  

    
                   
    //..........................myZ0, myZ1, myW defined : BEGINNING OF CUTS  ..........................//

 
    bool passLeptonsPt = false;                       // 3: Events with 2 well-defined Z bosons, 1 well-defined W boson, 1lepton pt>10 and 1lepton pt >20------------
    
    int count10 =0;
    int count20 = 0;
    
    for (int j = 0; j<=1; ++j) {
      if (myZ0.daughter(j).pt() > 10 ) ++count10;
      if (myZ1.daughter(j).pt() > 10 ) ++count10;	
    }
    
    
    if (myZ0.daughter(0).pt() > 20 || myZ0.daughter(1).pt() > 20 || myZ1.daughter(0).pt() > 20 || myZ1.daughter(1).pt() > 20 ) {
      ++count20;
    }
    
    if ( count10 >= 2 && count20 == 1 ) passLeptonsPt = true;
    
    if(passLeptonsPt == false) return -1;
    
   
    theHistograms->fill("Number of events", 10, 0, 10, 3, theWeight);   


    if ( met->pt() > 80 ) return -1;        // 4: Events with 2 well-defined Z bosons, 1 well-defined W boson, 1lepton pt>10 and 1lepton pt >20, MET < 80--------------

    theHistograms->fill("Number of events", 10, 0, 10, 4, theWeight);  


    TLorentzVector p_myZ0 = myZ0.p4();
    TLorentzVector p_myZ1 = myZ1.p4();
    TLorentzVector p_myW  = myW.p4();


    TLorentzVector p_myl1 = myZ0.daughter(0).p4();
    TLorentzVector p_myl2 = myZ0.daughter(1).p4();
    TLorentzVector p_myl3 = myZ1.daughter(0).p4();
    TLorentzVector p_myl4 = myZ1.daughter(1).p4();
    TLorentzVector p_myj1 = myW.daughter(0).p4();
    TLorentzVector p_myj2 = myW.daughter(1).p4();

    TLorentzVector p_4l = p_myZ0 + p_myZ1; 
    TLorentzVector p_6f = p_4l  + p_myj1 + p_myj2;


    if ( p_6f.M() < 300 ) return -1;         // 5: Events with 2 well-defined Z bosons, 1 well-defined W boson, no wrong leptons pairing, massll > 4 GeV, 1lepton pt>10 and 1lepton pt >20, MET < 80, m6f>300GeV--------------

    theHistograms->fill("Number of events", 10, 0, 10, 5, theWeight);  


    double ptj1 = max(p_myj1.Pt(), p_myj2.Pt());
    double ptj2 = min(p_myj1.Pt(), p_myj2.Pt());

    double deltaRJets11 = physmath::deltaR(phys::Particle(p_myj1), phys::Particle(p_myl1));
    double deltaRJets12 = physmath::deltaR(phys::Particle(p_myj1), phys::Particle(p_myl2));
    double deltaRJets13 = physmath::deltaR(phys::Particle(p_myj1), phys::Particle(p_myl3));
    double deltaRJets14 = physmath::deltaR(phys::Particle(p_myj1), phys::Particle(p_myl4));
    double deltaRJets21 = physmath::deltaR(phys::Particle(p_myj2), phys::Particle(p_myl1));
    double deltaRJets22 = physmath::deltaR(phys::Particle(p_myj2), phys::Particle(p_myl2));
    double deltaRJets23 = physmath::deltaR(phys::Particle(p_myj2), phys::Particle(p_myl3));
    double deltaRJets24 = physmath::deltaR(phys::Particle(p_myj2), phys::Particle(p_myl4));
    
    double DPhi_j1_j2   = physmath::deltaPhi(p_myj1, p_myj2);
    double DPhi_ZZ_j1j2 = physmath::deltaPhi((p_myj1 + p_myj2), (p_myZ0 + p_myZ1));
    double DPhi_Z_Z     = physmath::deltaPhi(p_myZ0, p_myZ1); 
    double DPhi_ZZ_W    = physmath::deltaPhi((p_myZ0 + p_myZ1), p_myW);
    double DPhi_Z0_W    = physmath::deltaPhi(p_myZ0, p_myW);
    double DPhi_Z1_W    = physmath::deltaPhi(p_myZ1, p_myW);


    cout << endl << Red("============Selected Event: ") << Red(event) << Red(" w= ") << Red(theWeight) << " " << Red(" ==================")  << endl;
    cout << "-----------reco Particles Masses ------------" << endl;
    cout << "Z0reco = " << Green(myZ0.p4().M()) <<endl;
    cout << "Z1reco = " << Green(myZ1.p4().M()) <<endl;
    cout << "Wreco  = " << Green(myW.p4().M())  <<endl;
    cout << "Lepton 1" << "\tPt = "   << p_myl1.Pt() << "\tEta = " <<  p_myl1.Eta() << "\tPhi = " << p_myl1.Phi() << endl; 
    cout << "Lepton 2" << "\tPt = "   << p_myl2.Pt() << "\tEta = " <<  p_myl2.Eta() << "\tPhi = " << p_myl2.Phi() << endl; 
    cout << "Lepton 3" << "\tPt = "   << p_myl3.Pt() << "\tEta = " <<  p_myl3.Eta() << "\tPhi = " << p_myl3.Phi() << endl; 
    cout << "Lepton 4" << "\tPt = "   << p_myl4.Pt() << "\tEta = " <<  p_myl4.Eta() << "\tPhi = " << p_myl4.Phi() << endl;   
    cout << "Jet 1"    << "\t\tPt = " << p_myj1.Pt() << "\tEta = " <<  p_myj1.Eta() << "\t\tPhi = " << p_myj1.Phi() << "\tdeltaR with leptons 1,2,3,4= " << deltaRJets11 << ", " << deltaRJets12 << ", " << deltaRJets13 << ", " << deltaRJets14 << endl;
    cout << "Jet 2"    << "\t\tPt = " << p_myj2.Pt() << "\tEta = " <<  p_myj2.Eta() << "\t\tPhi = " << p_myj2.Phi() << "\tdeltaR with leptons 1,2,3,4= " << deltaRJets21 << ", " << deltaRJets22 << ", " << deltaRJets23 << ", " << deltaRJets24 << endl;
     

    /////-----------------Histograms-------------------
  
    if(genParticles->size() == 9)  theHistograms->fill("0 jet", "0 jet", 3, 0, 3, 1, theWeight);
    
    if(genParticles->size() == 10) theHistograms->fill("1 jet", "1 jet", 3, 0, 3, 1, theWeight);
    
    
    //------------Mass-------------
    
    theHistograms->fill("Z0_Mass" , "Z0_Mass" , 200, 0, 200, p_myZ0.M() , theWeight);
    theHistograms->fill("Z1_Mass" , "Z1_Mass" , 200, 0, 200, p_myZ1.M() , theWeight);
    theHistograms->fill("Wjj_Mass", "Wjj_Mass", 200, 0, 200, p_myW.M()  , theWeight);  
    
      
    //------------Pt--------------
      
    theHistograms->fill("Z0_Pt"   , "Z0_Pt"   , 300, 0, 300, myZ0.pt() , theWeight);
    theHistograms->fill("Z1_Pt"   , "Z1_Pt"   , 300, 0, 300, myZ1.pt() , theWeight);
    theHistograms->fill("W_Pt"    , "W_Pt"    , 300, 0, 300, myW.pt()  , theWeight);  

    theHistograms->fill("j1_Pt"   , "j1_Pt"   , 300, 0, 300, ptj1      , theWeight);
    theHistograms->fill("j2_Pt"   , "j2_Pt"   , 300, 0, 300, ptj2      , theWeight);

    theHistograms->fill("TotalPt", "TotalPt", 100, 0, 50, (p_myZ0 + p_myZ1 + p_myW).Pt(), theWeight);

      
    //----------Fermions Masses--------
          
    theHistograms->fill("6f_Mass" , "6f_Mass" , 1500, 0, 1500, p_6f.M(), theWeight);

    theHistograms->fill("4l_Mass" , "4l_Mass" , 1500, 0 ,1500, p_4l.M() , theWeight);

    //-------------MET-------------

    theHistograms->fill("MET after cut", "MET after cut", 300, 0, 300, met->pt(), theWeight);

    //-----------delta_R-----------

    theHistograms->fill("Reco_DR_J1L1", "Reco_DR_J1L1", 100, 0, 8, deltaRJets11, theWeight);
    theHistograms->fill("Reco_DR_J1L2", "Reco_DR_J1L2", 100, 0, 8, deltaRJets12, theWeight);
    theHistograms->fill("Reco_DR_J1L3", "Reco_DR_J1L3", 100, 0, 8, deltaRJets13, theWeight);
    theHistograms->fill("Reco_DR_J1L4", "Reco_DR_J1L4", 100, 0, 8, deltaRJets14, theWeight);
    theHistograms->fill("Reco_DR_J2L1", "Reco_DR_J2L1", 100, 0, 8, deltaRJets21, theWeight);
    theHistograms->fill("Reco_DR_J2L2", "Reco_DR_J2L2", 100, 0, 8, deltaRJets22, theWeight);
    theHistograms->fill("Reco_DR_J2L3", "Reco_DR_J2L3", 100, 0, 8, deltaRJets23, theWeight);
    theHistograms->fill("Reco_DR_J2L4", "Reco_DR_J2L4", 100, 0, 8, deltaRJets24, theWeight);

    //------------Kinematics-------------

    theHistograms->fill("DeltaPhi_j1_j2"   , "DeltaPhi_j1_j2"   , 50, 0, 5, DPhi_j1_j2	, theWeight);
    theHistograms->fill("DeltaPhi_Z_Z"     , "DeltaPhi_Z_Z"     , 50, 0, 5, DPhi_Z_Z  	, theWeight);
    theHistograms->fill("DeltaPhi_ZZ_W"    , "DeltaPhi_ZZ_W"    , 50, 0, 5, DPhi_ZZ_W 	, theWeight);
    theHistograms->fill("DeltaPhi_Z0_W"    , "DeltaPhi_Z0_W"    , 50, 0, 5, DPhi_Z0_W 	, theWeight);
    theHistograms->fill("DeltaPhi_Z1_W"    , "DeltaPhi_Z1_W"    , 50, 0, 5, DPhi_Z1_W 	, theWeight);
    theHistograms->fill("DeltaPhi_ZZ_j1j2" , "DeltaPhi_ZZ_j1j2" , 50, 0, 5, DPhi_ZZ_j1j2 , theWeight);
      

    theCutCounter += theWeight; 

  
    return pass ? 1 : -1;   

  }
    
  else return -1;

}

 


void ZZWSRDefinition::analyze() {

  if (genCategory == 999) genCategory = 10;
  
  cout << Yellow("\nCategory= ") << Yellow(genCategory) << endl;
  
  theHistograms->fill("Event category", "Event category", 11, 0, 11, genCategory, theWeight);  //// Event Category check //////
  
 

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
    
      const Boson<Particle>& Z0gen = *GenZ.at(0);
      const Boson<Particle>& Z1gen = *GenZ.at(1);
      const Boson<Particle>& Wgen  = *GenW.at(0);
      
      double DR11 = physmath::deltaR(Wgen.daughter(0), Z0gen.daughter(0));
      double DR12 = physmath::deltaR(Wgen.daughter(0), Z0gen.daughter(1));
      double DR13 = physmath::deltaR(Wgen.daughter(0), Z1gen.daughter(0));
      double DR14 = physmath::deltaR(Wgen.daughter(0), Z1gen.daughter(1));
      double DR21 = physmath::deltaR(Wgen.daughter(1), Z0gen.daughter(0));
      double DR22 = physmath::deltaR(Wgen.daughter(1), Z0gen.daughter(1));
      double DR23 = physmath::deltaR(Wgen.daughter(1), Z1gen.daughter(0));
      double DR24 = physmath::deltaR(Wgen.daughter(1), Z1gen.daughter(1));
      
      
      
      /////-----------------Histograms-------------------
      
      //------------Mass--------------
  
      theHistograms->fill("Z0Gen_Mass", "Z0Gen_Mass", 200, 0, 200, Z0gen.p4().M(), theWeight);
      theHistograms->fill("Z1Gen_Mass", "Z1Gen_Mass", 200, 0, 200, Z1gen.p4().M(), theWeight);
      theHistograms->fill("WGen_Mass" , "WGen_Mass" , 200, 0, 200, Wgen.p4().M() , theWeight);
  
      //------------Pt--------------
  
      theHistograms->fill("Z0Gen_Pt"  , "Z0Gen_Pt"  , 300, 0, 300, Z0gen.pt()    , theWeight);
      theHistograms->fill("Z1Gen_Pt"  , "Z1Gen_Pt"  , 300, 0, 300, Z1gen.pt()    , theWeight);
      theHistograms->fill("WGen_Pt"   , "WGen_Pt"   , 300, 0, 300, Wgen.pt()     , theWeight);

      theHistograms->fill("j1Gen_Pt"  , "j1Gen_Pt"  , 300, 0, 300, Wgen.daughter(0).pt() , theWeight);
      theHistograms->fill("j2Gen_Pt"  , "j2Gen_Pt"  , 300, 0, 300, Wgen.daughter(1).pt() , theWeight);
   
      //------------DR--------------

      theHistograms->fill("Gen_DR_J1L1", "Gen_DR_J1L1", 100, 0, 8, DR11, theWeight);
      theHistograms->fill("Gen_DR_J1L2", "Gen_DR_J1L2", 100, 0, 8, DR12, theWeight);
      theHistograms->fill("Gen_DR_J1L3", "Gen_DR_J1L3", 100, 0, 8, DR13, theWeight);
      theHistograms->fill("Gen_DR_J1L4", "Gen_DR_J1L4", 100, 0, 8, DR14, theWeight);
      theHistograms->fill("Gen_DR_J2L1", "Gen_DR_J2L1", 100, 0, 8, DR21, theWeight);
      theHistograms->fill("Gen_DR_J2L2", "Gen_DR_J2L2", 100, 0, 8, DR22, theWeight);
      theHistograms->fill("Gen_DR_J2L3", "Gen_DR_J2L3", 100, 0, 8, DR23, theWeight);
      theHistograms->fill("Gen_DR_J2L4", "Gen_DR_J2L4", 100, 0, 8, DR24, theWeight);
      
      
      //%%%%%%%% Comparison genParticles - recoParticles %%%%%%%%//
      
      std::vector< std::pair<Boson<Particle> , Boson<Lepton> > > ZcomparatorVector;
      std::vector< std::pair<Boson<Particle> , Boson<Jet> > >    WcomparatorVector;  
  


      // foreach(const Boson<Lepton>& z, *Z) {
      // 	ZcomparatorVector.push_back(make_pair(Z0gen, z.clone<Lepton>()));
      // 	ZcomparatorVector.push_back(make_pair(Z1gen, z.clone<Lepton>()));
      // }
      
      foreach(const Boson<Jet> w, *Vhad) {  // FIXME: Z too!
	WcomparatorVector.push_back(make_pair(Wgen, w));
      }
  
      //std::stable_sort(ZcomparatorVector.begin(), ZcomparatorVector.end(), VdeltaRComparator());
      std::stable_sort(WcomparatorVector.begin(), WcomparatorVector.end(), VdeltaRComparator());

    
      
      //Boson<Lepton> Z0 = ZcomparatorVector.at(0).second;         // Definition of correctly matched bosons
      //Boson<Lepton> Z1 = ZcomparatorVector.at(1).second;         //
      Boson<Jet> W     = WcomparatorVector.at(0).second;         //
  
      
      // cout <<  "\n---------- MASSES COMPARISON: Z gen  ||  Z reco matched with gen  ||  Z reco ----------"   << endl;

      // cout << "Z0gen= " << Z0gen.p4().M() << "\t\tZ0matched = " << Z0.p4().M() << "\t\tZ0reco = " << Green(myZ0.p4().M()) <<endl;
      // cout << "Z1gen= " << Z1gen.p4().M() << "\t\tZ1matched = " << Z1.p4().M() << "\t\tZ1reco = " << Green(myZ1.p4().M()) <<endl;
      // cout << "Wgen=  " << Wgen.p4().M()  << "\t\tWmatched =  " << W.p4().M()  << "\t\tWreco = "  << Green(myW.p4().M())  <<endl;

      // cout <<  "\n------- DAUGHTERS IDS COMPARISON: Z gen  ||  Z reco matched with gen  ||  Z reco -------"   << endl;
      // cout << "daughter1 \tZ0gen = " << Z0gen.daughter(0).id() << "\tZ0matched = " << Z0.daughter(0).id() << "\tZ0reco = " << myZ0.daughter(0).id() << endl;
      // cout << "daughter2 \tZ0gen = " << Z0gen.daughter(1).id() << "\tZ0matched = " << Z0.daughter(1).id() << "\tZ0reco = " << myZ0.daughter(1).id() << endl;
      // cout << "daughter1 \tZ0gen = " << Z1gen.daughter(0).id() << "\tZ0matched = " << Z1.daughter(0).id() << "\tZ0reco = " << myZ1.daughter(0).id() << endl;
      // cout << "daughter2 \tZ0gen = " << Z1gen.daughter(1).id() << "\tZ0matched = " << Z1.daughter(1).id() << "\tZ0reco = " << myZ1.daughter(1).id() << endl;
      
      
      // bool ZcorrectMatch = (myZ0.p4() == Z0.p4() && myZ1.p4() == Z1.p4()) || (myZ0.p4() == Z1.p4() && myZ1.p4() == Z0.p4());
      // bool WcorrectMatch = myW.p4() == W.p4();
      
      // if ( ZcorrectMatch ) theHistograms->fill("Efficiency of Z definition", "Efficiency of Z definition", 3, 0, 3, 1, theWeight);
      
      // if ( WcorrectMatch ) theHistograms->fill("Efficiency of W definition", "Efficiency of W definition", 3, 0, 3, 1, theWeight);
      
      // if ( ZcorrectMatch && WcorrectMatch ) theHistograms->fill("Efficiency of signal definition", "Efficiency of signal definition", 3, 0, 3, 1, theWeight);

      // double deltaRJet1Lep1 = deltaR(W.daughter(0).p4().Rapidity(), W.daughter(0).p4().Phi(), Z0.daughter(0).p4().Rapidity(), Z0.daughter(0).p4().Phi());
      // double deltaRJet1Lep2 = deltaR(W.daughter(0).p4().Rapidity(), W.daughter(0).p4().Phi(), Z0.daughter(1).p4().Rapidity(), Z0.daughter(1).p4().Phi());
      // double deltaRJet1Lep3 = deltaR(W.daughter(0).p4().Rapidity(), W.daughter(0).p4().Phi(), Z1.daughter(0).p4().Rapidity(), Z1.daughter(0).p4().Phi());
      // double deltaRJet1Lep4 = deltaR(W.daughter(0).p4().Rapidity(), W.daughter(0).p4().Phi(), Z1.daughter(1).p4().Rapidity(), Z1.daughter(1).p4().Phi());
      // double deltaRJet2Lep1 = deltaR(W.daughter(1).p4().Rapidity(), W.daughter(1).p4().Phi(), Z0.daughter(0).p4().Rapidity(), Z0.daughter(0).p4().Phi());
      // double deltaRJet2Lep2 = deltaR(W.daughter(1).p4().Rapidity(), W.daughter(1).p4().Phi(), Z0.daughter(1).p4().Rapidity(), Z0.daughter(1).p4().Phi());
      // double deltaRJet2Lep3 = deltaR(W.daughter(1).p4().Rapidity(), W.daughter(1).p4().Phi(), Z1.daughter(0).p4().Rapidity(), Z1.daughter(0).p4().Phi());
      // double deltaRJet2Lep4 = deltaR(W.daughter(1).p4().Rapidity(), W.daughter(1).p4().Phi(), Z1.daughter(1).p4().Rapidity(), Z1.daughter(1).p4().Phi());

      // theHistograms->fill("Matched_DR_J1L1", "Matched_DR_J1L1", 100, 0, 8, deltaRJet1Lep1, theWeight);
      // theHistograms->fill("Matched_DR_J1L2", "Matched_DR_J1L2", 100, 0, 8, deltaRJet1Lep2, theWeight);
      // theHistograms->fill("Matched_DR_J1L3", "Matched_DR_J1L3", 100, 0, 8, deltaRJet1Lep3, theWeight);
      // theHistograms->fill("Matched_DR_J1L4", "Matched_DR_J1L4", 100, 0, 8, deltaRJet1Lep4, theWeight);
      // theHistograms->fill("Matched_DR_J2L1", "Matched_DR_J2L1", 100, 0, 8, deltaRJet2Lep1, theWeight);
      // theHistograms->fill("Matched_DR_J2L2", "Matched_DR_J2L2", 100, 0, 8, deltaRJet2Lep2, theWeight);
      // theHistograms->fill("Matched_DR_J2L3", "Matched_DR_J2L3", 100, 0, 8, deltaRJet2Lep3, theWeight);
      // theHistograms->fill("Matched_DR_J2L4", "Matched_DR_J2L4", 100, 0, 8, deltaRJet2Lep4, theWeight);
      
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



