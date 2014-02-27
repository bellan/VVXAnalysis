/**
 *
 *  Simple analyzer of the generator history.
 *
 */


#include <FWCore/Framework/interface/Frameworkfwd.h>
#include <FWCore/Framework/interface/EDAnalyzer.h>
#include <FWCore/Framework/interface/ESHandle.h>
#include <FWCore/Framework/interface/Event.h>
#include <FWCore/ParameterSet/interface/ParameterSet.h>
#include <FWCore/ServiceRegistry/interface/Service.h>
#include <CommonTools/UtilAlgos/interface/TFileService.h>
#include <DataFormats/HepMCCandidate/interface/GenParticle.h>
#include <Calibration/IsolatedParticles/plugins/IsolatedGenParticles.h>

#include <TH1F.h>
#include <TH2F.h>

#include <iostream>
#include <vector>
#include <iterator>
#include <string>
#include <cmath>

#include "H6f.h"
#include "Hbos.h"
#include "Hjets.h"
#include "VVXAnalysis/Producers/interface/Boson.h"

using namespace std;
using namespace edm;
using namespace reco;



class ZZWCombinedGenAnalyzer: public edm::EDAnalyzer {
public:
  ZZWCombinedGenAnalyzer(const ParameterSet& pset) {

    cout << "Type a number \n1: MC history \n2: Real signal, MadGraph pairing \n3:Real signal, real pairing" << endl;
    cin >> num;
    
      }
  
  void analyze(const edm::Event & event, const edm::EventSetup& eventSetup);
  
  virtual void beginJob();
  virtual void endJob(){}
 
  
  //-----------Check of the bosons partons composition-----
  
  int partonsComposition(const Candidate* j_0, const Candidate* j_1) {    
    int j0Id = j_0->pdgId();
    int j1Id = j_1->pdgId();    
    if ( abs(j0Id) < 6 && abs(j1Id ) < 6) {
      if( abs(j0Id + j1Id) == 1 || abs(j0Id + j1Id) == 3 ) {
	if( j0Id % 2 == 0 )  return copysign(24,j0Id);       // W
	else if( j1Id % 2 == 0 )  return copysign(24,j1Id);  // W
	else return 0;
      }
      else if( j0Id + j1Id == 0 ) return 23;                 // Z
      else return 0;                            
    }
    else return 0;
  } 
  //------------------------------------  
  
private:
  H6f* hAll6f;
  H6f* hZZW6f;
  H6f* hBackgr;
  H6f* hAllCut;
  H6f* hZZWCut;
  H6f* hBackgrCut;
  H6f* hZZZ6f;
  H6f* hZZX6f;   
  Hbos* hBosons;
  Hjets* hAllJets;
  Hjets* hBackgrJets;
  Hjets* hCutJets;
  Hjets* hBackgrCutJets;
  TH2F* WZ;
  TH1F* all6fMass;
  TH1F* all6fMassCut;
  TH1F* all6fMassBackgr;
  TH1F* category;
  TH1F* notEv;
  TH1F* lostEvEtaRange;
  int num; 
};



void ZZWCombinedGenAnalyzer::beginJob() {
  hAll6f         = new H6f("All6f");
  hZZW6f         = new H6f("ZZW6f");
  hBackgr        = new H6f("Backgr");
  hAllCut        = new H6f("AllCut");
  hZZWCut        = new H6f("ZZWCut");
  hBackgrCut     = new H6f("BackgrCut");
  hZZZ6f         = new H6f("ZZZ6f");
  hZZX6f         = new H6f("ZZX6f");
  hBosons        = new Hbos("Bosons");
  hAllJets       = new Hjets("AllJets");
  hBackgrJets    = new Hjets("BackgrJets");
  hCutJets       = new Hjets("CutJets");
  hBackgrCutJets = new Hjets("BackgrCutJets");

  edm::Service<TFileService> fileService;
  WZ              = fileService->make<TH2F>("WZ", "WZ", 4, 0., 4., 4, 0., 4.);
  all6fMass       = fileService->make<TH1F>("all6fMass", "all6fMass", 300, 0., 3000.);
  all6fMassCut    = fileService->make<TH1F>("all6fMassCut", "all6fMassCut", 300, 0., 3000.);
  all6fMassBackgr = fileService->make<TH1F>("all6fMassBackgr", "all6fMassBackgr", 300, 0., 3000.);
  category        = fileService->make<TH1F>("category", "category", 7, 0., 7.);
  notEv           = fileService->make<TH1F>("notEv", "notEv", 3, 0., 3.);
  lostEvEtaRange  = fileService->make<TH1F>("lostEvEtaRange", "lostEvEtaRange", 3, 0., 3.);
}




void ZZWCombinedGenAnalyzer::analyze(const Event & event, const EventSetup& eventSetup) { 

  typedef Candidate::LorentzVector LorentzVector;

  std::vector<const reco::Candidate *> theGenZ;
  std::vector<const reco::Candidate *> theGenW;
  std::vector<const reco::Candidate *> theGenl;
  std::vector<const reco::Candidate *> theGenlp;
  std::vector<const reco::Candidate *> theGenlm;
  std::vector<const reco::Candidate *> theGenj;
  std::vector<const reco::Candidate *> theGenq;

  int categoryNum;
  int numMu = 0;
  int numE  = 0;
  
  // Get the collection of gen particles
  edm::Handle<edm::View<reco::Candidate> > genParticles;
  event.getByLabel("genParticles", genParticles);

 
  //------------------ loop over genparticles ---------------------------------------------------------
  for (View<Candidate>::const_iterator p = genParticles->begin(); p != genParticles->end(); ++p) {
    
    if (p->status()==3) {                      // Status 3 particles are the generated ones.
  
      int s_id = p->pdgId();
      int id   = abs(s_id);
      int idx  = std::distance(genParticles->begin(),p);     
      
      bool printIt=true;
      if ( (id < 8 || id == 9) && idx > 5 ) {   // quarks and gluons
	theGenj.push_back(&*p);
      } 
      if ( id < 8 && idx > 5 ) {               // quarks
	theGenq.push_back(&*p);
      } else if ( id==23 ) {                   // Z
	theGenZ.push_back(&*p);
      } else if ( id==24 ) {                   // W
	theGenW.push_back(&*p);
      } else if ( id >= 11 && id <= 16 ) {     // leptons 
	numE  = id == 11 ? numE+1  : numE;
	numMu = id == 13 ? numMu+1 : numMu;
	theGenl.push_back(&*p);
	if(s_id>0) {theGenlp.push_back(&*p);
	} else {theGenlm.push_back(&*p);}
      } else {
	printIt=false;
      }
      if (printIt) {
	cout << "ZZWCombinedGenAnalyzer: idx= "  <<  idx << " \tid= " << p->pdgId() << "\tp4= " << p->p4() <<  endl; 
      }
    }    
  }   
  //------------------ end of loop over genparticles --------------------------------------------------
  

  int leptonCode = 0;
  if( numMu == 2 && numE == 2 ) leptonCode = 2;
  if( (numMu == 4 && numE == 0) || (numMu == 0 && numE == 4) ) leptonCode = 4;

  WZ->Fill(theGenZ.size(), theGenW.size());
  WZ->GetXaxis()->SetTitle("Z");
  WZ->GetYaxis()->SetTitle("W");

  //selection of EVENTS; 2jets, 4leptons
  if ( theGenq.size()!=2 || theGenl.size()!=4 ) { 
    cout << "WARNING: ZZWCombinedGenAnalyzer: " << "\tNumber of jets= " << theGenq.size() << "\tNumber of leptons= "  <<  theGenl.size() << endl;
    notEv->Fill(1);
  } 

  else if ( theGenq.size() == 2 && theGenl.size() == 4 && (leptonCode == 2 || leptonCode == 4) ) {
 
    const float mZ = 91.19;
    const float mW = 80.39;
    
    const Candidate* j0 = theGenq[0];
    const Candidate* j1 = theGenq[1];
    
    if (theGenq[0]->p4().pt() < theGenq[1]->p4().pt()) {
      j0 = theGenq[1];  
      j1 = theGenq[0];  
    }

    Boson *Z0     = new Boson();
    Boson *Z1     = new Boson();
    Boson *W      = new Boson();
    Boson *Z2     = new Boson();

    
    bool has3VCand = false;
    bool has3Z     = false;
    bool isWloose  = false;
    bool isZloose  = false;
    bool isWtight  = false;
    bool isZtight  = false;

    int bosonId = partonsComposition(j0,j1);


    //================Definition of loose particles (mass) 
    if ( fabs((theGenj[0]->p4() + theGenj[1]->p4()).mass() - mW) < 10. ) isWloose = true;

    if ( fabs((theGenj[0]->p4() + theGenj[1]->p4()).mass() - mZ) < 10. ) isZloose = true;

    
    //--------------------1: MC history------------------------------------
    if ( num==1 ) {              
      
      bool LeptonsMotherSelec = true;   
      for(int t=0; t<4; ++t) {
	LeptonsMotherSelec = LeptonsMotherSelec && theGenl[t]->mother()->pdgId() == 23;
      }
      
      if ( isWloose && theGenW.size() == 1 ) isWtight = true;      //definition of tight W (mass + cat)
      if ( isZloose && theGenZ.size() == 3 ) isZtight = true;      //definition of tight Z (mass + cat)
      
      if ( theGenZ.size() >= 2 && LeptonsMotherSelec ) {
	
 	Z0->Setdaughter1(theGenl[0]->p4());
	Z0->Setdaughter2(theGenl[1]->p4());
	Z0->SetbosonId(theGenZ[0]->pdgId());
	
	Z1->Setdaughter1(theGenl[2]->p4());
	Z1->Setdaughter2(theGenl[3]->p4());
	Z1->SetbosonId(theGenZ[1]->pdgId());
	
	if ( abs(theGenl[0]->pdgId()) == 11 ) Z0->SetdaughtersId(1); //u
	if ( abs(theGenl[0]->pdgId()) == 13 ) Z0->SetdaughtersId(2); //e    
	if ( abs(theGenl[2]->pdgId()) == 11 ) Z1->SetdaughtersId(1); //u
	if ( abs(theGenl[2]->pdgId()) == 13 ) Z1->SetdaughtersId(2); //e   
	
	if ( isWtight ) {       

	  W->Setdaughter1(j0->p4());
	  W->Setdaughter2(j1->p4());
	  W->SetbosonId(theGenW[0]->pdgId());
	  W->SetdaughtersId(3);
	  
	  has3VCand = true;
	}
	 
	else if ( isZtight ) {

	  Z2->Setdaughter1(j0->p4());
	  Z2->Setdaughter2(j1->p4());
	  Z2->SetbosonId(theGenZ[2]->pdgId());
	  Z2->SetdaughtersId(3);
	  
	  has3Z = true; 
	}
      }
    }
    

    // -----------------2: Real signal, MadGraph pairing------------------
    else if ( num==2 ) {         
      
      Z0->Setdaughter1(theGenl[0]->p4());            
      Z0->Setdaughter2(theGenl[1]->p4());
      Z0->SetbosonId(23);
      
      Z1->Setdaughter1(theGenl[2]->p4());
      Z1->Setdaughter2(theGenl[3]->p4());
      Z1->SetbosonId(23);
      
      if ( abs(theGenl[0]->pdgId()) == 11 ) Z0->SetdaughtersId(1); //u
      if ( abs(theGenl[0]->pdgId()) == 13 ) Z0->SetdaughtersId(2); //e    
      if ( abs(theGenl[2]->pdgId()) == 11 ) Z1->SetdaughtersId(1); //u
      if ( abs(theGenl[2]->pdgId()) == 13 ) Z1->SetdaughtersId(2); //e   
      
      if ( isWloose && fabs(bosonId) == 24 ) {      //definition of tight W (mass + cat)
	
    	W->Setdaughter1(j0->p4());
    	W->Setdaughter2(j1->p4());
    	W->SetbosonId(bosonId);
	W->SetdaughtersId(3);
   	
	isWtight = true;
	has3VCand = true;
	
      } else if ( isZloose && bosonId == 23 ) {     //definition of tight Z (mass + cat)
	
	Z2->Setdaughter1(j0->p4());
    	Z2->Setdaughter2(j1->p4());
    	Z2->SetbosonId(bosonId);
	Z2->SetdaughtersId(3);
   	
	isZtight = true;
	has3Z = true;     	
      } 	
    }
 
   
    //-----------------3: Real signal, real pairing-----------------------
     else if ( num==3 ) {         
      float minMDiff=99999.;
      if (leptonCode == 4) {
	
	for (int k=0; k<2; ++k) {
	  for (int j=0; j<2; ++j) {
	    float mDiff = fabs((theGenlp[k]->p4() + theGenlm[j]->p4()).mass() - mZ);
	    if ( mDiff < minMDiff ) {
	      minMDiff=mDiff;   
	      
	      Z0->Setdaughter1(theGenlp[k]->p4());            
	      Z0->Setdaughter2(theGenlm[j]->p4());
	      	      
	      Z1->Setdaughter1(theGenlp[(k+1)%2]->p4());            
	      Z1->Setdaughter2(theGenlm[(j+1)%2]->p4());
	        
	      if ( fabs(theGenl[0]->pdgId()) == 11 ) {
		Z0->SetdaughtersId(1); //u
		Z1->SetdaughtersId(1); //u
	      }
	      if ( fabs(theGenl[0]->pdgId()) == 13 ) {
		Z0->SetdaughtersId(2); //e   
		Z1->SetdaughtersId(2); //e   
	      }	
	    }      
	  } 	
	}
      }
      else { 

	for (int z=0; z<2; ++z) {
	  if ( fabs(theGenlp[z]->pdgId()) == fabs(theGenlm[0]->pdgId()) ) { 
	    
	    Z0->Setdaughter1(theGenlp[z]->p4());
	    Z0->Setdaughter2(theGenlm[0]->p4());	  
	    
	    Z1->Setdaughter1(theGenlp[(z+1)%2]->p4());
	    Z1->Setdaughter2(theGenlm[1]->p4());
	    
	    if ( fabs(theGenlm[0]->pdgId()) == 11 ) {
	      Z0->SetdaughtersId(1); //u
	      Z1->SetdaughtersId(2); //e
	    }
	    if ( fabs(theGenlm[0]->pdgId()) == 13 ) {
	      Z0->SetdaughtersId(2); //e
	      Z1->SetdaughtersId(1); //u
	    }

	  }
	}	
      }
    
      Z0->SetbosonId(23);
      Z1->SetbosonId(23);	
      
      if ( Z0->p4().mass() != 0 && Z1->p4().mass() != 0 ) {

	if ( isWloose && fabs(bosonId) == 24 ) {    //definition of tight W (mass + cat)
	  
	  W->Setdaughter1(j0->p4());
	  W->Setdaughter2(j1->p4());
	  W->SetbosonId(bosonId);
	  W->SetdaughtersId(3);
	  
	  isWtight = true;
	  has3VCand = true;   
	   
	} else if ( isZloose && bosonId == 23 ) {   //definition of tight Z (mass + cat)
	  
	  Z2->Setdaughter1(j0->p4());
	  Z2->Setdaughter2(j1->p4());
	  Z2->SetbosonId(bosonId);
	  Z2->SetdaughtersId(3);
	  
	  isZtight = true;
	  has3Z = true; 
	  
	}
      }     
        } 
    
    else {
      abort();
    }
    
    
    //=====================================================================================

    bool hasZZ4l    = fabs(Z0->p4().mass()-mZ) < 10. && fabs(Z1->p4().mass()-mZ) < 10.;
    bool isMySignal = hasZZ4l && fabs(W->p4().mass()-mW) < 10.;
      
    bool passEtaAccLep = true;
    bool passPtAccLep  = true;
    
    LorentzVector p_4l(0.,0.,0.,0.);  
    LorentzVector p_6f(0.,0.,0.,0.);
    for(int i=0; i<4; ++i) {
      p_4l += theGenl[i]->p4();
      passEtaAccLep = passEtaAccLep && fabs(theGenl[i]->p4().eta()) < 2.5;
      passPtAccLep  = passPtAccLep  && (theGenl[i]->p4()).pt() > 7.;
    }

    p_6f = p_4l + j0->p4() + j1->p4(); 
    LorentzVector p_6fBos = Z0->p4daughter1() + Z0->p4daughter2() + Z1->p4daughter1() + Z1->p4daughter2() + W->p4daughter1() + W->p4daughter2();
    float m_6fBos = p_6fBos.mass();
    float m_6f = p_6f.mass();
      
    bool pass6fMass  = m_6f > 300.;
    bool pass6fMassBos = m_6fBos > 300.;
    bool passPtCutJet = j0->p4().pt() > 40. && j1->p4().pt() > 25.;
      
    //////// eta cut for all leptons ///////

    if ( passEtaAccLep ) {
	
      //NO selection-----------------------------
      all6fMass->Fill(p_6f.mass());
      hAll6f->Fill(Z0, Z1, W);
      hAllJets->Filljet(theGenq);
      
      //Cuts-------------------------------------
	
      if ( passPtAccLep && passPtCutJet && pass6fMass ) {
	all6fMassCut->Fill(p_6f.mass());
	hCutJets->Filljet(theGenq);
      }

      if ( passPtAccLep && passPtCutJet && pass6fMassBos ) {
	hAllCut->Fill(Z0, Z1, W);     
      }
      
      //Signal-----------------------------------categoryNum=0----------------------
	
      if ( isMySignal ) {
	categoryNum = 0;
	cout << "SIGNAL: "  << event.id().event() << "\nEvent category: " << categoryNum << endl;
	category->Fill(0);
	hZZW6f->Fill(Z0, Z1, W);

	
	if (passPtAccLep && passPtCutJet && pass6fMass) {
	  hZZWCut->Fill(Z0, Z1, W); 
	}
	
      } 

      //=========Background=========//

      else {
	
	//Background (All-signal)------------------------------------------------------
	all6fMassBackgr->Fill(p_6f.mass());
	hBackgr->Fill(Z0, Z1, W);
	hBackgrJets->Filljet(theGenq);
	
	if ( hasZZ4l ) {
	  
	  //ZZZ events------------------------------ categoryNum= 1 ----------------------    
	  if ( has3Z ) {
	    categoryNum = 1;
	    cout << "ZZZ: " << event.id().event() << "\nEvent category: " << categoryNum << endl;
	    category->Fill(1);
	    hZZZ6f->Fill(Z0, Z1, Z2);
	  }
	}  else {
	  
	  //ZZjj+X events----------------------------- categoryNum= 6 ----------------------    
	  if ( has3Z ) {
	    categoryNum = 6;
	    category->Fill(6);
	    hZZX6f->Fill(Z0, Z1, Z2);
	  }
	  	  
	}
	
	if ( passPtAccLep && passPtCutJet && pass6fMass ) {
	  //Background with cuts
	  hBackgrCut->Fill(Z0, Z1, W);     
	  hBackgrCutJets->Filljet(theGenq);      
	} 
	
      }
      
      //------Bosons mass,Pt,DR-------//
      if ( has3VCand ) {   
	hBosons->FillBos(Z0, Z1, W);
      }  
     
    }
    
    else {
      cout << "LOST_EVENT (out of leptons eta range): " << event.id().event() << endl;
      lostEvEtaRange->Fill(1);
      
    }
    
  }
  
}

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(ZZWCombinedGenAnalyzer);

 
