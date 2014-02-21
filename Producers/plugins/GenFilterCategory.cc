/**
 *
 *  Filter of the generator history.
 *
 */


#include <FWCore/Framework/interface/Frameworkfwd.h>
#include <FWCore/Framework/interface/EDFilter.h>
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

#include "VVXAnalysis/Producers/interface/Boson.h"
#include "VVXAnalysis/Producers/interface/SignalDefinitionUtilities.h"

using namespace std;
using namespace edm;
using namespace reco;


class GenFilterCategory: public edm::EDFilter {

public:
  GenFilterCategory(const ParameterSet& pset)
    : sel_     (pset.getParameter<int>("Category"))
    , num      (pset.getParameter<int>("SignalDefinition"))
    , genLabel_(pset.getParameter<edm::InputTag>("src")) {
    produces<int>();

  }
  
  bool filter(edm::Event & event, const edm::EventSetup& eventSetup);
  
  virtual void beginJob();
  virtual void endJob(){}
  
  //-----------FUNCTION: Check of the bosons partons composition---------
  
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
  int sel_;
  int num;
  edm::InputTag genLabel_;
  TH1F* category;

};

void GenFilterCategory::beginJob() {

  //edm::Service<TFileService> fileService;
  // category        = fileService->make<TH1F>("category", "category", 10, 0., 10.);
}

 bool GenFilterCategory::filter(Event & event, const EventSetup& eventSetup) { 

  typedef Candidate::LorentzVector LorentzVector;

  std::vector<const reco::Candidate *> theGenZ;
  std::vector<const reco::Candidate *> theGenW;
  std::vector<const reco::Candidate *> theGenl;
  std::vector<const reco::Candidate *> theGenlp;
  std::vector<const reco::Candidate *> theGenlm;
  std::vector<const reco::Candidate *> theGenj;
  std::vector<const reco::Candidate *> theGenq;

  int categoryNum = 999; 
  int numMu       = 0;
  int numE        = 0;
  
  // Get the collection of gen particles
  edm::Handle<edm::View<reco::Candidate> > genParticles;
  event.getByLabel(genLabel_, genParticles);
 
  //------------------ loop over genparticles ---------------------------------------------------------
  for (View<Candidate>::const_iterator p = genParticles->begin(); p != genParticles->end(); ++p) {
    
    if (p->status()==3) {                      // Status 3 particles are the generated ones.
  
      int s_id = p->pdgId();
      int id   = abs(s_id);
      int idx  = std::distance(genParticles->begin(),p);     
      
      if ( (id < 7 || id == 21) && idx > 5 ) {   // quarks and gluons
	theGenj.push_back(&*p);
      }
      if ( id < 7  && idx > 5 ) {               // quarks
	theGenq.push_back(&*p);
      } else if ( id==23 ) {                    // Z
	theGenZ.push_back(&*p);
      } else if ( id==24 ) {                    // W
	theGenW.push_back(&*p);
      } else if ( id >= 11 && id <= 16 ) {      // leptons 
	numE  = id == 11 ? numE+1  : numE;
	numMu = id == 13 ? numMu+1 : numMu;
	theGenl.push_back(&*p);
	if(s_id>0) {theGenlp.push_back(&*p);
	} else {theGenlm.push_back(&*p);}
      } 
    }    
  }   
  //------------------ end of loop over genparticles --------------------------------------------------
  

  int leptonCode = 0;
  if( numMu == 2 && numE == 2 ) leptonCode = 2;
  if( (numMu == 4 && numE == 0) || (numMu == 0 && numE == 4) ) leptonCode = 4;

  if ( theGenq.size() == 2 && theGenl.size() == 4 && (leptonCode == 2 || leptonCode == 4) ) {

    const float mZ = 91.19;
    const float mW = 80.39;
    
    const Candidate* j0 = theGenq[0];
    const Candidate* j1 = theGenq[1];
    
    if ( theGenq[0]->p4().pt() < theGenq[1]->p4().pt() ) {
      j0 = theGenq[1];  
      j1 = theGenq[0];  
    }

    Boson *Z0 = new Boson();
    Boson *Z1 = new Boson();
    Boson *W  = new Boson();
    Boson *Z2 = new Boson();

    
    bool has3Z     = false;
    bool isWloose  = false;
    bool isZloose  = false;
    bool isWtight  = false;
    bool isZtight  = false;

    int bosonId = partonsComposition(j0,j1);
    bool qqPassMWwindow = fabs((j0->p4() + j1->p4()).mass() - mW) < 10;
    bool qqPassMZwindow = fabs((j0->p4() + j1->p4()).mass() - mZ) < 10;

    //================Definition of loose particles (mass) =======================

    
    for(uint i = 0;  i < theGenj.size()-1; ++i)
      for(uint j = i+1;  j < theGenj.size(); ++j){
	if ( fabs((theGenj[i]->p4() + theGenj[j]->p4()).mass() - mW) < 10. ) isWloose = true;  
	if ( fabs((theGenj[i]->p4() + theGenj[j]->p4()).mass() - mZ) < 10. ) isZloose = true;
      }
     

    
    //--------------------1: MC history------------------------------------
    if ( num==1 ) {              
      
      bool LeptonsMotherSelec = true;   
      for(int t=0; t<4; ++t) {
	LeptonsMotherSelec = LeptonsMotherSelec && theGenl[t]->mother()->pdgId() == 23;
      }
      
      if ( (isWloose || qqPassMWwindow) && theGenW.size() == 1) isWtight = true;      //definition of tight W (mass + cat)
      if ( (isZloose || qqPassMZwindow) && theGenZ.size() == 3) isZtight = true;      //definition of tight Z (mass + cat)
      
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
      
      if ( (isWloose || qqPassMWwindow) && fabs(bosonId) == 24 ) {      //definition of tight W (mass + cat)
	
    	W->Setdaughter1(j0->p4());
    	W->Setdaughter2(j1->p4());
    	W->SetbosonId(bosonId);
	W->SetdaughtersId(3);
   	
	isWtight = true;
	
      } else if ( (isZloose || qqPassMZwindow) && bosonId == 23 ) {     //definition of tight Z (mass + cat)
	
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

      std::pair<Boson*,Boson*> ZZ = makeZbosonsFromLeptons(theGenl, theGenlm, theGenlp, leptonCode, mZ);
      
      Z0 = ZZ.first;
      Z1 = ZZ.second;

      if ( Z0->p4().mass() != 0 && Z1->p4().mass() != 0 ) {

	if ( (isWloose || qqPassMWwindow)  && fabs(bosonId) == 24 ) {    //definition of tight W (mass + cat)
	  
	  W->Setdaughter1(j0->p4());
	  W->Setdaughter2(j1->p4());
	  W->SetbosonId(bosonId);
	  W->SetdaughtersId(3);
	  
	  if (qqPassMWwindow) isWtight = true;  
	  
	} else if ( (isZloose || qqPassMZwindow)  && bosonId == 23 ) {   //definition of tight Z (mass + cat)
	  
	  Z2->Setdaughter1(j0->p4());
	  Z2->Setdaughter2(j1->p4());
	  Z2->SetbosonId(bosonId);
	  Z2->SetdaughtersId(3);
	  
	  if(qqPassMZwindow){
	    isZtight = true;
	    has3Z = true; 
	  }
	}
      }     
    } 
    
    else {
      abort();
    }

   
    
    //=====================================================================================

    bool hasZZ4l    = fabs(Z0->p4().mass()-mZ) < 10. && fabs(Z1->p4().mass()-mZ) < 10.;
    bool isMySignal = hasZZ4l && isWtight;
      
    bool passEtaAccLep = true;
      
    for(int i=0; i<4; ++i) {
      passEtaAccLep = passEtaAccLep && fabs(theGenl[i]->p4().eta()) < 2.5;
    }
      
 
    //////// eta cut for all leptons ///////
    
    if ( passEtaAccLep ) {
      
      //cout << "CHECK: Z0 MASS =  " << Z0->p4().mass() << endl;

      //Signal: ZZW---------------------------------------categoryNum=0---------------------
      if ( isMySignal ) {
	categoryNum = 0;
	//cout << "SIGNAL: "  << event.id().event() << "\nEvent category: " << categoryNum << endl;
	//category->Fill(0);

      } 

      //=========Background=========//

      else {     
	   
	if( hasZZ4l ){

	  //ZZZ ----------------------------------------- categoryNum = 1 ------------------------ 
	  if ( has3Z ) {
	    categoryNum = 1;
	    //cout << "ZZZ: " << event.id().event() << "\nEvent category: " << categoryNum << endl;
	    //category->Fill(1);
	  }
	  
	  //ZZWloose ------------------------------------ categoryNum = 2 -----------------------
	  else if ( !has3Z && isWloose ) {
	    categoryNum = 2;
	    //category->Fill(2);
	  }
	  
	  //ZZZloose ------------------------------------ categoryNum = 3 -----------------------
	  else if ( !has3Z && !isWloose && isZloose ) {
	    categoryNum = 3;
	    //category->Fill(3);
	  }

	  //ZZ+X ---------------------------------------- categoryNum = 4 ------------------------
	  else {	
	    categoryNum = 4;
	    //category->Fill(4);
	  }
	}

	else {
	  
	  //WZ+X ---------------------------------------- categoryNum = 5 ------------------------
	  if ( isWtight ) {
	    categoryNum = 5;
	    //category->Fill(5);
	  }
	  
	  //ZZjj+X -------------------------------------- categoryNum = 6 ------------------------
	  else if ( has3Z ) {
	    categoryNum = 6;
	    //category->Fill(6);
	  }
	  
	  //ZWloose+X ----------------------------------- categoryNum = 7 ------------------------
	  else if ( !isWtight && !has3Z && isWloose ) {
	    categoryNum = 7;
	    //category->Fill(7);
	  }
	  
	  //ZZjj+X -------------------------------------- categoryNum = 8 ------------------------
	  else if ( !has3Z && !isWloose && isZloose ) {
	    categoryNum = 8;
	    //category->Fill(8);
	  }
	  
	  //Z+X+Y --------------------------------------- categoryNum = 9 ------------------------
	  else {
	    categoryNum = 9;
	    //category->Fill(9);
	  }
	  
	}		
      }      
    }
  }
  
  std::auto_ptr<int> output(new int(categoryNum));
  event.put(output);
  
  if(sel_ >= 0)
    return sel_ == categoryNum;
  else
    return true;
 }

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(GenFilterCategory);

 
