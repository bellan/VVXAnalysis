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

#include "VVXAnalysis/DataFormats/interface/Boson.h"
#include "VVXAnalysis/Producers/interface/SignalDefinitionUtilities.h"

using namespace std;
using namespace edm;
using namespace reco;


const float MZ = 91.19;
const float MW = 80.39;


class GenFilterCategory: public edm::EDFilter {

public:
  GenFilterCategory(const ParameterSet& pset)
    : sel_     (pset.getParameter<int>("Category"))
    , num      (pset.getParameter<int>("SignalDefinition"))
    , genLabel_(pset.getParameter<edm::InputTag>("src")) {
    produces<int>();
    produces<std::vector<reco::GenParticle> >();

  }
 

  bool filter(edm::Event & event, const edm::EventSetup& eventSetup);
  std::auto_ptr<std::vector<reco::GenParticle> > loadGenBoson(const phys::Boson<phys::Particle> &vb, const GenParticleRefProd &genRefs, std::auto_ptr<std::vector<reco::GenParticle> > outputGenColl);
  
  virtual void beginJob();
  virtual void endJob(){}
    
  //------------------------------------  


private:
  int sel_;
  int num;
  edm::InputTag genLabel_;
  TH1F* category;

};

void GenFilterCategory::beginJob() {}

bool GenFilterCategory::filter(Event & event, const EventSetup& eventSetup) { 

  std::vector<phys::Particle> theGenZ, theGenW;
  std::vector<phys::Particle> theGenl, theGenlp, theGenlm, theGenj, theGenq;

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
      
      
      if ( (id < 7 || id == 21) && idx > 5 ) { theGenj.push_back(convert(*p)); } // quarks and gluons

      if ( id < 7  && idx > 5 )              { theGenq.push_back(convert(*p)); } // quarks
      
      else if ( id == 23 )                   { theGenZ.push_back(convert(*p)); } // Z
   
      else if ( id == 24 )                   { theGenW.push_back(convert(*p)); } // W 

      else if ( id >= 11 && id <= 16 )       {                                   // leptons  
	numE  = id == 11 ? numE+1  : numE; numMu = id == 13 ? numMu+1 : numMu;
	                                       theGenl.push_back (convert(*p));
	if (s_id>0)                          { theGenlp.push_back(convert(*p));} // positive leptons                                          
	else                                 { theGenlm.push_back(convert(*p));} // negative leptons                                          
      } 
    }
    
  }   
  //------------------ end of loop over genparticles --------------------------------------------------
  
  int leptonCode = 0;
  if(  numMu == 2 && numE == 2)                                leptonCode = 2;
  if( (numMu == 4 && numE == 0) || (numMu == 0 && numE == 4) ) leptonCode = 4;

  phys::Boson<phys::Particle> Z0, Z1, Z2, W;

  // FIXME!!!!
  if ( theGenq.size() == 2 && theGenl.size() == 4 && (leptonCode == 2 || leptonCode == 4) ) {

    phys::Particle q0 = theGenq[0], q1 = theGenq[1];
    
    if ( q0.pt() < q1.pt() ) { q0 = theGenq[1];  q1 = theGenq[0]; }
    
    bool has3Z = false, isWloose  = false, isZloose  = false, isWtight  = false, isZtight  = false;

    // FIXME!!!!
    int bosonId = makeVBosonsFromIds(q0.id(), q1.id());

    // FIXME!!!!
    bool qqPassMWwindow = fabs((q0.p4() + q1.p4()).M() - MW) < 10;
    bool qqPassMZwindow = fabs((q0.p4() + q1.p4()).M() - MZ) < 10;

    //================Definition of loose particles (mass) =======================

    
    for(uint i = 0;  i < theGenj.size()-1; ++i)
      for(uint j = i+1;  j < theGenj.size(); ++j){
	if ( fabs((theGenj[i].p4() + theGenj[j].p4()).M() - MW) < 10. ) isWloose = true;  
	if ( fabs((theGenj[i].p4() + theGenj[j].p4()).M() - MZ) < 10. ) isZloose = true;
      }
     
    
    //--------------------1: MC history------------------------------------
    if ( num==1 ) {              
      
      bool LeptonsMotherSelec = true;   
      for(int t=0; t<4; ++t) LeptonsMotherSelec = LeptonsMotherSelec && theGenl[t].motherId() == 23;
      
      
      if ( (isWloose || qqPassMWwindow) && theGenW.size() == 1) isWtight = true;      //definition of tight W (mass + cat)
      if ( (isZloose || qqPassMZwindow) && theGenZ.size() == 3) isZtight = true;      //definition of tight Z (mass + cat)
      
      if ( theGenZ.size() >= 2 && LeptonsMotherSelec ) {

	Z0.setDaughter(0, theGenl[0]);
	Z0.setDaughter(1, theGenl[1]);
	Z0.setId(theGenZ[0].id());

 	Z1.setDaughter(0, theGenl[2]);
	Z1.setDaughter(1, theGenl[3]);
	Z1.setId(theGenZ[1].id());


	if ( isWtight ) {       

	  W.setDaughter(0, q0);
	  W.setDaughter(1, q1);
	  W.setId(theGenW[0].id());
	  
	}
	 
	else if ( isZtight ) {

	  Z2.setDaughter(0, q0);
	  Z2.setDaughter(1, q1);
	  Z2.setId(theGenZ[2].id());
	  
	  has3Z = true; 
	}
      }
    }
    

    // -----------------2: Real signal, MadGraph pairing------------------
    else if ( num==2 ) {         
      
      Z0.setDaughter(0, theGenl[0]);
      Z0.setDaughter(1, theGenl[1]);
      Z0.setId(theGenZ[0].id());
      
      Z1.setDaughter(0, theGenl[2]);
      Z1.setDaughter(1, theGenl[3]);
      Z1.setId(theGenZ[1].id());
      
      if ( (isWloose || qqPassMWwindow) && fabs(bosonId) == 24 ) {      //definition of tight W (mass + cat)
	
    	W.setDaughter(0, q0);
	W.setDaughter(1, q1);
	W.setId(bosonId);
   	
	isWtight = true;
	
      } else if ( (isZloose || qqPassMZwindow) && bosonId == 23 ) {     //definition of tight Z (mass + cat)
	
	Z2.setDaughter(0, q0);
	Z2.setDaughter(1, q1);
	Z2.setId(bosonId);
   	
	isZtight = true;
	has3Z = true;     	
      } 	
    }
 
   
    //-----------------3: Real signal, real pairing-----------------------
    else if ( num==3 ) {         
      std::pair<phys::Boson<phys::Particle>, phys::Boson<phys::Particle> > ZZ = makeZBosonsFromLeptons(theGenlm, theGenlp, leptonCode, MZ);

      Z0 = ZZ.first;
      Z1 = ZZ.second;

      if ( Z0.p4().M() != 0 && Z1.p4().M() != 0 ) {

	if ( (isWloose || qqPassMWwindow)  && fabs(bosonId) == 24 ) {    //definition of tight W (mass + cat)
	  
	  W.setDaughter(0, q0);
	  W.setDaughter(1, q1);
	  W.setId(bosonId);
  
	  if (qqPassMWwindow) isWtight = true;  
	  
	} else if ( (isZloose || qqPassMZwindow)  && bosonId == 23 ) {   //definition of tight Z (mass + cat)
	  
	  Z2.setDaughter(0, q0);
	  Z2.setDaughter(1, q1);
	  Z2.setId(bosonId);

	  if(qqPassMZwindow){
	    isZtight = true;
	    has3Z = true; 
	  }
	}
      }     
    } 
    
    else { cout << "*** Signal definition not found! ***" << endl; abort(); }

    //=====================================================================================

    bool hasZZ4l    = fabs(Z0.p4().M()-MZ) < 10. && fabs(Z1.p4().M()-MZ) < 10.;
    bool isMySignal = hasZZ4l && isWtight;
      
    bool passEtaAccLep = true;
      
    for(int i=0; i<4; ++i) {
      passEtaAccLep = passEtaAccLep && fabs(theGenl[i].eta()) < 2.5;
    }
      
 
    //////// eta cut for all leptons ///////
    
    if ( passEtaAccLep) {

      //Signal: ZZW---------------------------------------categoryNum=0---------------------
      if ( isMySignal ) {
	categoryNum = 0;
	//cout << "SIGNAL: "  << event.id().event() << "\nEvent category: " << categoryNum << endl;
      } 

      //=========Background=========//

      else {     
	   
	if( hasZZ4l ){

	  //ZZZ ----------------------------------------- categoryNum = 1 ------------------------ 
	  if ( has3Z ) {
	    categoryNum = 1;
	    //cout << "ZZZ: " << event.id().event() << "\nEvent category: " << categoryNum << endl;
	  }
	  
	  //ZZWloose ------------------------------------ categoryNum = 2 -----------------------
	  else if ( !has3Z && isWloose ) {
	    categoryNum = 2;
	  }
	  
	  //ZZZloose ------------------------------------ categoryNum = 3 -----------------------
	  else if ( !has3Z && !isWloose && isZloose ) {
	    categoryNum = 3;
	  }

	  //ZZ+X ---------------------------------------- categoryNum = 4 ------------------------
	  else {	
	    categoryNum = 4;
	  }
	}

	else {
	  
	  //WZ+X ---------------------------------------- categoryNum = 5 ------------------------
	  if ( isWtight ) {
	    categoryNum = 5;
	  }
	  
	  //ZZjj+X -------------------------------------- categoryNum = 6 ------------------------
	  else if ( has3Z ) {
	    categoryNum = 6;
	  }
	  
	  //ZWloose+X ----------------------------------- categoryNum = 7 ------------------------
	  else if ( !isWtight && !has3Z && isWloose ) {
	    categoryNum = 7;
	  }
	  
	  //ZZjj+X -------------------------------------- categoryNum = 8 ------------------------
	  else if ( !has3Z && !isWloose && isZloose ) {
	    categoryNum = 8;
	  }
	  
	  //Z+X+Y --------------------------------------- categoryNum = 9 ------------------------
	  else {
	    categoryNum = 9;
	  }
	  
	}		
      }      
    }
  }
 
  std::auto_ptr<int> output(new int(categoryNum));
  event.put(output);


  std::auto_ptr<std::vector<reco::GenParticle> > outputGenColl(new std::vector<reco::GenParticle>());

  GenParticleRefProd genRefs = event.getRefBeforePut<std::vector<reco::GenParticle> >();

  if(Z0.id() > 0) outputGenColl = loadGenBoson(Z0, genRefs, outputGenColl);
  if(Z1.id() > 0) outputGenColl = loadGenBoson(Z1, genRefs, outputGenColl);
  if(Z2.id() > 0) outputGenColl = loadGenBoson(Z2, genRefs, outputGenColl);
  if(fabs(W.id()) > 0)  outputGenColl = loadGenBoson(W, genRefs, outputGenColl);
  
  event.put(outputGenColl);
  
  if(sel_ >= 0)
    return sel_ == categoryNum;
  else
    return true;
 }


std::auto_ptr<std::vector<reco::GenParticle> > GenFilterCategory::loadGenBoson(const phys::Boson<phys::Particle> &vb, const GenParticleRefProd &genRefs, std::auto_ptr<std::vector<reco::GenParticle> > outputGenColl){

   size_t initialSize = outputGenColl->size();

   outputGenColl->push_back(reco::GenParticle(vb.id() == 23 ? 0 : copysign(1, vb.id()), phys::Particle::convert(vb.p4())            , GenParticle::Point(0.,0.,0.), vb.id()            , 3, true));
   outputGenColl->push_back(reco::GenParticle(copysign(1,-1*vb.daughter(0).id())      , phys::Particle::convert(vb.daughter(0).p4()), GenParticle::Point(0.,0.,0.), vb.daughter(0).id(), 3, true));
   outputGenColl->push_back(reco::GenParticle(copysign(1,-1*vb.daughter(1).id())      , phys::Particle::convert(vb.daughter(1).p4()), GenParticle::Point(0.,0.,0.), vb.daughter(1).id(), 3, true));
   outputGenColl->at(initialSize).addDaughter(GenParticleRef(genRefs,initialSize+1));
   outputGenColl->at(initialSize).addDaughter(GenParticleRef(genRefs,initialSize+2));
   
   return outputGenColl;
 }






#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(GenFilterCategory);

 
