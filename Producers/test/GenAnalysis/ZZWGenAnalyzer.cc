/**
 *
 *  Gen Analyzer
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
#include <FWCore/Utilities/interface/InputTag.h>

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

#include "VVXAnalysis/DataFormats/interface/Boson.h"
#include "VVXAnalysis/Producers/interface/SignalDefinitionUtilities.h"

using namespace std;
using namespace edm;
using namespace reco;



class ZZWGenAnalyzer: public edm::EDAnalyzer {

public:
  ZZWGenAnalyzer(const ParameterSet& pset):categoryLabel_(pset.getParameter<edm::InputTag>("Category"))  {}

 void analyze(const edm::Event & event, const edm::EventSetup& eventSetup);
  
  virtual void beginJob();
  virtual void endJob(){}

private:

  edm::InputTag categoryLabel_;

  Hbos* hBosons;
  Hbos* hBosonsCut;
  TH1F* all6fMass;
  TH1F* all6fMassCut;  
  TH1F* notEv;
  TH1F* lostEvEtaRange;
};

void ZZWGenAnalyzer::beginJob() {

  hBosons    = new Hbos("Bosons");
  hBosonsCut = new Hbos("BosonsCut");

  edm::Service<TFileService> fileService;
  all6fMass       = fileService->make<TH1F>("all6fMass", "all6fMass", 300, 0., 3000.);
  all6fMassCut    = fileService->make<TH1F>("all6fMassCut", "all6fMassCut", 300, 0., 3000.);
  notEv           = fileService->make<TH1F>("notEv", "notEv", 3, 0., 3.);
  lostEvEtaRange  = fileService->make<TH1F>("lostEvEtaRange", "lostEvEtaRange", 3, 0., 3.);
}

void ZZWGenAnalyzer::analyze(const Event & event, const EventSetup& eventSetup) { 
  
  typedef Candidate::LorentzVector LorentzVector;
  
  std::vector<const reco::Candidate *> theGenZ;
  std::vector<const reco::Candidate *> theGenW;
  std::vector<const reco::Candidate *> theGenl;
  std::vector<const reco::Candidate *> theGenlp;
  std::vector<const reco::Candidate *> theGenlm;
  std::vector<const reco::Candidate *> theGenj;
  std::vector<const reco::Candidate *> theGenq;

  int numMu = 0;
  int numE  = 0;
  
  edm::Handle<int> category;
  event.getByLabel(categoryLabel_, category);

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
      }

 else {
	printIt=false;
      }
      if (printIt) {
	cout << "ZZWGenAnalyzer: idx= "  <<  idx << " \tid= " << p->pdgId() << "\tp4= " << p->p4() <<  endl; 
      }
    }    
  }   
  //------------------ end of loop over genparticles --------------------------------------------------
  
  int leptonCode = 0;
  if( numMu == 2 && numE == 2 ) leptonCode = 2;
  if( (numMu == 4 && numE == 0) || (numMu == 0 && numE == 4) ) leptonCode = 4; 
 
  //selection of EVENTS; 2jets, 4leptons
  if ( theGenq.size()!=2 || theGenl.size()!=4 ) {  
    cout << "WARNING: ZZWGenAnalyzer: " << "\tNumber of jets= " << theGenq.size() << "\tNumber of leptons= "  <<  theGenl.size() << endl;
    notEv->Fill(1);
  } 

  else if ( theGenq.size() == 2 && theGenl.size() == 4 && (leptonCode == 2 || leptonCode == 4) ) {
      
    const float mZ = 91.19;
    
    const Candidate* j0 = theGenq[0];
    const Candidate* j1 = theGenq[1];
    
    if (theGenq[0]->p4().pt() < theGenq[1]->p4().pt()) {
      j0 = theGenq[1];  
      j1 = theGenq[0];  
    }   
   
 
    phys::Boson<phys::Particle> Z0;
    phys::Boson<phys::Particle> Z1;
    phys::Boson<phys::Particle> V ;

    bool passPtAccLep  = true;
    LorentzVector p_4l(0.,0.,0.,0.);  
    LorentzVector p_6f(0.,0.,0.,0.);
    
    for(int i=0; i<4; ++i) {
      p_4l += theGenl[i]->p4();
      passPtAccLep  = passPtAccLep  && (theGenl[i]->p4()).pt() > 7.;
    }  
    p_6f = p_4l + j0->p4() + j1->p4(); 
    
    float m_6f = p_6f.mass();
      
    bool pass6fMass  = m_6f > 300.;
    bool passPtCutJet = j0->p4().pt() > 40. && j1->p4().pt() > 25.;
    
    bool passEtaAccLep = true;
    
    for(int i=0; i<4; ++i) {
      passEtaAccLep = passEtaAccLep && fabs(theGenl[i]->p4().eta()) < 2.5;
    }   
      
    //--------Eta Cut-------------------------
    
    if ( passEtaAccLep ) {
      
      //==========================CATEGORY SELECTION=====================//
      
      if ( *category == 2 ) {
	
	std::pair<phys::Boson<phys::Particle>, phys::Boson<phys::Particle> > ZZ = makeZBosonsFromLeptons(theGenlm, theGenlp, leptonCode, mZ);
	
	Z0 = ZZ.first;
	Z1 = ZZ.second;
	
	V.setDaughter(0,phys::Particle(j0->p4(),phys::Particle::computeCharge(j0->pdgId()), j0->pdgId()));
	V.setDaughter(1,phys::Particle(j1->p4(),phys::Particle::computeCharge(j1->pdgId()), j1->pdgId()));
	
	all6fMass->Fill(m_6f);
	hBosons->FillBos(Z0,Z1,V);
	
	//Cuts-----------------
	if ( passPtAccLep && passPtCutJet && pass6fMass ) {
	  all6fMassCut->Fill(p_6f.mass());
	  hBosonsCut->FillBos(Z0,Z1,V);
	} //-------------------  

      } //==============================================================//
    }
    //--------------------------------------
    
    else {
	cout << "LOST_EVENT (out of leptons eta range): " << event.id().event() << endl;
	lostEvEtaRange->Fill(1);
      }  
      
  }
}


#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(ZZWGenAnalyzer);


