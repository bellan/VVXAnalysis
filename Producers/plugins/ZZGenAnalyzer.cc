/**
 *
 *  ZZ Gen Analyzer
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
#include <FWCore/ServiceRegistry/interface/Service.h>

#include <TH1F.h>
#include <TH2F.h>

#include <iostream>
#include <vector>
#include <iterator>
#include <string>
#include <cmath>
#include <tuple> 

#include "VVXAnalysis/DataFormats/interface/Boson.h"
#include "VVXAnalysis/Commons/interface/SignalDefinitions.h"
#include "VVXAnalysis/Commons/interface/PhysTools.h"

using namespace std;
using namespace edm;
using namespace reco;

#include <boost/foreach.hpp>
#define foreach BOOST_FOREACH



class ZZGenAnalyzer: public edm::EDAnalyzer {

public:
  ZZGenAnalyzer(const ParameterSet& pset)
    : topologyLabel_(pset.getParameter<edm::InputTag>("Topology")) {}
					 
  void analyze(const edm::Event & event, const edm::EventSetup& eventSetup);
  
  virtual void beginJob();
  virtual void endJob(){}

private:
  edm::InputTag topologyLabel_;


  TH1I* toptest;
};

void ZZGenAnalyzer::beginJob() {

  edm::Service<TFileService> fileService;

  toptest = fileService->make<TH1I>("toptest", "toptest", 10, 0, 10);
}

void ZZGenAnalyzer::analyze(const Event & event, const EventSetup& eventSetup) { 
  
  typedef Candidate::LorentzVector LorentzVector;
  
  std::vector<phys::Particle> theGenZ;
  std::vector<phys::Particle> theGenW;
  std::vector<phys::Particle> theGenl;
  std::vector<phys::Particle> theGenlp;
  std::vector<phys::Particle> theGenlm;
  std::vector<phys::Particle> theGenj;
  std::vector<phys::Particle> theGenq;

 
  edm::Handle<int> topology;
  event.getByLabel(topologyLabel_, topology);
  
  edm::Handle<edm::View<reco::Candidate> > genParticles;
  event.getByLabel("genParticlesPruned", genParticles);  
  
  //------------------ loop over genparticles ---------------------------------------------------------
  for (View<reco::Candidate>::const_iterator p = genParticles->begin(); p != genParticles->end(); ++p) {
    
    if (p->status()==3) {                      // Status 3 particles are the generated ones.
      
      int id   = abs(p->pdgId());
      int idx  = std::distance(genParticles->begin(),p);     
      
      
      if ( (id < 7 || id == 21) && idx > 5 ) { theGenj.push_back(phys::convert(*p)); } // quarks and gluons
      
      if ( id < 7  && idx > 5 )              { theGenq.push_back(phys::convert(*p)); } // quarks
      
      else if ( id == 23 )                   { theGenZ.push_back(phys::convert(*p)); } // Z
      
      else if ( id == 24 )                   { theGenW.push_back(phys::convert(*p)); } // W 
      
      else if ( id == 11 && id == 13 )       { theGenl.push_back(phys::convert(*p)); }// leptons        
    }
  } //------------------ end of loop over genparticles --------------------------------------------------
  

  // if (theGenj.size() != 2) cout << "WARNING!!!!!!! theGenj size != 2" << endl;

//   else {

//     float mZ = 91.19;
  
//     int numMu = 0; 
//     int numE = 0;
  
//     foreach(const phys::Particle &p, theGenl){
//       numE = abs(p.id()) == 11 ? numE+1 : numE;
//       numMu = abs(p.id()) == 13 ? numMu+1 : numMu;
    
//       if (p.id() > 0) theGenlm.push_back(p); // negative leptons
//       else theGenlp.push_back(p); // positive leptons
//     }
  
  
//     int leptonCode = 0;
//     if( numMu == 2 && numE == 2) leptonCode = 2;
//     if( (numMu == 4 && numE == 0) || (numMu == 0 && numE == 4) ) leptonCode = 4;
 
    
//     std::vector<phys::Boson<phys::Particle> > Z;
    
    
//     foreach(const phys::Particle &p, theGenlp){
//       foreach(const phys::Particle &m, theGenlm){
	
	
// 	phys::Boson<phys::Particle> z;
	
// 	if(abs(p.id()) == abs(m.id())) {
// 	  z.setDaughter(0,p);
// 	  z.setDaughter(1,m);
// 	  z.setId(23);
// 	  Z.push_back(z);
// 	}    
//       }
//     }
      

   //  Z0 = ZZ.first;
//     Z1 = ZZ.second;
//     if( fabs((((ZZ.first).p4()).M()) - mZ) > fabs((((ZZ.second).p4()).M()) - mZ)){
//       Z1 = ZZ.first;
//       Z0 = ZZ.second;
//     }



  //   W.setDaughter(0, phys::Particle(theGenj[0].p4()));
//     W.setDaughter(1, phys::Particle(theGenj[1].p4()));
    
//     phys::Particle j0 = theGenj[0];  
//     phys::Particle j1 = theGenj[1];  
    
//     TLorentzVector p_4l;
//     TLorentzVector p_6f;
//     TLorentzVector p_jj;
    
//     for(int i=0; i<4; ++i) {
//       p_4l += theGenl[i].p4();
//     }  

//     p_jj = j0.p4() + j1.p4();

//     p_6f = p_4l + p_jj; 

//     p_6f = p_4l + j0.p4() + j1.p4(); 
    
//     float m_6f = p_6f.M();
    //float m_4l = p_4l.M();
     
  //   bool passEtaAccLep = true;
    
//     for(int i=0; i<4; ++i) {
//       passEtaAccLep = passEtaAccLep && fabs(theGenl[i].p4().Eta()) < 2.5;
//     }   
    
    
//     if ( passEtaAccLep ) {
      
//       hAll->Fill(Z0,Z1,W);
//       hBosons->FillBos(Z0,Z1,W);
      
//       if (*category == 0)    all6fMass_0->Fill(m_6f);
//       if (*category == 1)    all6fMass_1->Fill(m_6f);
//       if (*category == 2)    all6fMass_2->Fill(m_6f);
//       if (*category == 3)    all6fMass_3->Fill(m_6f);
//       if (*category == 4)    all6fMass_4->Fill(m_6f);
//       if (*category == 5)    all6fMass_5->Fill(m_6f);
//       if (*category == 6)    all6fMass_6->Fill(m_6f);
//       if (*category == 7)    all6fMass_7->Fill(m_6f);
//       if (*category == 8)    all6fMass_8->Fill(m_6f);
//       if (*category == 9)    all6fMass_9->Fill(m_6f);

//       categoryCheck->Fill(*category);
      
//     }  
    
//   }
  
  
 //  if ( theGenq.size()!=2 || theGenl.size()!=4 ) {
//     cout << "WARNING: ZZWGenAnalyzer: " << "\tNumber of jets= " << theGenq.size() << "\tNumber of leptons= " << theGenl.size() << endl;
//     notEv->Fill(1);
//   } 
  
//   else {

//     phys::Particle j0 = theGenq[0];
//     phys::Particle j1 = theGenq[1];
   
//     if (theGenq[0].pt() < theGenq[1].pt()) {
//       j0 = theGenq[1];  
//       j1 = theGenq[0];  
//     }
    
   
    zz::SignalTopology zzSignalTopology = zz::getSignalTopology(theGenl, theGenj);

    phys::Boson<phys::Particle> Z0  = std::get<1>(zzSignalTopology);
    phys::Boson<phys::Particle> Z1  = std::get<2>(zzSignalTopology);
    phys::Boson<phys::Particle> Z2  = std::get<3>(zzSignalTopology);
    phys::Boson<phys::Particle> Z3  = std::get<3>(zzSignalTopology);
    phys::Boson<phys::Particle> W0  = std::get<4>(zzSignalTopology);
    phys::Boson<phys::Particle> W1  = std::get<5>(zzSignalTopology);

    toptest->Fill(*topology);


  //   TLorentzVector p_4l;//(0.,0.,0.,0.);  
//     TLorentzVector p_6f;//(0.,0.,0.,0.);
//     TLorentzVector p_jj;
    
//     for(int i=0; i<4; ++i) {
//       p_4l += theGenl[i].p4();
//     }  
//     p_6f = p_4l + j0.p4() + j1.p4(); 
    
//     p_jj = theGenj[0].p4() + theGenj[1].p4();


//     float m_6f = p_6f.M();
//     float m_4l = p_4l.M();
//     float m_jj = p_jj.M();
     
//     bool passEtaAccLep = true;
    
//     for(int i=0; i<4; ++i) {
//       passEtaAccLep = passEtaAccLep && fabs(theGenl[i].p4().Eta()) < 2.5;
//     }   
   
//     float l1_eta = theGenl[0].p4().Eta();
//     float l2_eta = theGenl[1].p4().Eta();
//     float l3_eta = theGenl[2].p4().Eta();
//     float l4_eta = theGenl[3].p4().Eta();
//     float j1_eta = theGenj[0].p4().Eta();
//     float j2_eta = theGenj[1].p4().Eta();
    

//     hl1_eta->Fill(l1_eta);
//     hl2_eta->Fill(l2_eta);
//     hl3_eta->Fill(l3_eta);
//     hl4_eta->Fill(l4_eta);
//     hj1_eta->Fill(j1_eta);
//     hj2_eta->Fill(j2_eta);
 
 
 
//     //--------Eta Cut-------------------------
    
//         if ( passEtaAccLep ) {
      
//       //==========================HISTOS FILLING=====================//
      
// 	 //	 if (*category == 0)   {
// 	  //all6fMass_0->Fill(m_6f);
// 	 //   all4lMass->Fill(m_4l);
// 	//    if(W.id() != 0) {
// // 	     hAll->Fill(Z0,Z1,W);
// // 	     hBosons->FillBos(Z0,Z1,W);
// 	//      //	categoryCheck->Fill(*category);
// // 	   } else if (Zjj.id() != 0) {
// // 	     hAll->Fill(Z0,Z1,Zjj);
// // 	     hBosons->FillBos(Z0,Z1,Zjj);
// 	     //	categoryCheck->Fill(*category);
// 	  //	   }
// 	   //	 }
// 	 //	 else {
// 	   if (*category == 0)    all6fMass_0->Fill(m_6f);
// 	   if (*category == 1)    all6fMass_1->Fill(m_6f);
// 	   if (*category == 2)    all6fMass_2->Fill(m_6f);
// 	   if (*category == 3)    all6fMass_3->Fill(m_6f);
// 	   if (*category == 4)    all6fMass_4->Fill(m_6f);
// 	   if (*category == 5)    all6fMass_5->Fill(m_6f);
// 	   if (*category == 6)    all6fMass_6->Fill(m_6f);
// 	   if (*category == 7)    all6fMass_7->Fill(m_6f);
// 	   if (*category == 8)    all6fMass_8->Fill(m_6f);
// 	   if (*category == 9)    all6fMass_9->Fill(m_6f);
// 	   //	 }
	 
	 
// 	   hmjj->Fill(m_jj);
// 	   all4lMass->Fill(m_4l);
// 	   h6f->Fill(m_6f);

// 	       }
       
       
//        else {
// 	 cout << "LOST_EVENT (out of leptons eta range): " << event.id().event() << endl;
// 	 lostEvEtaRange->Fill(1);
//     }  
       
//   }
  
}

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(ZZGenAnalyzer);


