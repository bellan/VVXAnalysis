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
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include <CommonTools/UtilAlgos/interface/TFileService.h>
#include <DataFormats/HepMCCandidate/interface/GenParticle.h>
#include "DataFormats/JetReco/interface/GenJet.h"

#include "AnalysisDataFormats/CMGTools/interface/PhysicsObject.h"

#include <iostream>
#include <vector>
#include <tuple> 

#include "VVXAnalysis/DataFormats/interface/Boson.h"
#include "VVXAnalysis/Commons/interface/SignalDefinitions.h"
#include "VVXAnalysis/Commons/interface/PhysTools.h"
#include <ZZAnalysis/AnalysisStep/interface/bitops.h>
#include "DataFormats/Math/interface/deltaR.h"

#include <boost/foreach.hpp>
#define foreach BOOST_FOREACH

using namespace std;
using namespace edm;


class ZZGenFilterCategory: public edm::EDFilter {

public:
  

  ZZGenFilterCategory(const ParameterSet& pset)
    : sel_             (pset.getParameter<int>("Topology"))
    , particleStatus_  (pset.getParameter<int>("ParticleStatus"))
    , genLabel_        (pset.getParameter<edm::InputTag>("src"))
    , genJetsLabel_    (pset.getParameter<edm::InputTag>("GenJets"))
  {
    if(particleStatus_ != 1 && particleStatus_ != 3)
      edm::LogError("ZZGenFilterCategory") << "ZZGenFilterCategory: input for signal definition non available!";

    produces<int>();
    produces<std::vector<reco::GenParticle> >("vectorBosons");
    produces<std::vector<reco::GenParticle> >("genJets");
  }
 

  bool filter(edm::Event & event, const edm::EventSetup& eventSetup);
  std::auto_ptr<std::vector<reco::GenParticle> > loadGenBoson(const phys::Boson<phys::Particle> &vb, const reco::GenParticleRefProd &genRefs, std::auto_ptr<std::vector<reco::GenParticle> > outputGenColl);
  
  virtual void beginJob();
  virtual void endJob(){}
    
  //------------------------------------  


private:
  int sel_;
  int particleStatus_;
  edm::InputTag genLabel_;
  edm::InputTag genJetsLabel_;
};

void ZZGenFilterCategory::beginJob() {}

bool ZZGenFilterCategory::filter(Event & event, const EventSetup& eventSetup) { 
  
  // cout << "Run: " << event.id().run() << " event: " << event.id().event() << " LS: " << event.luminosityBlock() << endl;

  std::vector<phys::Particle> genLeptons, genJets;
  
  // Get the collection of gen particles
  edm::Handle<edm::View<reco::Candidate> > genParticles;
  event.getByLabel(genLabel_, genParticles);
  

  switch(particleStatus_){

  case 1:{
    // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    std::vector<phys::Particle> genPhotons;

    //------------------ loop over genparticles ---------------------------------------------------------
    for (View<reco::Candidate>::const_iterator p = genParticles->begin(); p != genParticles->end(); ++p) {
      
      if (false && p->status() == 1){

	if(p->p4().P() != p->p4().P()){
	  cout << "This particle, " << p->pdgId() << ", as NaN in p4 components: " << p->p4().P() << endl;
	  continue;
	}

	int id   = abs(p->pdgId());  
	
	if( id == 11 || id == 13) { 
	  bool fromTau = false;
	  for(unsigned int i = 0; i < p->numberOfMothers(); ++i)
	    if(abs(p->mother(i)->pdgId() == 15)) fromTau = true;
	  if(!fromTau)  genLeptons.push_back(phys::convert(*p)); } // leptons     
	if( id == 22) 
	  genPhotons.push_back(phys::convert(*p)); 
      }
      if (p->status() == 3){
	int id   = abs(p->pdgId());     
	if ( id == 11 || id == 13 ) { genLeptons.push_back(phys::convert(*p)); }// leptons   
      }
    }
    
    //foreach(phys::Particle &lep, genLeptons)
    //  foreach(const phys::Particle &pho, genPhotons)
    //  if(reco::deltaR(lep,pho) < 0.1) lep.setP4(lep.p4()+pho.p4()); // update iteratively, so if there is more than one photon to be associated it will mange it
                                                                    // Note: right no it does not look at the order in which the photons are added, so there could be some
                                                                    // imperfection in the association. One way to overcome this would be to add the nearest one first.
    
    
    //------------------ end of loop over genparticles --------------------------------------------------


    
    // Get gen jets
    edm::Handle<std::vector<cmg::PhysicsObjectWithPtr<edm::Ptr<reco::GenJet> > > > genJetsH;
    event.getByLabel(genJetsLabel_,  genJetsH);
    
    foreach(const cmg::PhysicsObjectWithPtr<edm::Ptr<reco::GenJet> > & jet, *genJetsH)
      if(jet.pt() > 30 && fabs(jet.eta()) < 4.7)
	genJets.push_back(phys::Particle(jet.p4(), phys::Particle::computeCharge(jet.pdgId()), jet.pdgId()));
    
    break;
    // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    }
    
  case 3:
    {
    //------------------ loop over genparticles ---------------------------------------------------------
    for (View<reco::Candidate>::const_iterator p = genParticles->begin(); p != genParticles->end(); ++p) {
      
      if (p->status() == 3) {                      // Status 3 particles are the generated ones.
	
	int id   = abs(p->pdgId());     
	int idx  = std::distance(genParticles->begin(),p);     
	
	if ( (id < 7 || id == 21) && idx > 5 ) { genJets.push_back(phys::convert(*p)); } // quarks and gluons
	
	else if ( id == 11 || id == 13 )       { genLeptons.push_back(phys::convert(*p)); }// leptons        
      }     
    }
    //------------------ end of loop over genparticles --------------------------------------------------
  
    break;
    }
  }

  zz::SignalTopology zzSignalTopology;
  if (particleStatus_ == 1) zzSignalTopology = zz::getSignalTopology(genLeptons, genJets);
  if (particleStatus_ == 3) zzSignalTopology = zz::getSignalTopologyStatus3(genLeptons, genJets);

  // cout << "Topology " <<  std::get<0>(zzSignalTopology) << endl;
  // cout << "Z0 " <<  std::get<1>(zzSignalTopology) << " " << std::get<1>(zzSignalTopology).mass() << endl;
  //if(std::get<1>(zzSignalTopology).id() == 23){
  // cout << std::get<1>(zzSignalTopology).daughter(0) << endl;
  // cout << std::get<1>(zzSignalTopology).daughter(1) << endl;
  // }
  // cout << "Z1 " <<  std::get<2>(zzSignalTopology) << " " << std::get<2>(zzSignalTopology).mass() << endl;
  //if(std::get<2>(zzSignalTopology).id() == 23){
  // cout << std::get<2>(zzSignalTopology).daughter(0) << endl;
  // cout << std::get<2>(zzSignalTopology).daughter(1) << endl;
  // }


  std::auto_ptr<int> output(new int(std::get<0>(zzSignalTopology))); //Topology
  
  event.put(output); 
  
  
  std::auto_ptr<std::vector<reco::GenParticle> > outputGenColl(new std::vector<reco::GenParticle>());
  
  reco::GenParticleRefProd genRefs = event.getRefBeforePut<std::vector<reco::GenParticle> >("vectorBosons");

  // a simple, beautiful loop over the tuple elements is not possible without writing more code than what implemented in this class, for curiosity see
  // http://stackoverflow.com/questions/1198260/iterate-over-tuple, http://stackoverflow.com/questions/18155533/how-to-iterate-through-stdtuple, http://www.metaxcompile.com/blog/2013/03/14/A_cleaner_way_to_do_tuple_iteration.html
  if(     std::get<1>(zzSignalTopology).id()  > 0) outputGenColl = loadGenBoson(std::get<1>(zzSignalTopology), genRefs, outputGenColl); // Z0 --> ll
  if(     std::get<2>(zzSignalTopology).id()  > 0) outputGenColl = loadGenBoson(std::get<2>(zzSignalTopology), genRefs, outputGenColl); // Z1 --> ll
  if(     std::get<3>(zzSignalTopology).id()  > 0) outputGenColl = loadGenBoson(std::get<3>(zzSignalTopology), genRefs, outputGenColl); // Z2 --> ll
  if(     std::get<4>(zzSignalTopology).id()  > 0) outputGenColl = loadGenBoson(std::get<4>(zzSignalTopology), genRefs, outputGenColl); // Z3 --> jj
  if(fabs(std::get<5>(zzSignalTopology).id()) > 0) outputGenColl = loadGenBoson(std::get<5>(zzSignalTopology), genRefs, outputGenColl); // W0 --> lv
  if(fabs(std::get<6>(zzSignalTopology).id()) > 0) outputGenColl = loadGenBoson(std::get<6>(zzSignalTopology), genRefs, outputGenColl); // W1 --> jj

  event.put(outputGenColl,"vectorBosons");


  // load the cleaned GenJets too
  std::auto_ptr<std::vector<reco::GenParticle> > outputGenJetColl(new std::vector<reco::GenParticle>());
  foreach(const phys::Particle& jet, genJets)
    outputGenJetColl->push_back(reco::GenParticle(0, phys::Particle::convert(jet.p4()), reco::GenParticle::Point(0.,0.,0.), jet.id(), 1, true));

  event.put(outputGenJetColl,"genJets");
  
  if (sel_ >= 0) {

    if( ((std::get<0>(zzSignalTopology) ^ sel_) & sel_) == 0)   return true;

    else return false;

  }

  else
    return true;

}

  std::auto_ptr<std::vector<reco::GenParticle> > ZZGenFilterCategory::loadGenBoson(const phys::Boson<phys::Particle> &vb, const reco::GenParticleRefProd &genRefs, std::auto_ptr<std::vector<reco::GenParticle> > outputGenColl){
    
    size_t initialSize = outputGenColl->size();
    
   outputGenColl->push_back(reco::GenParticle(vb.id() == 23 ? 0 : copysign(1, vb.id()), phys::Particle::convert(vb.p4())            , reco::GenParticle::Point(0.,0.,0.), vb.id()            , 3, true));
   outputGenColl->push_back(reco::GenParticle(copysign(1,-1*vb.daughter(0).id())      , phys::Particle::convert(vb.daughter(0).p4()), reco::GenParticle::Point(0.,0.,0.), vb.daughter(0).id(), 3, true));
   outputGenColl->push_back(reco::GenParticle(copysign(1,-1*vb.daughter(1).id())      , phys::Particle::convert(vb.daughter(1).p4()), reco::GenParticle::Point(0.,0.,0.), vb.daughter(1).id(), 3, true));
   outputGenColl->at(initialSize).addDaughter(reco::GenParticleRef(genRefs,initialSize+1));
   outputGenColl->at(initialSize).addDaughter(reco::GenParticleRef(genRefs,initialSize+2));
   
   return outputGenColl;
 }



#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(ZZGenFilterCategory);

 
