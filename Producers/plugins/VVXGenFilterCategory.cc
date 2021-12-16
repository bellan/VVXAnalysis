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
#include "FWCore/Utilities/interface/EDGetToken.h"
#include <CommonTools/UtilAlgos/interface/TFileService.h>
#include <DataFormats/HepMCCandidate/interface/GenParticle.h>
#include <DataFormats/HepMCCandidate/interface/GenParticleFwd.h>
#include "DataFormats/JetReco/interface/GenJet.h"
#include "DataFormats/JetReco/interface/GenJetCollection.h"

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

class VVXGenFilterCategory: public edm::EDFilter {

public:
  
  
  VVXGenFilterCategory(const ParameterSet& pset)
    : sel_             (pset.getParameter<int>("Topology"))
    , genToken_        (consumes<reco::GenParticleCollection>(pset.getParameter<edm::InputTag>("src")))
    , genJetsToken_    (consumes<reco::GenJetCollection>     (pset.getParameter<edm::InputTag>("GenJets")))
    , genJetsAK8Token_ (consumes<reco::GenJetCollection>     (pset.getParameter<edm::InputTag>("GenJetsAK8")))
  {
    
    produces<int>();
    produces<reco::GenParticleCollection>("vectorBosons");
    produces<reco::GenParticleCollection>("genJets");
    produces<reco::GenParticleCollection>("genJetsAK8");
  }
 

  bool filter(edm::Event & event, const edm::EventSetup& eventSetup);
  void loadGenBoson(const phys::Boson<phys::Particle> &vb, const reco::GenParticleRefProd &genRefs, std::unique_ptr<std::vector<reco::GenParticle> > &outputGenColl);
  
  virtual void beginJob();
  virtual void endJob(){}
    
  //------------------------------------  


private:
  int sel_;
  edm::EDGetTokenT<reco::GenParticleCollection> genToken_;
  edm::EDGetTokenT<reco::GenJetCollection>      genJetsToken_;
  edm::EDGetTokenT<reco::GenJetCollection>      genJetsAK8Token_;
};

void VVXGenFilterCategory::beginJob() {}

bool VVXGenFilterCategory::filter(Event & event, const EventSetup& eventSetup) { 
  
  std::vector<phys::Particle> genLeptons, genJets, genJetsAK8;
  
  // Get the collection of gen particles from hard process (plus prompt photons). It contains leptons (charged and neutrinos) and quarks.
  edm::Handle<reco::GenParticleCollection> genParticles;
  event.getByToken(genToken_, genParticles);
  
  // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  // make the two categories of photon mutually exclusive?
  std::vector<phys::Particle> genPhotons; // this are used for FS correction
  
  
  //------------------ loop over genParticles ---------------------------------------------------------
  foreach (const reco::GenParticle& p, *genParticles) {
    
    if(p.p4().P() != p.p4().P()) continue;
    
    int id = abs(p.pdgId());
    phys::Particle php = phys::convert(p,p.statusFlags().flags_);
    
    // Photons for FSR correction
    if(id == 22)  genPhotons.push_back(php);
    
    // Stable leptons
    if(id == 11 || id == 13 || id == 12 || id == 14) genLeptons.push_back(php);
  }
  
  
  foreach(phys::Particle &lep, genLeptons)
    foreach(const phys::Particle &pho, genPhotons)
    if((abs(lep.id()) == 11 || abs(lep.id()) == 13) && reco::deltaR(lep,pho) < 0.1) lep.setP4(lep.p4()+pho.p4()); // update iteratively, so if there is more than one photon to be associated it will mange it
  // Note: right no it does not look at the order in which the photons are added, so there could be some
  // imperfection in the association. One way to overcome this would be to add the nearest one first.
  
    
  //------------------ end of loop over genparticles --------------------------------------------------
  
  
  
  // Get gen jets
  edm::Handle<reco::GenJetCollection> genJetsH;      event.getByToken(genJetsToken_,  genJetsH);
  
  foreach(const reco::GenJet& jet, *genJetsH)
      genJets.push_back(phys::Particle(jet.p4(), phys::Particle::computeCharge(jet.pdgId()), jet.pdgId()));
  
 
  edm::Handle<reco::GenJetCollection> genJetsAK8H;   event.getByToken(genJetsAK8Token_,  genJetsAK8H);
  
  foreach(const reco::GenJet& jet, *genJetsAK8H)
      genJetsAK8.push_back(phys::Particle(jet.p4(), phys::Particle::computeCharge(jet.pdgId()), jet.pdgId()));
  
  
  
  
  zz::SignalTopology zzSignalTopology;
  zzSignalTopology = zz::getSignalTopology(genLeptons, genJets, genJetsAK8);
  
  
  auto output = std::make_unique<int>(std::get<0>(zzSignalTopology)); //Topology
  
  event.put(std::move(output)); 
  
  
  auto outputGenColl = std::make_unique<reco::GenParticleCollection>();

  
  reco::GenParticleRefProd genRefs = event.getRefBeforePut<reco::GenParticleCollection>("vectorBosons");
  
  // a simple, beautiful loop over the tuple elements is not possible without writing more code than what implemented in this class, for curiosity see
  // http://stackoverflow.com/questions/1198260/iterate-over-tuple, http://stackoverflow.com/questions/18155533/how-to-iterate-through-stdtuple, http://www.metaxcompile.com/blog/2013/03/14/A_cleaner_way_to_do_tuple_iteration.html
  if(     std::get<1>(zzSignalTopology).id()  > 0) loadGenBoson(std::get<1>(zzSignalTopology), genRefs, outputGenColl); // Z0 --> ll
  if(     std::get<2>(zzSignalTopology).id()  > 0) loadGenBoson(std::get<2>(zzSignalTopology), genRefs, outputGenColl); // Z1 --> ll
  if(     std::get<3>(zzSignalTopology).id()  > 0) loadGenBoson(std::get<3>(zzSignalTopology), genRefs, outputGenColl); // Z2 --> ll
  if(     std::get<4>(zzSignalTopology).id()  > 0) loadGenBoson(std::get<4>(zzSignalTopology), genRefs, outputGenColl); // Z3 --> jj
  if(fabs(std::get<5>(zzSignalTopology).id()) > 0) loadGenBoson(std::get<5>(zzSignalTopology), genRefs, outputGenColl); // W0 --> lv
  if(fabs(std::get<6>(zzSignalTopology).id()) > 0) loadGenBoson(std::get<6>(zzSignalTopology), genRefs, outputGenColl); // W1 --> lv
  if(fabs(std::get<7>(zzSignalTopology).id()) > 0) loadGenBoson(std::get<7>(zzSignalTopology), genRefs, outputGenColl); // W2 --> jj

  event.put(std::move(outputGenColl), "vectorBosons");
  
  
  // load the cleaned GenJets too
  auto outputGenJetColl = std::make_unique<reco::GenParticleCollection>();
  foreach(const phys::Particle& jet, genJets)
    outputGenJetColl->push_back(reco::GenParticle(0, phys::Particle::convert(jet.p4()), reco::GenParticle::Point(0.,0.,0.), jet.id(), 1, true));

  event.put(std::move(outputGenJetColl),"genJets");


  auto outputGenJetAK8Coll = std::make_unique<reco::GenParticleCollection>();
  foreach(const phys::Particle& jet, genJetsAK8)
    outputGenJetAK8Coll->push_back(reco::GenParticle(0, phys::Particle::convert(jet.p4()), reco::GenParticle::Point(0.,0.,0.), jet.id(), 1, true));

  event.put(std::move(outputGenJetAK8Coll),"genJetsAK8");

  
  
  if (sel_ >= 0) {
    
    if( ((std::get<0>(zzSignalTopology) ^ sel_) & sel_) == 0)   return true;
    
    else return false;
    
  }
  
  else
    return true;
  
}

void VVXGenFilterCategory::loadGenBoson(const phys::Boson<phys::Particle> &vb, const reco::GenParticleRefProd &genRefs, std::unique_ptr<reco::GenParticleCollection> &outputGenColl){
    
  size_t initialSize = outputGenColl->size();
  
  outputGenColl->push_back(reco::GenParticle(vb.id() == 23 ? 0 : copysign(1, vb.id()), phys::Particle::convert(vb.p4())            , reco::GenParticle::Point(0.,0.,0.), vb.id()            , 3, true));
  outputGenColl->push_back(reco::GenParticle(copysign(1,-1*vb.daughter(0).id())      , phys::Particle::convert(vb.daughter(0).p4()), reco::GenParticle::Point(0.,0.,0.), vb.daughter(0).id(), 3, true));
  outputGenColl->push_back(reco::GenParticle(copysign(1,-1*vb.daughter(1).id())      , phys::Particle::convert(vb.daughter(1).p4()), reco::GenParticle::Point(0.,0.,0.), vb.daughter(1).id(), 3, true));
  outputGenColl->at(initialSize).addDaughter(reco::GenParticleRef(genRefs,initialSize+1));
  outputGenColl->at(initialSize).addDaughter(reco::GenParticleRef(genRefs,initialSize+2));
}



#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(VVXGenFilterCategory);

 
