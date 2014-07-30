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

#include <iostream>
#include <vector>
#include <tuple> 

#include "VVXAnalysis/DataFormats/interface/Boson.h"
#include "VVXAnalysis/Commons/interface/SignalDefinitions.h"
#include "VVXAnalysis/Commons/interface/PhysTools.h"
#include <ZZAnalysis/AnalysisStep/interface/bitops.h>

using namespace std;
using namespace edm;


class ZZGenFilterCategory: public edm::EDFilter {

public:
  

  ZZGenFilterCategory(const ParameterSet& pset)
    : sel_             (pset.getParameter<int>("Topology"))
    , genLabel_        (pset.getParameter<edm::InputTag>("src")) {
    produces<int>();
    produces<std::vector<reco::GenParticle> >();

  }
 

  bool filter(edm::Event & event, const edm::EventSetup& eventSetup);
  std::auto_ptr<std::vector<reco::GenParticle> > loadGenBoson(const phys::Boson<phys::Particle> &vb, const reco::GenParticleRefProd &genRefs, std::auto_ptr<std::vector<reco::GenParticle> > outputGenColl);
  
  virtual void beginJob();
  virtual void endJob(){}
    
  //------------------------------------  


private:
  int sel_;
  edm::InputTag genLabel_;
  TH1F* topology;

};

void ZZGenFilterCategory::beginJob() {}

bool ZZGenFilterCategory::filter(Event & event, const EventSetup& eventSetup) { 

  std::vector<phys::Particle> theGenZ, theGenW;
  std::vector<phys::Particle> theGenl, theGenj;
  
  // Get the collection of gen particles
  edm::Handle<edm::View<reco::Candidate> > genParticles;
  event.getByLabel(genLabel_, genParticles);
 

  //------------------ loop over genparticles ---------------------------------------------------------
  for (View<reco::Candidate>::const_iterator p = genParticles->begin(); p != genParticles->end(); ++p) {
    
    if (p->status()==3) {                      // Status 3 particles are the generated ones.
      
      int id   = abs(p->pdgId());
      int idx  = std::distance(genParticles->begin(),p);     
      
      
      if ( (id < 7 || id == 21) && idx > 5 ) { theGenj.push_back(phys::convert(*p)); } // quarks and gluons
      
      else if ( id == 23 )                   { theGenZ.push_back(phys::convert(*p)); } // Z
      
      else if ( id == 24 )                   { theGenW.push_back(phys::convert(*p)); } // W 
      
      else if ( id == 11 || id == 13 )       { theGenl.push_back(phys::convert(*p)); }// leptons        
    }
    
  }
  //------------------ end of loop over genparticles --------------------------------------------------

  
  zz::SignalTopology zzSignalTopology = zz::getSignalTopology(theGenl, theGenj);

  // std::auto_ptr<int> output(new int(std::get<0>(zzSignalTopology))); //Topology //FIXME
  std::auto_ptr<int> output(new int(std::get<0>(zzSignalTopology))); //Topology
  
  event.put(output); //To understand...maybe
  
  
  std::auto_ptr<std::vector<reco::GenParticle> > outputGenColl(new std::vector<reco::GenParticle>());
  
  reco::GenParticleRefProd genRefs = event.getRefBeforePut<std::vector<reco::GenParticle> >();

  // a simple, beautiful loop over the tuple elements is not possible without writing more code than what implemented in this class, for curiosity see
  // http://stackoverflow.com/questions/1198260/iterate-over-tuple, http://stackoverflow.com/questions/18155533/how-to-iterate-through-stdtuple, http://www.metaxcompile.com/blog/2013/03/14/A_cleaner_way_to_do_tuple_iteration.html
  if(     std::get<1>(zzSignalTopology).id()  > 0) outputGenColl = loadGenBoson(std::get<1>(zzSignalTopology), genRefs, outputGenColl); // Z0 --> ll
  if(     std::get<2>(zzSignalTopology).id()  > 0) outputGenColl = loadGenBoson(std::get<2>(zzSignalTopology), genRefs, outputGenColl); // Z1 --> ll
  if(     std::get<3>(zzSignalTopology).id()  > 0) outputGenColl = loadGenBoson(std::get<3>(zzSignalTopology), genRefs, outputGenColl); // Z2 --> ll
  if(     std::get<4>(zzSignalTopology).id()  > 0) outputGenColl = loadGenBoson(std::get<4>(zzSignalTopology), genRefs, outputGenColl); // Z3 --> jj
  if(fabs(std::get<5>(zzSignalTopology).id()) > 0) outputGenColl = loadGenBoson(std::get<5>(zzSignalTopology), genRefs, outputGenColl); // W0 --> lv
  if(fabs(std::get<6>(zzSignalTopology).id()) > 0) outputGenColl = loadGenBoson(std::get<6>(zzSignalTopology), genRefs, outputGenColl); // W1 --> jj

  
  event.put(outputGenColl);
  

  if(sel_ >= 0)
    return sel_ == std::get<0>(zzSignalTopology);
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

 
