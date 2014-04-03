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

using namespace std;
using namespace edm;


class GenFilterCategory: public edm::EDFilter {

public:
  

  GenFilterCategory(const ParameterSet& pset)
    : sel_             (pset.getParameter<int>("Category"))
    , signalDefinition_(pset.getParameter<int>("SignalDefinition"))
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
  int signalDefinition_;
  edm::InputTag genLabel_;
  TH1F* category;

};

void GenFilterCategory::beginJob() {}

bool GenFilterCategory::filter(Event & event, const EventSetup& eventSetup) { 

  std::vector<phys::Particle> theGenZ, theGenW;
  std::vector<phys::Particle> theGenl, theGenj, theGenq;
  
  // Get the collection of gen particles
  edm::Handle<edm::View<reco::Candidate> > genParticles;
  event.getByLabel(genLabel_, genParticles);
 
  //------------------ loop over genparticles ---------------------------------------------------------
  for (View<reco::Candidate>::const_iterator p = genParticles->begin(); p != genParticles->end(); ++p) {
    
    if (p->status()==3) {                      // Status 3 particles are the generated ones.
  
      int id   = abs(p->pdgId());
      int idx  = std::distance(genParticles->begin(),p);     
      
      
      if ( (id < 7 || id == 21) && idx > 5 ) { theGenj.push_back(phys::convert(*p)); } // quarks and gluons

      if ( id < 7  && idx > 5 )              { theGenq.push_back(phys::convert(*p)); } // quarks
      
      else if ( id == 23 )                   { theGenZ.push_back(phys::convert(*p)); } // Z
   
      else if ( id == 24 )                   { theGenW.push_back(phys::convert(*p)); } // W 

      else if ( id >= 11 && id <= 16 )       { theGenl.push_back(phys::convert(*p)); }// leptons        
    }
    //     if(p->status()==1 && p->mother()->pdgId() == 23  && (abs(p->pdgId()) == 11 || abs(p->pdgId()) == 13)){
    //       numE_st1  = abs(p->pdgId()) == 11 ? numE_st1+1  : numE_st1; numMu_st1 = abs(p->pdgId()) == 13 ? numMu_st1+1 : numMu_st1;
    //       theGenl_st1.push_back (convert(*p));
    //       if (p->pdgId() > 0)                  { theGenlp_st1.push_back(phys::convert(*p));} // positive leptons                                          
    //       else                                 { theGenlm_st1.push_back(phys::convert(*p));} // negative leptons      
    //     }
  }
  //   if(theGenl.size() == 0){
  //     theGenl  = theGenl_st1;
  //     theGenlp = theGenlp_st1;
  //     theGenlm = theGenlm_st1;
  //     numMu    = numMu_st1;
  //     numE    = numE_st1;
  //   }
  
  //------------------ end of loop over genparticles --------------------------------------------------

  
  zzw::GenTopology zzwGenTopology = zzw::getGenTopology(signalDefinition_, theGenl, theGenj, theGenZ, theGenW);

  std::auto_ptr<int> output(new int(std::get<0>(zzwGenTopology))); // GenCategory
  event.put(output);
  
  
  std::auto_ptr<std::vector<reco::GenParticle> > outputGenColl(new std::vector<reco::GenParticle>());
  
  reco::GenParticleRefProd genRefs = event.getRefBeforePut<std::vector<reco::GenParticle> >();

  // a simple, beautiful loop over the tuple elements is not possible without writing more code than what implemented in this class, for curiosity see
  // http://stackoverflow.com/questions/1198260/iterate-over-tuple, http://stackoverflow.com/questions/18155533/how-to-iterate-through-stdtuple, http://www.metaxcompile.com/blog/2013/03/14/A_cleaner_way_to_do_tuple_iteration.html
  if(     std::get<1>(zzwGenTopology).id()  > 0) outputGenColl = loadGenBoson(std::get<1>(zzwGenTopology), genRefs, outputGenColl); // Z0 --> ll
  if(     std::get<2>(zzwGenTopology).id()  > 0) outputGenColl = loadGenBoson(std::get<2>(zzwGenTopology), genRefs, outputGenColl); // Z1 --> ll
  if(     std::get<3>(zzwGenTopology).id()  > 0) outputGenColl = loadGenBoson(std::get<3>(zzwGenTopology), genRefs, outputGenColl); // Z2 --> jj
  if(fabs(std::get<4>(zzwGenTopology).id()) > 0) outputGenColl = loadGenBoson(std::get<4>(zzwGenTopology), genRefs, outputGenColl); //  W --> jj
  
  event.put(outputGenColl);
  
  if(sel_ >= 0)
    return sel_ == std::get<0>(zzwGenTopology);
  else
    return true;
}


std::auto_ptr<std::vector<reco::GenParticle> > GenFilterCategory::loadGenBoson(const phys::Boson<phys::Particle> &vb, const reco::GenParticleRefProd &genRefs, std::auto_ptr<std::vector<reco::GenParticle> > outputGenColl){

   size_t initialSize = outputGenColl->size();

   outputGenColl->push_back(reco::GenParticle(vb.id() == 23 ? 0 : copysign(1, vb.id()), phys::Particle::convert(vb.p4())            , reco::GenParticle::Point(0.,0.,0.), vb.id()            , 3, true));
   outputGenColl->push_back(reco::GenParticle(copysign(1,-1*vb.daughter(0).id())      , phys::Particle::convert(vb.daughter(0).p4()), reco::GenParticle::Point(0.,0.,0.), vb.daughter(0).id(), 3, true));
   outputGenColl->push_back(reco::GenParticle(copysign(1,-1*vb.daughter(1).id())      , phys::Particle::convert(vb.daughter(1).p4()), reco::GenParticle::Point(0.,0.,0.), vb.daughter(1).id(), 3, true));
   outputGenColl->at(initialSize).addDaughter(reco::GenParticleRef(genRefs,initialSize+1));
   outputGenColl->at(initialSize).addDaughter(reco::GenParticleRef(genRefs,initialSize+2));
   
   return outputGenColl;
 }






#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(GenFilterCategory);

 
