/** \class WZCandidateFiller
 *
 *  No description available.
 *
 *  $Date: 2014/01/28 00:00:00 $
 *  $Revision: 1.1 $
 *  \author R. Bellan (Torino)
 */

#include <FWCore/Framework/interface/Frameworkfwd.h>
#include <FWCore/Framework/interface/EDProducer.h>
#include <FWCore/Framework/interface/Event.h>
#include <FWCore/ParameterSet/interface/ParameterSet.h>
#include <DataFormats/PatCandidates/interface/CompositeCandidate.h>
#include "CommonTools/Utils/interface/StringCutObjectSelector.h"

#include "DataFormats/PatCandidates/interface/MET.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Electron.h"

#include "VVXAnalysis/DataFormats/interface/Particle.h"

#include <TLorentzVector.h>

#include <boost/foreach.hpp>
#define foreach BOOST_FOREACH


class WZCandidateFiller : public edm::EDProducer {

public:
  /// Constructor
  explicit WZCandidateFiller(const edm::ParameterSet& config)
    : theZLToken (consumes<edm::View<pat::CompositeCandidate> >(config.getParameter<edm::InputTag>("ZL")))
    , theMETToken(consumes<pat::METCollection>(config.getParameter<edm::InputTag>("MET")))
  {
    produces<pat::CompositeCandidateCollection>(); 

  }
  
  /// Destructor
  ~WZCandidateFiller(){};  
  
  virtual void beginJob(){};  
  virtual void produce(edm::Event&, const edm::EventSetup&);
  virtual void endJob(){};
  
private:

  typedef pat::CompositeCandidate ZLCompositeCandidate;
  typedef std::vector<ZLCompositeCandidate> ZLCompositeCandidates;

  edm::EDGetTokenT<edm::View<pat::CompositeCandidate> > theZLToken;
  edm::EDGetTokenT<pat::METCollection> theMETToken;

  bool isMC_;
 

};


void WZCandidateFiller::produce(edm::Event& event, const edm::EventSetup& setup){  

  
  edm::Handle<edm::View<pat::CompositeCandidate> > ZLs; event.getByToken(theZLToken ,  ZLs);
  edm::Handle<pat::METCollection> MET                 ;      event.getByToken(theMETToken    , MET);
  
   std::auto_ptr<pat::CompositeCandidateCollection> result(new pat::CompositeCandidateCollection);

  
  if(MET->at(0).pt() < 30) return;
  
  foreach(const ZLCompositeCandidate ZL, *ZLs){
    
    const pat::CompositeCandidate* edmV0  = dynamic_cast<const pat::CompositeCandidate*>(ZL.daughter(0)->masterClone().get());

    TLorentzVector Lepp4;
    bool isGoodLep = false;

    if     (abs(ZL.daughter(1)->pdgId()) == 11){ 
      Lepp4 = phys::Particle::convert((dynamic_cast<const pat::Electron*>(ZL.daughter(1)->masterClone().get()))->p4());
      isGoodLep = (dynamic_cast<const pat::Electron*>(ZL.daughter(1)->masterClone().get()))->userFloat("isGood");
    }
    else if(abs(ZL.daughter(1)->pdgId()) == 13) {
      Lepp4 = phys::Particle::convert((dynamic_cast<const pat::Muon*>(ZL.daughter(1)->masterClone().get()))->p4());
    isGoodLep = (dynamic_cast<const pat::Muon*>(ZL.daughter(1)->masterClone().get()))->userFloat("isGood");}
    else{
      edm::LogError("WZCandidateFiller") << "Do not know what to cast in fillZLCandidates, LEP part" << std::endl;
      abort();
    }


    
    if(edmV0->mass() > 60 && edmV0->mass() < 120 && isGoodLep){
      
      TLorentzVector Zp4   = phys::Particle::convert(edmV0->p4());
      TLorentzVector W     = Zp4 + Lepp4;
      
      if(W.Mt() > 30 && W.Mt() < 500) result->push_back(ZL);
    }
    
  }

  event.put(result);

  
}


#include <FWCore/Framework/interface/MakerMacros.h>
DEFINE_FWK_MODULE(WZCandidateFiller);
