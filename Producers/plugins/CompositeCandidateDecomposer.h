#ifndef VVXAnalysis_Producers_plugins_CompositeCandidateDecomposer_h
#define VVXAnalysis_Producers_plugins_CompositeCandidateDecomposer_h

#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "DataFormats/PatCandidates/interface/CompositeCandidate.h"


#include <boost/foreach.hpp>
#define foreach BOOST_FOREACH


namespace pat{
  template<class PATObjType>
    class  CompositeCandidateDecomposer: public edm::EDProducer {
  public:
    explicit CompositeCandidateDecomposer(const edm::ParameterSet & iConfig);
  virtual ~CompositeCandidateDecomposer() { }
  
  virtual void produce(edm::Event & iEvent, const edm::EventSetup & iSetup);
  
  private:
  /// Label for input collection
  edm::EDGetTokenT<edm::View<pat::CompositeCandidate> > srcToken_;
  int splitLevel_;
  };
  
  
  template<class PATObjType>
    pat::CompositeCandidateDecomposer<PATObjType>::CompositeCandidateDecomposer(const edm::ParameterSet & iConfig)
    : srcToken_(consumes<edm::View<pat::CompositeCandidate> >(iConfig.getParameter<edm::InputTag>("src")))
    , splitLevel_    (iConfig.getParameter<int>("SplitLevel")){
    
    produces<std::vector<PATObjType> >();
    
  }

  template <class PATObjType>
    void  pat::CompositeCandidateDecomposer<PATObjType>::produce(edm::Event & event, const edm::EventSetup & iSetup) {
    
    
    edm::Handle<edm::View<pat::CompositeCandidate> > CC   ; event.getByToken(srcToken_, CC);
    auto out = std::make_unique<std::vector<PATObjType> >();
      

    if(splitLevel_ != 0 && splitLevel_ != 1 && splitLevel_ != 2)
      edm::LogError("CompositeCandidateDecomposer") << "CompositeCandidateDecomposer: Unknown split level: " << splitLevel_;
    

    // FIXME: remove the split level and do it auto

    foreach(const pat::CompositeCandidate &cc, *CC){
      
      for(unsigned int i = 0; i < cc.numberOfDaughters(); ++i){
	
	if(splitLevel_ == 0){
	  if(cc.daughter(i)->hasMasterClone()){
	    const PATObjType* d = dynamic_cast<const PATObjType*>(cc.daughter(i)->masterClone().get());
	    if(d) out->push_back(*d);
	  }
	}
	
	else if(splitLevel_ == 1){ 
	  const CompositeCandidate* d = dynamic_cast<const CompositeCandidate*>(cc.daughter(i)->masterClone().get());	  
	  if(!d)  edm::LogError("CompositeCandidateDecomposer") << "Split level is 1, but the daughter is not a composite candidate";

	  for(unsigned int j = 0; j < d->numberOfDaughters(); ++j){
	    if(d->daughter(j)->hasMasterClone()){
	      const PATObjType* gd = dynamic_cast<const PATObjType*>(d->daughter(j)->masterClone().get());
	      if(gd) out->push_back(*gd);
	    }
	  }
	}
      }

      if(splitLevel_ == 2){
	const CompositeCandidate* d = dynamic_cast<const CompositeCandidate*>(cc.daughter(0)->daughter(0)->masterClone().get());	  
	if(!d)  edm::LogError("CompositeCandidateDecomposer") << "Split level is 2, but the daughter is not a composite candidate";
	
	for(unsigned int j = 0; j < d->numberOfDaughters(); ++j){
	  if(d->daughter(j)->hasMasterClone()){
	    const PATObjType* gd = dynamic_cast<const PATObjType*>(d->daughter(j)->masterClone().get());
	    if(gd) out->push_back(*gd);
	  }
	}
	if(cc.daughter(0)->daughter(1)->hasMasterClone()){
	  const PATObjType* gd = dynamic_cast<const PATObjType*>(cc.daughter(0)->daughter(1)->masterClone().get());	  
	  if(gd) out->push_back(*gd);
	}
	
      }      
    }    
    event.put(std::move(out));
  }
}
#endif
