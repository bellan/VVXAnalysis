#ifndef VVXAnalyzer_Producers_Merger_h
#define VVXAnalyzer_Producers_Merger_h
/** \class Merger
 *
 * Merges an arbitrary number of collections 
 * into a single collection.
 * 
 * Template parameters:
 * - C : collection type
 * - P : policy class that specifies how objects 
 *       in the collection are are cloned
 *
 * \author Luca Lista, INFN
 *
 * \version $Revision: 1.2 $
 *
 * $Id: Merger.h,v 1.2 2010/02/20 20:55:21 wmtan Exp $
 *
 */
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/transform.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "DataFormats/Common/interface/CloneTrait.h"
#include <vector>
#include <iterator>

namespace ZZAnalysis {
template<typename InputCollection, 
  typename OutputCollection = InputCollection,
  typename P = typename edm::clonehelper::CloneTrait<InputCollection>::type>
  class Merger : public edm::EDProducer {
 public:
 /// constructor from parameter set
 explicit Merger( const edm::ParameterSet& );
 /// destructor
 ~Merger();
 
 private:
 /// process an event
 virtual void produce( edm::Event&, const edm::EventSetup&) override;
 
 typedef std::vector<edm::EDGetTokenT<InputCollection> > vtoken;
 /// labels of the collections to be merged
 vtoken srcToken_;
 
 std::vector<int> priority_;
};

  template<typename InputCollection, typename OutputCollection, typename P>
  Merger<InputCollection, OutputCollection, P>::Merger( const edm::ParameterSet& par )
  : srcToken_( edm::vector_transform(par.template getParameter<std::vector<edm::InputTag> >( "src" ), [this](edm::InputTag const & tag){return consumes<InputCollection>(tag);} ) ) 
    , priority_( par.template getParameter<std::vector<int> >( "priority" ) ) {
    produces<OutputCollection>();
}

template<typename InputCollection, typename OutputCollection, typename P>
  Merger<InputCollection, OutputCollection, P>::~Merger() {
}

template<typename InputCollection, typename OutputCollection, typename P>
  void Merger<InputCollection, OutputCollection, P>::produce( edm::Event& evt, const edm::EventSetup&) {
  auto coll = std::make_unique<OutputCollection>();

  // Fill the collection with labels with priority
  for( typename vtoken::iterator s = srcToken_.begin(); s != srcToken_.end(); ++ s ) {
    if(!priority_.at(std::distance(srcToken_.begin(),s))) continue;
    edm::Handle<InputCollection> h;
    evt.getByToken( * s, h );
    for( typename InputCollection::const_iterator c = h->begin(); c != h->end(); ++c ) {
      coll->push_back( P::clone( * c ) );
    }
  }
  
  // if there are not candidates in the priority list, then proceed with the other labels
  if(coll->empty()){
    for( typename vtoken::iterator s = srcToken_.begin(); s != srcToken_.end(); ++ s ) {
      if(priority_.at(std::distance(srcToken_.begin(),s))) continue;
      edm::Handle<InputCollection> h;
      evt.getByToken( * s, h );
      for( typename InputCollection::const_iterator c = h->begin(); c != h->end(); ++c ) {
	coll->push_back( P::clone( * c ) );
      }
    }
  }
  
  evt.put( std::move(coll));
}
}

#endif
