#ifndef VVXAnalyzer_Producers_Merger_h
#define VVXAnalyzer_Producers_Merger_h
/** \class Merger
 *
 * Merges an arbitrary number of collections 
 * into a single collection.
 * It only uses the collections for which a priority is set, if any.
 * Otherwise it falls back to using all the collections.
 * 
 * Template parameters:
 * - C : collection type
 * - P : policy class that specifies how objects 
 *       in the collection are are cloned
 *
 * \author Luca Lista, INFN (original class "Merger")
 * \author R. Bellan <riccardo.bellan@cern.ch> (modified class with priority)
 * \author A. Mecca <alberto.mecca@cern.ch> (changed name to distingish from the original class, update to edm::global::EDProducer)
 *
 * \version $Revision: 1.2 $
 *
 * $Id: Merger.h,v 1.2 2010/02/20 20:55:21 wmtan Exp $
 *
 */
#include "FWCore/Framework/interface/global/EDProducer.h"
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
class PriorityMerger : public edm::global::EDProducer<> {
 public:
 /// constructor from parameter set
 explicit PriorityMerger( const edm::ParameterSet& );
 /// destructor
 ~PriorityMerger() override;
 
 private:
 /// process an event
 virtual void produce(edm::StreamID, edm::Event&, const edm::EventSetup&) const override;
 
 typedef std::vector<edm::EDGetTokenT<InputCollection> > vtoken;
 /// labels of the collections to be merged
 vtoken srcToken_;
 
 std::vector<int> priority_;
};

  template<typename InputCollection, typename OutputCollection, typename P>
  PriorityMerger<InputCollection, OutputCollection, P>::PriorityMerger( const edm::ParameterSet& par )
  : srcToken_( edm::vector_transform(par.template getParameter<std::vector<edm::InputTag> >( "src" ), [this](edm::InputTag const & tag){return consumes<InputCollection>(tag);} ) ) 
    , priority_( par.template getParameter<std::vector<int> >( "priority" ) ) {
    produces<OutputCollection>();
}

template<typename InputCollection, typename OutputCollection, typename P>
  PriorityMerger<InputCollection, OutputCollection, P>::~PriorityMerger() {
}

template<typename InputCollection, typename OutputCollection, typename P>
void PriorityMerger<InputCollection, OutputCollection, P>::produce(edm::StreamID, edm::Event& evt, const edm::EventSetup&) const {
  auto coll = std::make_unique<OutputCollection>();

  // Fill the collection with labels with priority
  for( typename vtoken::const_iterator s = srcToken_.begin(); s != srcToken_.end(); ++ s ) {
    if(!priority_.at(std::distance(srcToken_.begin(),s))) continue;
    edm::Handle<InputCollection> h;
    evt.getByToken( * s, h );
    for( typename InputCollection::const_iterator c = h->begin(); c != h->end(); ++c ) {
      coll->push_back( P::clone( * c ) );
    }
  }
  
  // if there are not candidates in the priority list, then proceed with the other labels
  if(coll->empty()){
    for( typename vtoken::const_iterator s = srcToken_.begin(); s != srcToken_.end(); ++ s ) {
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
