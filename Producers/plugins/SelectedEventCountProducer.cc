/** \class SelectedEventCountProducer
 *   An event counter that can store the number of events in the lumi block, counting events that pass a certain path/filter
 *
 *  $Date: 2013/03/15 13:37:31 $
 *  $Revision: 1.3 $
 *  \author R. Bellan - UNITO <riccardo.bellan@cern.ch>
 */

#include <memory>
#include <vector>
#include <algorithm>
#include <iostream>

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/LuminosityBlock.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Utilities/interface/transform.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/Common/interface/MergeableCounter.h"
#include <DataFormats/Common/interface/TriggerResults.h>
#include <FWCore/Common/interface/TriggerNames.h>

#include <boost/foreach.hpp>
#define foreach BOOST_FOREACH


class SelectedEventCountProducer : public edm::one::EDProducer<edm::one::WatchLuminosityBlocks,
							       edm::EndLuminosityBlockProducer> {
  
public:
  explicit SelectedEventCountProducer(const edm::ParameterSet&);
  ~SelectedEventCountProducer(){};

private:
  virtual void produce(edm::Event &, const edm::EventSetup &) override;
  virtual void beginLuminosityBlock(const edm::LuminosityBlock &, const edm::EventSetup&) override;
  virtual void endLuminosityBlock(edm::LuminosityBlock const&, const edm::EventSetup&) {};
  virtual void endLuminosityBlockProduce(edm::LuminosityBlock &, const edm::EventSetup&) override;

 
  // ----------member data ---------------------------
  unsigned int eventsProcessedInLumi_;
  edm::EDGetTokenT<edm::TriggerResults> triggerToken_;
  
  typedef std::pair<std::string, edm::EDGetTokenT<bool> > TokenWihtName;
  typedef std::vector<TokenWihtName> vtoken;
  /// labels of the collections to be merged
  vtoken filterNamesToken_;
 
};

using namespace edm;
using namespace std;

SelectedEventCountProducer::SelectedEventCountProducer(const edm::ParameterSet& iConfig)
  : eventsProcessedInLumi_(0)
  , triggerToken_(consumes<edm::TriggerResults>(edm::InputTag("TriggerResults")))
  , filterNamesToken_( edm::vector_transform(iConfig.template getParameter<std::vector<std::string> >( "names" ), [this](std::string const & tag){return std::make_pair(tag,consumes<bool>(edm::InputTag(tag)));} ) ) 
{
  
  //produces<edm::MergeableCounter, edm::InLumi>();
  produces<edm::MergeableCounter, edm::Transition::EndLuminosityBlock>();
}


void SelectedEventCountProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup){
  
  edm::Handle<edm::TriggerResults> triggerResults;
  edm::TriggerNames triggerNames;
  
  if (iEvent.getByToken(triggerToken_, triggerResults)) {
    triggerNames = iEvent.triggerNames(*triggerResults);
  } else {
    cout << "ERROR: failed to get TriggerResults" << endl;
  }
  
  //for(uint i = 0; i < triggerNames.size(); ++i)
  //  cout<<  triggerNames.triggerName(i)<<endl;
  
  bool passAllRequirements = true;

  foreach(const TokenWihtName &filterNameToken, filterNamesToken_){
    unsigned i = triggerNames.triggerIndex(filterNameToken.first);
    // Search also if a producer put the result inside the event, insteadof into the trigger results
    edm::Handle<bool> filterResult;
    bool foundIntoTheEvent = iEvent.getByToken(filterNameToken.second, filterResult);
    
    if (i == triggerNames.size() && !foundIntoTheEvent){
      cout << "ERROR: SelectedEventCountProducer::isTriggerBit: path does not exist anywhere! " << filterNameToken.first << endl;
      abort();
    }

    if (i != triggerNames.size()) passAllRequirements &=  triggerResults->accept(i);
    else                          passAllRequirements &= *filterResult;

  }
  
  if(passAllRequirements) ++eventsProcessedInLumi_;
}

void SelectedEventCountProducer::beginLuminosityBlock(const edm::LuminosityBlock &theLuminosityBlock, const edm::EventSetup& theSetup) {
  
  eventsProcessedInLumi_ = 0;
}

void SelectedEventCountProducer::endLuminosityBlockProduce(LuminosityBlock & theLuminosityBlock, const EventSetup & theSetup) {

  
  LogTrace("SelectedEventCounting") << "endLumi: adding " << eventsProcessedInLumi_ << " events" << endl;
  auto numEventsPtr = std::make_unique<edm::MergeableCounter>();
  numEventsPtr->value = eventsProcessedInLumi_;
  theLuminosityBlock.put(std::move(numEventsPtr));
}

//define this as a plug-in
DEFINE_FWK_MODULE(SelectedEventCountProducer);
