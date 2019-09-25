/** \class ZZTriggerFilter
 *
 *  No description available.
 *
 *  $Date: 2014/01/28 00:00:00 $
 *  $Revision: 1.1 $
 *  \author R. Bellan (Torino)
 */

#include <FWCore/Framework/interface/Frameworkfwd.h>
#include "FWCore/Framework/interface/ConsumesCollector.h"
#include <FWCore/Framework/interface/EDFilter.h>
#include <FWCore/Framework/interface/Event.h>
#include <FWCore/ParameterSet/interface/ParameterSet.h>
#include <DataFormats/PatCandidates/interface/CompositeCandidate.h>
#include <DataFormats/Candidate/interface/ShallowCloneCandidate.h>

#include "VVXAnalysis/Producers/interface/FilterController.h"
#include <ZZAnalysis/AnalysisStep/interface/FinalStates.h>

class ZZTriggerFilter : public edm::EDFilter {
public:
  /// Constructor
  explicit ZZTriggerFilter(const edm::ParameterSet& config) 
    : filterController_(config, consumesCollector())
    , srcToken_(consumes<edm::View<pat::CompositeCandidate> >(config.getParameter<edm::InputTag>("src")))
  {
    consumes<edm::TriggerResults>(edm::InputTag("TriggerResults"));
    produces<bool>();}
  
  /// Destructor
  ~ZZTriggerFilter(){};  
  
  virtual void beginJob(){};  
  virtual bool filter(edm::Event&, const edm::EventSetup&);
  virtual void endJob(){};
  
private:
  FilterController filterController_;
  edm::EDGetTokenT<edm::View<pat::CompositeCandidate> > srcToken_;
};


bool ZZTriggerFilter::filter(edm::Event& event, const edm::EventSetup& setup){  
  edm::Handle<edm::View<pat::CompositeCandidate> > edmVVs   ; event.getByToken(srcToken_      ,        edmVVs);
  Short_t triggerWord(0);
  bool passTrigger = filterController_.passTrigger(NONE, event, triggerWord);

  auto output = std::make_unique<bool>(false); //Topology

  if(!passTrigger) {event.put(std::move(output));return false;}

  // HACK HERE!! Do not consider cases where there is more than 1 ZZ candidate! 
  // The selection is as the same as in TreePlanter, but here it is not the right place where put the
  // requirement of having 1 and only 1 candidate
  // if(edmVVs->size() != 1) {event.put(output);return false;}
  // const pat::CompositeCandidate& edmVV = edmVVs->front();
  
  // int finalStateZ1 = abs(edmVV.daughter(0)->daughter(0)->pdgId())+abs(edmVV.daughter(0)->daughter(1)->pdgId());
  // int finalStateZ2 = abs(edmVV.daughter(1)->daughter(0)->pdgId())+abs(edmVV.daughter(1)->daughter(1)->pdgId());

  // int rawchannel = finalStateZ1+finalStateZ2;

  Channel effectiveChannel = ZZ;
  // if      (rawchannel == 44) effectiveChannel = EEEE;  // ZZ->4e
  // else if (rawchannel == 48) effectiveChannel = EEMM;  // ZZ->2e2mu
  // else if (rawchannel == 52) effectiveChannel = MMMM;  // ZZ->4mu
  // else {std::cout << "Do not know what to do when setting trigger bit in ZZTriggerFilter. Unknown ZZ id: " << rawchannel << std::endl; abort();}
  // //std::cout <<"channel "<<effectiveChannel<< "  Event: " << event.id().event() << std::endl;

  bool result = filterController_.passTrigger(effectiveChannel, triggerWord); 
  *output = result;
  event.put(std::move(output));
  return result; 
}



#include <FWCore/Framework/interface/MakerMacros.h>
DEFINE_FWK_MODULE(ZZTriggerFilter);
