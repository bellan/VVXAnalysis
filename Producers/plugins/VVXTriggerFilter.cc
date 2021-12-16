/** \class VVXTriggerFilter
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
#include "VVXAnalysis/DataFormats/interface/RegionTypes.h"


class VVXTriggerFilter : public edm::EDFilter {
public:
  /// Constructor
  explicit VVXTriggerFilter(const edm::ParameterSet& config) 
    : filterController_(config, consumesCollector()){
    consumes<edm::TriggerResults>(edm::InputTag("TriggerResults"));
    produces<bool>();
    channel_ = phys::channelType(config.getParameter<std::string>("channelType"));
  }
  
  /// Destructor
  ~VVXTriggerFilter(){};  
  
  virtual void beginJob(){};  
  virtual bool filter(edm::Event&, const edm::EventSetup&);
  virtual void endJob(){};
  
private:
  FilterController filterController_;
  phys::Channel channel_;
};


bool VVXTriggerFilter::filter(edm::Event& event, const edm::EventSetup& setup){  

  Short_t triggerWord(0);

  bool passTrigger = filterController_.passTrigger(phys::UNDEF, event, triggerWord);

  auto output = std::make_unique<bool>(false);

  if(!passTrigger) {event.put(std::move(output));return false;}
  
  bool result = filterController_.passTrigger(channel_, triggerWord); 
  *output = result;
  event.put(std::move(output));
  return result; 
}



#include <FWCore/Framework/interface/MakerMacros.h>
DEFINE_FWK_MODULE(VVXTriggerFilter);
