/** \class ZLFilter
 *
 *  No description available.
 *
 *  $Date: 2014/01/28 00:00:00 $
 *  $Revision: 1.1 $
 *  \author R. Bellan (Torino)
 */

#include <FWCore/Framework/interface/Frameworkfwd.h>
#include <FWCore/Framework/interface/EDFilter.h>
#include <FWCore/Framework/interface/Event.h>
#include <FWCore/ParameterSet/interface/ParameterSet.h>
#include <DataFormats/PatCandidates/interface/CompositeCandidate.h>


class ZLFilter : public edm::EDFilter {
public:
  /// Constructor
  explicit ZLFilter(const edm::ParameterSet& config) 
    : theZLLLabel(config.getParameter<edm::InputTag>("ZLL"))
    , theZLLabel(config.getParameter<edm::InputTag>("ZL"))
  {}
  
  /// Destructor
  ~ZLFilter(){};  
  
  virtual void beginJob(){};  
  virtual bool filter(edm::Event&, const edm::EventSetup&);
  virtual void endJob(){};
  
private:
  edm::InputTag theZLLLabel;
  edm::InputTag theZLLabel;
};


bool ZLFilter::filter(edm::Event& event, const edm::EventSetup& setup){  


  edm::Handle<edm::View<pat::CompositeCandidate> > ZLL  ; event.getByLabel(theZLLLabel     ,        ZLL);
  edm::Handle<edm::View<pat::CompositeCandidate> > ZL   ; event.getByLabel(theZLLabel      ,        ZL);


  // at least a ZZ or a ZLL candidate in the event --> want to keep it
  if(ZLL->size() > 0) return true;
  // if there is only one ZL candidate --> want to keep it if it passes trigger requirements
  if(ZL->size() == 1) 
    if((ZL->at(0).daughter(0)->daughter(0)->pt() > 20 && ZL->at(0).daughter(0)->daughter(1)->pt() > 10) || (ZL->at(0).daughter(0)->daughter(0)->pt() > 10 && ZL->at(0).daughter(0)->daughter(1)->pt() > 20))
      return true;

  return false;


}



#include <FWCore/Framework/interface/MakerMacros.h>
DEFINE_FWK_MODULE(ZLFilter);
