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
#include "CommonTools/Utils/interface/StringCutObjectSelector.h"
#include <ZZAnalysis/AnalysisStep/interface/bitops.h>

class ZLFilter : public edm::EDFilter {
public:
  /// Constructor
  explicit ZLFilter(const edm::ParameterSet& config) 
    : theZLLToken(consumes<edm::View<pat::CompositeCandidate> >(config.getParameter<edm::InputTag>("ZLL")))
    , theZLToken(consumes<edm::View<pat::CompositeCandidate> >(config.getParameter<edm::InputTag>("ZL")))
    , selectionZL_(config.getParameter<std::string>("ZLSelection"))
    , isMC_ (config.getUntrackedParameter<bool>("isMC",false))
  {

    theGenCategoryToken = consumes<int> (config.getUntrackedParameter<edm::InputTag>("GenCategory" , edm::InputTag("genCategory")));

  }
  
  /// Destructor
  ~ZLFilter(){};  
  
  virtual void beginJob(){};  
  virtual bool filter(edm::Event&, const edm::EventSetup&);
  virtual void endJob(){};
  
private:

  edm::EDGetTokenT<edm::View<pat::CompositeCandidate> > theZLLToken;
  edm::EDGetTokenT<edm::View<pat::CompositeCandidate> > theZLToken;
  edm::EDGetTokenT<int> theGenCategoryToken;

  StringCutObjectSelector<pat::CompositeCandidate> selectionZL_;

  bool isMC_;

};


bool ZLFilter::filter(edm::Event& event, const edm::EventSetup& setup){  

  edm::Handle<int> genCategory;  event.getByToken(theGenCategoryToken, genCategory);

  edm::Handle<edm::View<pat::CompositeCandidate> > ZLL  ; event.getByToken(theZLLToken, ZLL);
  edm::Handle<edm::View<pat::CompositeCandidate> > ZL   ; event.getByToken(theZLToken ,  ZL);

  if(ZLL->size() > 0) return true;
  // if there is only one ZL candidate --> want to keep it if it passes trigger requirements
  if(ZL->size() == 1 && selectionZL_(ZL->at(0))) return true; 
  if(isMC_ && test_bit_16(*genCategory,0)) return true;
    
  return false;
}



#include <FWCore/Framework/interface/MakerMacros.h>
DEFINE_FWK_MODULE(ZLFilter);
