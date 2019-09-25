/** \class WCandidateFiller
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
#include <FWCore/Framework/interface/ESHandle.h>
#include <CommonTools/Utils/interface/StringCutObjectSelector.h>
#include <CommonTools/Utils/interface/StringObjectFunction.h>
#include <DataFormats/PatCandidates/interface/CompositeCandidate.h>
#include <DataFormats/Candidate/interface/ShallowCloneCandidate.h>
#include <DataFormats/PatCandidates/interface/Muon.h>
#include <DataFormats/PatCandidates/interface/Electron.h>
#include <ZZAnalysis/AnalysisStep/interface/PhotonFwd.h>

#include <ZZAnalysis/AnalysisStep/interface/CutSet.h>
#include <ZZAnalysis/AnalysisStep/interface/LeptonIsoHelper.h>
#include <ZZAnalysis/AnalysisStep/interface/DaughterDataHelpers.h>

#include <string>
#include <Math/VectorUtil.h>

class WCandidateFiller : public edm::EDProducer {
 public:
  /// Constructor
  explicit WCandidateFiller(const edm::ParameterSet&);
    
  /// Destructor
  ~WCandidateFiller(){};  

 private:
  virtual void beginJob(){};  
  virtual void produce(edm::Event&, const edm::EventSetup&);
  virtual void endJob(){};
  
  edm::EDGetTokenT<edm::View<reco::CompositeCandidate> > srcToken_;
};


WCandidateFiller::WCandidateFiller(const edm::ParameterSet& iConfig)
  : srcToken_ (consumes<edm::View<reco::CompositeCandidate> >(iConfig.getParameter<edm::InputTag>("src"))){
  produces<pat::CompositeCandidateCollection>();  
}


void
WCandidateFiller::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{  
  using namespace edm;
  using namespace std;
  using namespace reco;
  
  auto result = std::make_unique<pat::CompositeCandidateCollection>();

  //-- Get jj candidates
  Handle<View<reco::CompositeCandidate> > jjCands;
  iEvent.getByToken(srcToken_, jjCands);


  //--- Fill user info

  const float WmassValue = 80.385;

  float closestWjjMassDiff = 99999.;
  //int bestWjjIdx = -1;

  //--- Loop over LL Candidates
  for(unsigned int i = 0; i < jjCands->size(); ++i) {
    const CompositeCandidate& c = (*jjCands)[i];
    pat::CompositeCandidate myCand(c); 
    
    //--- Find "best Z" (closest to mZ) among those passing the "bestZAmong" selection (2011 PRL logic)
    float diffWmass = fabs(WmassValue - myCand.mass());
    if (diffWmass < closestWjjMassDiff) { // Best among any ll in the collection
      //bestWjjIdx = i;
      closestWjjMassDiff = diffWmass;
    }
   
    result->push_back(myCand);
  }

 
  iEvent.put(std::move(result));
}



#include <FWCore/Framework/interface/MakerMacros.h>
DEFINE_FWK_MODULE(WCandidateFiller);
