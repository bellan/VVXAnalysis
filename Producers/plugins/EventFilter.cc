/** \class EventFilter
 *  No description available.
 *
 *  $Date: $
 *  $Revision: $
 *  \author R. Bellan - UNITO <riccardo.bellan@cern.ch>
 */


#include <FWCore/Framework/interface/Frameworkfwd.h>
#include <FWCore/Framework/interface/stream/EDFilter.h>
#include <FWCore/Framework/interface/Event.h>
#include <FWCore/ParameterSet/interface/ParameterSet.h>
#include "FWCore/Utilities/interface/EDGetToken.h"
#include "CommonTools/Utils/interface/StringCutObjectSelector.h"
#include "DataFormats/PatCandidates/interface/Photon.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/MET.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Electron.h"

#include <iostream>
#include <vector>
#include <tuple> 

#include <ZZAnalysis/AnalysisStep/interface/bitops.h>
#include "DataFormats/Math/interface/deltaR.h"

#include <boost/foreach.hpp>
#define foreach BOOST_FOREACH

using namespace std;
using namespace edm;


class EventFilter: public edm::stream::EDFilter<> {

public:
  

  EventFilter(const ParameterSet& pset) 
    : muonsToken_             (consumes<pat::MuonCollection>(pset.getParameter<edm::InputTag>("muons")))
    , tightMuonSelection_     (pset.getParameter<std::string>("tightMuonSelection"))
    , electronsToken_         (consumes<pat::ElectronCollection>(pset.getParameter<edm::InputTag>("electrons")))
    , tightElectronSelection_ (pset.getParameter<std::string>("tightElectronSelection"))
    , minLooseLeptons_      (pset.getParameter<int>("minLooseLeptons"))
    , maxLooseLeptons_      (pset.getParameter<int>("maxLooseLeptons"))
    , minTightLeptons_      (pset.getParameter<int>("minTightLeptons"))
    , maxTightLeptons_      (pset.getParameter<int>("maxTightLeptons"))
    , photonsToken_         (consumes<pat::PhotonCollection>(pset.getParameter<edm::InputTag>("photons")))
    , photonSelection_      (pset.getParameter<std::string>("photonSelection"))
    , minPhotons_           (pset.getParameter<int>("minPhotons"))
    , maxPhotons_           (pset.getParameter<int>("maxPhotons"))
    , jetsAK4Token_         (consumes<pat::JetCollection>(pset.getParameter<edm::InputTag>("jetsAK4")))
    , jetAK4Selection_      (pset.getParameter<std::string>("jetAK4Selection"))
    , minAK4s_              (pset.getParameter<int>("minAK4s"))
    , jetsAK8Token_         (consumes<pat::JetCollection>(pset.getParameter<edm::InputTag>("jetsAK8")))
    , jetAK8Selection_      (pset.getParameter<std::string>("jetAK8Selection"))
    , minAK8s_              (pset.getParameter<int>("minAK8s"))
    , minAK4orMinAK8_       (pset.getParameter<bool>("minAK4orMinAK8"))
    , METToken_             (consumes<pat::METCollection>(pset.getParameter<edm::InputTag>("MET")))
    , METSelection_         (pset.getParameter<std::string>("METSelection")) {


    
}
 
  
  bool filter(edm::Event & event, const edm::EventSetup& eventSetup);
  
  
private:
  
  template <typename T>
  std::pair<int,int> countLooseAndTight(const edm::Event & event,
					const edm::EDGetTokenT<std::vector<T> >& token,
					const StringCutObjectSelector<T>& selection);
  
  edm::EDGetTokenT<pat::MuonCollection> muonsToken_;
  StringCutObjectSelector<pat::Muon> tightMuonSelection_;
  edm::EDGetTokenT<pat::ElectronCollection> electronsToken_;
  StringCutObjectSelector<pat::Electron> tightElectronSelection_;

  int minLooseLeptons_;
  int maxLooseLeptons_;
  int minTightLeptons_;
  int maxTightLeptons_;
  
  edm::EDGetTokenT<pat::PhotonCollection> photonsToken_;
  StringCutObjectSelector<pat::Photon> photonSelection_;
  int minPhotons_;
  int maxPhotons_;

  
  edm::EDGetTokenT<pat::JetCollection> jetsAK4Token_;
  StringCutObjectSelector<pat::Jet> jetAK4Selection_;
  int minAK4s_;

  edm::EDGetTokenT<pat::JetCollection> jetsAK8Token_;
  StringCutObjectSelector<pat::Jet> jetAK8Selection_;
  int minAK8s_;

  bool minAK4orMinAK8_;

  edm::EDGetTokenT<pat::METCollection> METToken_;
  StringCutObjectSelector<pat::MET> METSelection_;
 
};

template <typename T>
std::pair<int,int> EventFilter::countLooseAndTight(const edm::Event & event,
						   const edm::EDGetTokenT<std::vector<T> >& token,
						   const StringCutObjectSelector<T>& selection){
  
  // Get the collection of
  edm::Handle<std::vector<T> > candidates;
  event.getByToken(token, candidates);
  
  int loose  = candidates->size(); int tight = 0;
  
  foreach(const T& cand, *candidates) if(selection(cand)) ++tight;
  
  return std::make_pair(loose,tight);
}



bool EventFilter::filter(Event & event, const EventSetup& eventSetup) { 
  
  std::pair<int,int> nmuons     = countLooseAndTight(event, muonsToken_, tightMuonSelection_);
  std::pair<int,int> nelectrons = countLooseAndTight(event, electronsToken_, tightElectronSelection_);
  std::pair<int,int> nphotons = countLooseAndTight(event, photonsToken_, photonSelection_);
  std::pair<int,int> njetsAK4 = countLooseAndTight(event, jetsAK4Token_, jetAK4Selection_);
  std::pair<int,int> njetsAK8 = countLooseAndTight(event, jetsAK8Token_, jetAK8Selection_);

  std::pair<int,int> nleptons = std::make_pair(nmuons.first+nelectrons.first,nmuons.second+nelectrons.second);

  edm::Handle<pat::METCollection>  met;      event.getByToken(METToken_ , met);
  
  bool hasLooseL  = (nleptons.first  >= minLooseLeptons_ && nleptons.first  <= maxLooseLeptons_);
  bool hasTightL  = (nleptons.second >= minTightLeptons_ && nleptons.second <= maxTightLeptons_);

  bool hasPhotons = (nphotons.second >= minPhotons_      && nphotons.second <= maxPhotons_);
  bool hasAK4     = (njetsAK4.second >= minAK4s_);
  bool hasAK8     = (njetsAK8.second >= minAK8s_);
  
  bool hasAKX     = minAK4orMinAK8_ ? (hasAK4 || hasAK8) : (hasAK4 && hasAK8);
  bool hasMET     = METSelection_(met->front());
  
  return (hasLooseL  && hasTightL &&  
	  hasPhotons &&
	  hasAKX     &&
	  hasMET);
}


#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(EventFilter);

 
