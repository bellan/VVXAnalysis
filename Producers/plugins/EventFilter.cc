/** \class EventFilter
 *  No description available.
 *
 *  $Date: $
 *  $Revision: $
 *  \author R. Bellan - UNITO <riccardo.bellan@cern.ch>
 */


#include <FWCore/Framework/interface/Frameworkfwd.h>
#include <FWCore/Framework/interface/EDFilter.h>
#include <FWCore/Framework/interface/ESHandle.h>
#include <FWCore/Framework/interface/Event.h>
#include <FWCore/ParameterSet/interface/ParameterSet.h>
#include <FWCore/ServiceRegistry/interface/Service.h>
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/Utilities/interface/EDGetToken.h"
#include <CommonTools/UtilAlgos/interface/TFileService.h>
#include "CommonTools/Utils/interface/StringCutObjectSelector.h"
#include <DataFormats/Candidate/interface/Candidate.h>
#include <DataFormats/PatCandidates/interface/CompositeCandidate.h>

#include <iostream>
#include <vector>
#include <tuple> 

#include <ZZAnalysis/AnalysisStep/interface/bitops.h>
#include "DataFormats/Math/interface/deltaR.h"

#include <boost/foreach.hpp>
#define foreach BOOST_FOREACH

using namespace std;
using namespace edm;


class EventFilter: public edm::EDFilter {

public:
  

  EventFilter(const ParameterSet& pset) 
    : leptonsToken_         (consumes<edm::View<reco::Candidate> >(pset.getParameter<edm::InputTag>("leptons")))
    , tightLeptonSelection_ (pset.getParameter<std::string>("tightLeptonSelection"))
    , minLooseLeptons_      (pset.getParameter<int>("minLooseLeptons"))
    , maxLooseLeptons_      (pset.getParameter<int>("maxLooseLeptons"))
    , minTightLeptons_      (pset.getParameter<int>("minTightLeptons"))
    , maxTightLeptons_      (pset.getParameter<int>("maxTightLeptons"))
    , photonsToken_         (consumes<edm::View<reco::Candidate> >(pset.getParameter<edm::InputTag>("photons")))
    , photonSelection_      (pset.getParameter<std::string>("photonSelection"))
    , minPhotons_           (pset.getParameter<int>("minPhotons"))
    , maxPhotons_           (pset.getParameter<int>("maxPhotons"))
    , jetsAK4Token_         (consumes<edm::View<reco::Candidate> >(pset.getParameter<edm::InputTag>("jetsAK4")))
    , jetAK4Selection_      (pset.getParameter<std::string>("jetAK4Selection"))
    , minAK4s_              (pset.getParameter<int>("minAK4s"))
    , jetsAK8Token_         (consumes<edm::View<reco::Candidate> >(pset.getParameter<edm::InputTag>("jetsAK8")))
    , jetAK8Selection_      (pset.getParameter<std::string>("jetAK8Selection"))
    , minAK8s_              (pset.getParameter<int>("minAK8s"))
    , minAK4orMinAK8_       (pset.getParameter<bool>("minAK4orMinAK8"))
    , METToken_             (consumes<edm::View<reco::Candidate> >(pset.getParameter<edm::InputTag>("MET")))
    , METSelection_         (pset.getParameter<std::string>("METSelection")) {}
 
  
  bool filter(edm::Event & event, const edm::EventSetup& eventSetup);
  
  
private:
  
  std::pair<int,int> countLooseAndTight(const edm::Event & event,
					const edm::EDGetTokenT<edm::View<reco::Candidate> >& token,
					const StringCutObjectSelector<reco::Candidate>& selection);
  
  edm::EDGetTokenT<edm::View<reco::Candidate> > leptonsToken_;
  StringCutObjectSelector<reco::Candidate> tightLeptonSelection_;
  int minLooseLeptons_;
  int maxLooseLeptons_;
  int minTightLeptons_;
  int maxTightLeptons_;
  
  edm::EDGetTokenT<edm::View<reco::Candidate> > photonsToken_;
  StringCutObjectSelector<reco::Candidate> photonSelection_;
  int minPhotons_;
  int maxPhotons_;

  
  edm::EDGetTokenT<edm::View<reco::Candidate> > jetsAK4Token_;
  StringCutObjectSelector<reco::Candidate> jetAK4Selection_;
  int minAK4s_;

  edm::EDGetTokenT<edm::View<reco::Candidate> > jetsAK8Token_;
  StringCutObjectSelector<reco::Candidate> jetAK8Selection_;
  int minAK8s_;

  bool minAK4orMinAK8_;

  edm::EDGetTokenT<edm::View<reco::Candidate> > METToken_;
  StringCutObjectSelector<reco::Candidate> METSelection_;
 
};

std::pair<int,int> EventFilter::countLooseAndTight(const edm::Event & event,
						   const edm::EDGetTokenT<edm::View<reco::Candidate> >& token,
						   const StringCutObjectSelector<reco::Candidate>& selection){
  
  // Get the collection of
  edm::Handle<edm::View<reco::Candidate> > candidates;
  event.getByToken(token, candidates);
  
  int loose  = candidates->size(); int tight = 0;
  
  foreach(const reco::Candidate& cand, *candidates) if(selection(cand)) ++tight;
  
  return std::make_pair(loose,tight);
}



bool EventFilter::filter(Event & event, const EventSetup& eventSetup) { 
  
  std::pair<int,int> nleptons = countLooseAndTight(event, leptonsToken_, tightLeptonSelection_);
  std::pair<int,int> nphotons = countLooseAndTight(event, photonsToken_, photonSelection_);
  std::pair<int,int> njetsAK4 = countLooseAndTight(event, jetsAK4Token_, jetAK4Selection_);
  std::pair<int,int> njetsAK8 = countLooseAndTight(event, jetsAK8Token_, jetAK8Selection_);

  edm::Handle<edm::View<reco::Candidate> >  met;      event.getByToken(METToken_ , met);
  
  bool hasLooseL  = (nleptons.first  >= minLooseLeptons_ && nleptons.first  <= maxLooseLeptons_);
  bool hasTightL  = (nleptons.second >= minTightLeptons_ && nleptons.second <= maxTightLeptons_);
  bool hasPhotons = (nphotons.second >= minPhotons_      && nphotons.second <= maxPhotons_);
  bool hasAK4     = (njetsAK4.second >= minAK4s_);
  bool hasAK8     = (njetsAK8.second >= minAK8s_);
  bool hasMET     = METSelection_(met->front());
  
  return (hasLooseL  && hasTightL &&  
	  hasPhotons &&
	  hasAK4     && hasAK8    &&
	  hasMET);
}


#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(EventFilter);

 
