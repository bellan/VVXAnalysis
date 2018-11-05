/** \class EventTagger
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


class EventTagger: public edm::EDFilter {

public:
  

  EventTagger(const ParameterSet& pset) 
    : sel_            (pset.getParameter<int>("Topology"))
    , leptonsToken_   (consumes<edm::View<reco::Candidate> >(pset.getParameter<edm::InputTag>("src")))
    , tightSelection_ (pset.getParameter<std::string>("TightSelection")){
    produces<int>();
  }
 

  bool filter(edm::Event & event, const edm::EventSetup& eventSetup);
    
  
private:
  int sel_;
  edm::EDGetTokenT<edm::View<reco::Candidate> > leptonsToken_;
  StringCutObjectSelector<reco::Candidate> tightSelection_;
  
};


bool EventTagger::filter(Event & event, const EventSetup& eventSetup) { 
  

  // Get the collection of gen particles
  edm::Handle<edm::View<reco::Candidate> > leptons;
  event.getByToken(leptonsToken_, leptons);
  
  int allLep  = leptons->size(); int tightLep  = 0; int looseNotTightLep = 0;
  int sumid = 0;
  
  bitset<5>  topology;


  foreach(const reco::Candidate& lep, *leptons){
    
    if (tightSelection_(lep)) topology.set(5-tightLep++);
    else                      topology.set(looseNotTightLep++);   
    sumid += lep.pdgId();
  }
  
 
  
  if (tightLep < 1 || allLep < 2                                // At least 1T+1LnotT in the event
      || allLep > 4                                             // At most 4L in the event
      || looseNotTightLep > 3                                   // At most 3LnotT in the event
      || (allLep == 4 && sumid != 0)                            // 4l -> only ZZ like events
      || (allLep == 3 && abs(sumid) != 11 && abs(sumid) != 13)) // 3l -> only WZ like events
    topology.reset();
  

  std::auto_ptr<int> output(new int(topology.to_ulong())); //Topology
  event.put(output); 

  if (sel_ >= 0 && topology.any()) {

    if( ((topology.to_ulong() ^ sel_) & sel_) == 0)   return true;

    else return false;

  }

  else
    return true;

}


#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(EventTagger);

 
