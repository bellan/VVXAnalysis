#include "VVXAnalysis/Producers/plugins/CompositeCandidateDecomposer.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/CompositeCandidate.h"

namespace pat {
  typedef pat::CompositeCandidateDecomposer<pat::Electron> PATElectronsFromCompositeCandidates;
  typedef pat::CompositeCandidateDecomposer<pat::Muon>     PATMuonsFromCompositeCandidates;
  typedef pat::CompositeCandidateDecomposer<pat::Jet>      PATJetsFromCompositeCandidates;
  typedef pat::CompositeCandidateDecomposer<pat::Jet>      PATCompositeCandidatesFromCompositeCandidates;
}
using namespace pat;
DEFINE_FWK_MODULE(PATElectronsFromCompositeCandidates);
DEFINE_FWK_MODULE(PATMuonsFromCompositeCandidates);
DEFINE_FWK_MODULE(PATJetsFromCompositeCandidates);
DEFINE_FWK_MODULE(PATCompositeCandidatesFromCompositeCandidates);

