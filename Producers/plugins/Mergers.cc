#include "FWCore/Framework/interface/MakerMacros.h"
#include "VVXAnalysis/Producers/plugins/Merger.h"
#include "DataFormats/Candidate/interface/Candidate.h"

#include "DataFormats/Candidate/interface/Candidate.h"
#include <DataFormats/PatCandidates/interface/CompositeCandidate.h>

// A merger which produces a collection of pat::CompositeCandidate
// (instead of a collection of reco::CompositeCandidate)
typedef ZZAnalysis::Merger<pat::CompositeCandidateCollection, pat::CompositeCandidateCollection> PATCompositeCandidateMergerWithPriority;


DEFINE_FWK_MODULE( PATCompositeCandidateMergerWithPriority );
