#include "FWCore/Framework/interface/MakerMacros.h"
#include "VVXAnalysis/Producers/interface/PriorityMerger.h"
#include "DataFormats/Candidate/interface/Candidate.h"

#include "DataFormats/Candidate/interface/Candidate.h"
#include <DataFormats/PatCandidates/interface/CompositeCandidate.h>

// A merger which produces a collection of pat::CompositeCandidate
// (instead of a collection of reco::CompositeCandidate)
typedef ZZAnalysis::PriorityMerger<pat::CompositeCandidateCollection, pat::CompositeCandidateCollection> PATCompositeCandidateMergerWithPriority;
DEFINE_FWK_MODULE( PATCompositeCandidateMergerWithPriority );

#include "CommonTools/UtilAlgos/interface/Merger.h"
#include <DataFormats/PatCandidates/interface/Muon.h>
typedef Merger<pat::MuonCollection, pat::MuonCollection> PATMuonCollectionMerger;
DEFINE_FWK_MODULE( PATMuonCollectionMerger );

#include <DataFormats/PatCandidates/interface/Electron.h>
typedef Merger<pat::ElectronCollection, pat::ElectronCollection> PATElectronCollectionMerger;
DEFINE_FWK_MODULE( PATElectronCollectionMerger );

