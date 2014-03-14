#include "FWCore/Framework/interface/MakerMacros.h"
#include "DataFormats/Common/interface/View.h"
#include "CommonTools/UtilAlgos/interface/SingleObjectSelector.h"
#include "CommonTools/UtilAlgos/interface/PtMinSelector.h"
#include "CommonTools/UtilAlgos/interface/EtaRangeSelector.h"
#include "CommonTools/UtilAlgos/interface/AndSelector.h"
#include "AnalysisDataFormats/CMGTools/interface/PFJet.h"

typedef SingleObjectSelector<
  std::vector<cmg::PFJet>,
  AndSelector<
    PtMinSelector,
    EtaRangeSelector
    >,
  std::vector<cmg::PFJet>
  > EtaPtMinCMGPFJetSelector;

DEFINE_FWK_MODULE( EtaPtMinCMGPFJetSelector );
