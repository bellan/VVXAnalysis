#include "VVXAnalysis/TreeAnalysis/interface/AnalysisFactory.h"
#include "VVXAnalysis/TreeAnalysis/interface/WlllnuAnalyzer.h"
#include "VVXAnalysis/TreeAnalysis/interface/VVXAnalyzer.h"
#include "VVXAnalysis/TreeAnalysis/interface/ZZWAnalyzer.h"
#include "VVXAnalysis/TreeAnalysis/interface/WZAnalyzer.h"
#include "VVXAnalysis/TreeAnalysis/interface/ZZSAnalyzer.h"
#include "VVXAnalysis/TreeAnalysis/interface/ZZjAnalyzer.h"
#include "VVXAnalysis/TreeAnalysis/interface/ZZjGenAnalyzer.h"
#include "VVXAnalysis/TreeAnalysis/interface/ZZWSRDefinition.h"
#include "VVXAnalysis/TreeAnalysis/interface/WWosAnalyzer.h"
#include "VVXAnalysis/TreeAnalysis/interface/ZZMCAnalyzer.h"
#include "VVXAnalysis/TreeAnalysis/interface/ZZRecoAnalyzer.h"
#include "VVXAnalysis/TreeAnalysis/interface/VBSMCAnalyzer.h"
#include "VVXAnalysis/TreeAnalysis/interface/VBSRecoAnalyzer.h"
#include "VVXAnalysis/TreeAnalysis/interface/ZVAnalyzer.h"
#include "VVXAnalysis/TreeAnalysis/interface/FakeRateAnalyzer.h"
#include "VVXAnalysis/TreeAnalysis/interface/WZZAnalyzer.h"
#include "VVXAnalysis/TreeAnalysis/interface/VZZAnalyzer.h"
#include "VVXAnalysis/TreeAnalysis/interface/VZZaQGCAnalyzer.h"
#include "VVXAnalysis/TreeAnalysis/interface/VVXnocutsAnalyzer.h"
#include "VVXAnalysis/TreeAnalysis/interface/DBGenerator.h"
#include "VVXAnalysis/TreeAnalysis/interface/ZZjjAnalyzer.h"
#include "VVXAnalysis/TreeAnalysis/interface/VVGammaAnalyzer.h"
#include "VVXAnalysis/TreeAnalysis/interface/VZGammaAnalyzer.h"
#include "VVXAnalysis/TreeAnalysis/interface/VZGAnalyzer.h"

AnalysisFactory::AnalysisFactory(){
  Register("WlllnuAnalyzer"   , &RegistrableAnalysis<WlllnuAnalyzer>::create    );
  Register("VVXAnalyzer"      , &RegistrableAnalysis<VVXAnalyzer>::create       );
  Register("ZZWAnalyzer"      , &RegistrableAnalysis<ZZWAnalyzer>::create       );
  Register("WZZAnalyzer"      , &RegistrableAnalysis<WZZAnalyzer>::create       );
  Register("WZAnalyzer"       , &RegistrableAnalysis<WZAnalyzer>::create        );
  Register("ZZSAnalyzer"      , &RegistrableAnalysis<ZZSAnalyzer>::create       );
  Register("ZZjAnalyzer"      , &RegistrableAnalysis<ZZjAnalyzer>::create       );
  Register("ZZjGenAnalyzer"   , &RegistrableAnalysis<ZZjGenAnalyzer>::create    );
  Register("ZZWSRDefinition"  , &RegistrableAnalysis<ZZWSRDefinition>::create   );
  Register("ZZMCAnalyzer"     , &RegistrableAnalysis<ZZMCAnalyzer>::create      );
  Register("ZZRecoAnalyzer"   , &RegistrableAnalysis<ZZRecoAnalyzer>::create    );
  Register("VBSMCAnalyzer"    , &RegistrableAnalysis<VBSMCAnalyzer>::create     );
  Register("VBSRecoAnalyzer"  , &RegistrableAnalysis<VBSRecoAnalyzer>::create   );
  Register("FakeRateAnalyzer" , &RegistrableAnalysis<FakeRateAnalyzer>::create  );
  Register("ZVAnalyzer"       , &RegistrableAnalysis<ZVAnalyzer>::create        );
  Register("WWosAnalyzer"     , &RegistrableAnalysis<WWosAnalyzer>::create      );
  Register("VZZAnalyzer"      , &RegistrableAnalysis<VZZAnalyzer>::create       );
  Register("VZZaQGCAnalyzer"  , &RegistrableAnalysis<VZZaQGCAnalyzer>::create   );
  Register("VVXnocutsAnalyzer", &RegistrableAnalysis<VVXnocutsAnalyzer>::create );
  Register("DBGenerator"      , &RegistrableAnalysis<DBGenerator>::create       );
  Register("ZZjjAnalyzer"     , &RegistrableAnalysis<ZZjjAnalyzer>::create      );
  Register("VVGammaAnalyzer"  , &RegistrableAnalysis<VVGammaAnalyzer>::create      );
  Register("VZGammaAnalyzer"  , &RegistrableAnalysis<VZGammaAnalyzer>::create      );
  Register("VZGAnalyzer"  , &RegistrableAnalysis<VZGAnalyzer>::create      );
}

void AnalysisFactory::Register(const std::string &analysisName, CreateAnFn pfnCreate)
{
    m_FactoryMap[analysisName] = pfnCreate;
}

//EventAnalyzer *AnalysisFactory::createAnalysis(const std::string &analysisName, const std::string& region, const std::string& filename, const double& lumi, const double& externalXSection, bool doBasicPlots)
EventAnalyzer *AnalysisFactory::createAnalysis(const AnalysisConfiguration &analysisConfiguration)
{
  FactoryMap::iterator it = m_FactoryMap.find(analysisConfiguration.getParameter<std::string>("analysisName"));
  if( it != m_FactoryMap.end() ){
    return it->second(analysisConfiguration);
  }
  return NULL;
}

