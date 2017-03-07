#include "VVXAnalysis/TreeAnalysis/interface/AnalysisFactory.h"
#include "VVXAnalysis/TreeAnalysis/interface/VVXAnalyzer.h"
#include "VVXAnalysis/TreeAnalysis/interface/ZZWAnalyzer.h"
#include "VVXAnalysis/TreeAnalysis/interface/ZZSAnalyzer.h"
#include "VVXAnalysis/TreeAnalysis/interface/ZZjAnalyzer.h"
#include "VVXAnalysis/TreeAnalysis/interface/ZZjGenAnalyzer.h"
#include "VVXAnalysis/TreeAnalysis/interface/ZZWSRDefinition.h"
//#include "VVXAnalysis/TreeAnalysis/interface/JetAnalyzer.h"
//#include "VVXAnalysis/TreeAnalysis/interface/ZZSDataAnalyzer.h"
#include "VVXAnalysis/TreeAnalysis/interface/ZZMCAnalyzer.h"
#include "VVXAnalysis/TreeAnalysis/interface/ZZRecoAnalyzer.h"
#include "VVXAnalysis/TreeAnalysis/interface/VBSMCAnalyzer.h"
#include "VVXAnalysis/TreeAnalysis/interface/VBSRecoAnalyzer.h"
#include "VVXAnalysis/TreeAnalysis/interface/FakeRateAnalyzer.h"
AnalysisFactory::AnalysisFactory(){
  Register("VVXAnalyzer", &RegistrableAnalysis<VVXAnalyzer>::create);
  Register("ZZWAnalyzer", &RegistrableAnalysis<ZZWAnalyzer>::create);
  Register("ZZSAnalyzer", &RegistrableAnalysis<ZZSAnalyzer>::create);
  Register("ZZjAnalyzer", &RegistrableAnalysis<ZZjAnalyzer>::create);
  Register("ZZjGenAnalyzer", &RegistrableAnalysis<ZZjGenAnalyzer>::create);
  Register("ZZWSRDefinition", &RegistrableAnalysis<ZZWSRDefinition>::create);
  //Register("JetAnalyzer", &RegistrableAnalysis<JetAnalyzer>::create);
  //Register("ZZSDataAnalyzer", &RegistrableAnalysis<ZZSDataAnalyzer>::create); 
  Register("ZZMCAnalyzer", &RegistrableAnalysis<ZZMCAnalyzer>::create);
  Register("ZZRecoAnalyzer", &RegistrableAnalysis<ZZRecoAnalyzer>::create);
  Register("VBSMCAnalyzer", &RegistrableAnalysis<VBSMCAnalyzer>::create);
  Register("VBSRecoAnalyzer", &RegistrableAnalysis<VBSRecoAnalyzer>::create);
  Register("FakeRateAnalyzer", &RegistrableAnalysis<FakeRateAnalyzer>::create);
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

