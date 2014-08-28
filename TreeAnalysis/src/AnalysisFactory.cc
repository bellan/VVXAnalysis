#include "VVXAnalysis/TreeAnalysis/interface/AnalysisFactory.h"
#include "VVXAnalysis/TreeAnalysis/interface/VVXAnalyzer.h"
#include "VVXAnalysis/TreeAnalysis/interface/ZZWAnalyzer.h"
#include "VVXAnalysis/TreeAnalysis/interface/ZZSAnalyzer.h"

AnalysisFactory::AnalysisFactory(){
  Register("VVXAnalyzer", &RegistrableAnalysis<VVXAnalyzer>::create);
  Register("ZZWAnalyzer", &RegistrableAnalysis<ZZWAnalyzer>::create);
  Register("ZZSAnalyzer", &RegistrableAnalysis<ZZSAnalyzer>::create);
}

void AnalysisFactory::Register(const std::string &analysisName, CreateAnFn pfnCreate)
{
    m_FactoryMap[analysisName] = pfnCreate;
}

//EventAnalyzer *AnalysisFactory::createAnalysis(const std::string &analysisName, const std::string& region, const std::string& filename, const double& lumi, const double& externalXSection, bool doBasicPlots)
EventAnalyzer *AnalysisFactory::createAnalysis(const AnalysisConfiguration &analysisConfiguration)
{
  FactoryMap::iterator it = m_FactoryMap.find(analysisConfiguration.getParameter<std::string>("analysisName"));
  if( it != m_FactoryMap.end() )
    return it->second(analysisConfiguration);
 
  return NULL;
}
