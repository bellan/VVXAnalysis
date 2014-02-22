#include "VVXAnalysis/TreeAnalysis/interface/AnalysisFactory.h"
#include "VVXAnalysis/TreeAnalysis/interface/VVXAnalyzer.h"

AnalysisFactory::AnalysisFactory(){
  Register("VVXAnalyzer", &VVXAnalyzer::create);
}

void AnalysisFactory::Register(const std::string &analysisName, EventAnalyzer::CreateAnFn pfnCreate)
{
    m_FactoryMap[analysisName] = pfnCreate;
}

EventAnalyzer *AnalysisFactory::createAnalysis(const std::string &analysisName, std::string filename, double lumi, double externalXSection, bool doBasicPlots)
{
  FactoryMap::iterator it = m_FactoryMap.find(analysisName);
  if( it != m_FactoryMap.end() )
    return it->second(filename,lumi, externalXSection, doBasicPlots);
 
  return NULL;
}
