#ifndef VVXAnalysis_TreeAnalysis_AnalysisFactory_h
#define VVXAnalysis_TreeAnalysis_AnalysisFactory_h

/** \class AnalysisFactory
 *  Creates concrete class for VVX analysis
 *
 *  $Date: 2013/03/15 13:37:31 $
 *  $Revision: 1.4 $
 *  \author R. Bellan - UNITO <riccardo.bellan@cern.ch>
 */

#include "AnalysisConfiguration.h"
#include "EventAnalyzer.h"
#include <map>
#include <string>

class AnalysisFactory{
private:
    AnalysisFactory();
    AnalysisFactory(const AnalysisFactory &) { }
    AnalysisFactory &operator=(const AnalysisFactory &) { return *this; }

    typedef EventAnalyzer* (*CreateAnFn)(const AnalysisConfiguration& analysisConfiguration);
    typedef std::map<std::string, CreateAnFn> FactoryMap;
    FactoryMap m_FactoryMap;
public:
    ~AnalysisFactory() { m_FactoryMap.clear(); }

    static AnalysisFactory *get()
    {
        static AnalysisFactory instance;
        return &instance;
    }

    void Register(const std::string &analysisName, CreateAnFn pfnCreate);
    //EventAnalyzer *createAnalysis(const std::string &analysisName, const std::string& region, const std::string& filename, const double& lumi = 1., const double& externalXSection = -1., bool doBasicPlots = true);
    EventAnalyzer *createAnalysis(const AnalysisConfiguration &analysisConfiguration);

};
#endif
