#ifndef VVXAnalysis_TreeAnalysis_RegistrableAnalysis_h
#define VVXAnalysis_TreeAnalysis_RegistrableAnalysis_h

class EventAnalyzer;

template <typename T>
class RegistrableAnalysis{
 public:

  RegistrableAnalysis(){}
  
  static EventAnalyzer* create(const std::string& region, const std::string& filename, const double& lumi, const double& externalXSection, bool doBasicPlots) {  
    return new T(region, filename, lumi, externalXSection, doBasicPlots);
  }
  
    
  
};


#endif
