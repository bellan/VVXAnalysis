#ifndef pippo
#define pippo

class EventAnalyzer;

template <typename T>
class RegistrableAnalysis{
 public:

  RegistrableAnalysis(){}
  
  static EventAnalyzer* create(std::string filename, double lumi, double externalXSection, bool doBasicPlots) {  
    return new T(filename, lumi, externalXSection, doBasicPlots);
  }
  
    
  
};


#endif
