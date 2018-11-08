#ifndef VVXAnalysis_TreeAnalysis_RegistrableAnalysis_h
#define VVXAnalysis_TreeAnalysis_RegistrableAnalysis_h

#include "AnalysisConfiguration.h"


class EventAnalyzer;


template <typename T>
class RegistrableAnalysis{
 public:
  
  RegistrableAnalysis(){}
  
  static EventAnalyzer* create(const AnalysisConfiguration& configuration) {  
    return new T(configuration);
  }
  
    
  
};


#endif
