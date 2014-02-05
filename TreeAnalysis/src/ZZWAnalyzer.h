#ifndef ZZWAnalyzer_h
#define ZZWAnalyzer_h

#include "EventAnalyzer.h"

class ZZWAnalyzer: public EventAnalyzer{

public:
  ZZWAnalyzer(std::string filename, double lumi = 1., double externalXSection = -1.)
    : EventAnalyzer(filename, lumi, externalXSection){}
  // virtual ~ZZWAnalyzer(){}
  virtual void analyze();
};
#endif

