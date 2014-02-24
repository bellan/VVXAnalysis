#ifndef ZZSAnalyzer_h
#define ZZSAnalyzer_h

/** \class ZZSAnalyzer
 *  Concrete class for ZZS analysis
 *
 *  $Date: 2013/03/15 13:37:31 $
 *  $Revision: 1.4 $
 *  \author R. Bellan - UNITO <riccardo.bellan@cern.ch>
 */


#include "EventAnalyzer.h"

class ZZSAnalyzer: public EventAnalyzer{

public:
 ZZSAnalyzer(std::string filename, double lumi = 1., double externalXSection = -1., bool doBasicPlots = false)
    : EventAnalyzer(filename, lumi, externalXSection, doBasicPlots){}

  virtual ~ZZSAnalyzer(){}

  static EventAnalyzer* create(std::string filename, double lumi, double externalXSection, bool doBasicPlots) {  
    return new ZZSAnalyzer(filename, lumi, externalXSection, doBasicPlots);
  }

  virtual void analyze();
};
#endif

