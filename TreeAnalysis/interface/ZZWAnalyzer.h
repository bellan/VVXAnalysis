#ifndef ZZWAnalyzer_h
#define ZZWAnalyzer_h

/** \class ZZWAnalyzer
 *  Concrete class for ZZW analysis
 *
 *  $Date: 2013/03/15 13:37:31 $
 *  $Revision: 1.4 $
 *  \author R. Bellan - UNITO <riccardo.bellan@cern.ch>
 */


#include "EventAnalyzer.h"

class ZZWAnalyzer: public EventAnalyzer{

public:
 ZZWAnalyzer(std::string filename, double lumi = 1., double externalXSection = -1., bool doBasicPlots = false)
    : EventAnalyzer(filename, lumi, externalXSection, doBasicPlots){}

  virtual ~ZZWAnalyzer(){}

  static EventAnalyzer* create(std::string filename, double lumi, double externalXSection, bool doBasicPlots) {  
    return new ZZWAnalyzer(filename, lumi, externalXSection, doBasicPlots);
  }

  virtual void analyze();

  virtual Int_t cut();
};
#endif

