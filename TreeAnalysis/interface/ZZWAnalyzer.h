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
  ZZWAnalyzer(std::string filename, double lumi = 1., double externalXSection = -1.)
    : EventAnalyzer(filename, lumi, externalXSection, true){}
  virtual ~ZZWAnalyzer(){}
  virtual void analyze();
};
#endif

