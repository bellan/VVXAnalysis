#ifndef VVXAnalyzer_h
#define VVXAnalyzer_h

/** \class VVXAnalyzer
 *  Concrete class for VVX analysis
 *
 *  $Date: 2013/03/15 13:37:31 $
 *  $Revision: 1.4 $
 *  \author R. Bellan - UNITO <riccardo.bellan@cern.ch>
 */


#include "EventAnalyzer.h"
#include "AnalysisFactory.h"

class VVXAnalyzer: public EventAnalyzer{

public:
 VVXAnalyzer(std::string filename, double lumi = 1., double externalXSection = -1., bool doBasicPlots = false)
    : EventAnalyzer(filename, lumi, externalXSection, doBasicPlots){}
  virtual ~VVXAnalyzer(){}
  virtual void analyze();


  static EventAnalyzer* create(std::string filename, double lumi, double externalXSection, bool doBasicPlots) {
    return new VVXAnalyzer(filename, lumi, externalXSection, doBasicPlots);
  }

  virtual void Register(std::string analyisName) {
    AnalysisFactory::get()->Register("VVXAnalyzer", &create);
  }

};
#endif

