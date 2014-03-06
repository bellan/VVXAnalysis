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
#include "RegistrableAnalysis.h"

class ZZSAnalyzer: public EventAnalyzer, RegistrableAnalysis<ZZSAnalyzer>{

public:
 ZZSAnalyzer(std::string filename, double lumi = 1., double externalXSection = -1., bool doBasicPlots = false)
   : EventAnalyzer(*(new Selector<ZZSAnalyzer>(*this)),
		   filename, lumi, externalXSection, doBasicPlots){}
  
  virtual ~ZZSAnalyzer(){}

  virtual void analyze();

 private:
  template< class PAR >
    bool bosonDefinition(phys::Boson<PAR> vb) const { 
    return true;
  }
  friend class Selector<ZZSAnalyzer>;

};
#endif

