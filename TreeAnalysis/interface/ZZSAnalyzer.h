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
#include "VVXAnalysis/Commons/interface/Constants.h"

class ZZSAnalyzer: public EventAnalyzer, RegistrableAnalysis<ZZSAnalyzer>{

public:
 ZZSAnalyzer(std::string filename, double lumi = 1., double externalXSection = -1., bool doBasicPlots = false)
   : EventAnalyzer(*(new Selector<ZZSAnalyzer>(*this)),
		   filename, lumi, externalXSection, doBasicPlots){}
  
  virtual ~ZZSAnalyzer(){}

  virtual void analyze();

 private:
  friend class Selector<ZZSAnalyzer>;
  template< class PAR >
    bool ZBosonDefinition(phys::Boson<PAR> cand) const{
    return fabs(cand.p4().M() - phys::ZMASS) < 20;
  }
  template< class PAR >
    bool WBosonDefinition(phys::Boson<PAR> cand) const{
    return fabs(cand.p4().M() - phys::WMASS) < 40;
  }
};
#endif

