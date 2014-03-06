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
#include "RegistrableAnalysis.h"

class VVXAnalyzer: public EventAnalyzer, RegistrableAnalysis<VVXAnalyzer>{

public:

 VVXAnalyzer(std::string filename, double lumi = 1., double externalXSection = -1., bool doBasicPlots = false)
   : EventAnalyzer(*(new Selector<VVXAnalyzer>(*this)),
		   filename, lumi, externalXSection, doBasicPlots){}

  virtual ~VVXAnalyzer(){}

  virtual void analyze();

 private:
  friend class Selector<VVXAnalyzer>;
  template< class PAR >
    bool ZBosonDefinition(phys::Boson<PAR> cand) const{
    return fabs(cand.p4().M() - ZMASS) < 20;
  }
  template< class PAR >
    bool WBosonDefinition(phys::Boson<PAR> cand) const{
    return fabs(cand.p4().M() - WMASS) < 40;
  }
};
#endif

