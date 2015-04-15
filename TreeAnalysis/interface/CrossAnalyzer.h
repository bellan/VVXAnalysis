#ifndef CrossAnalyzer_h
#define CrossAnalyzer_h

/** \class CrossAnalyzer
 *  Concrete class for Cross analysis
 *
 *  $Date: 2013/03/15 13:37:31 $
 *  $Revision: 1.4 $
 *  \author R. Bellan - UNITO <riccardo.bellan@cern.ch>
 */


#include "EventAnalyzer.h"
#include "RegistrableAnalysis.h"
#include "VVXAnalysis/Commons/interface/Constants.h"

class CrossAnalyzer: public EventAnalyzer, RegistrableAnalysis<CrossAnalyzer>{

public:
 CrossAnalyzer(const AnalysisConfiguration& configuration)
   : EventAnalyzer(*(new Selector<CrossAnalyzer>(*this)),
		   configuration){}
  
  virtual ~CrossAnalyzer(){}

  void ZZplots(int id = -1);

  virtual void analyze();

  virtual void begin();

  virtual void end( TFile &);

 private:

  friend class Selector<CrossAnalyzer>;

  std::vector<double> Xbins; 
  std::vector<double> Ybins; 

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

