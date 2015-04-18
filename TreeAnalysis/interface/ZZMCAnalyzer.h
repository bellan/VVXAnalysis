#ifndef ZZMCAnalyzer_h
#define ZZMCAnalyzer_h

/** \class ZZMCAnalyzer
 *  Concrete class for ZZMC analysis
 *
 *  $Date: 2013/03/15 13:37:31 $
 *  $Revision: 1.4 $
 *  \author R. Bellan - UNITO <riccardo.bellan@cern.ch>
 */


#include "EventAnalyzer.h"
#include "RegistrableAnalysis.h"
#include "VVXAnalysis/Commons/interface/Constants.h"
#include <TTree.h>

class ZZMCAnalyzer: public EventAnalyzer, RegistrableAnalysis<ZZMCAnalyzer>{

public:
 ZZMCAnalyzer(const AnalysisConfiguration& configuration)
   : EventAnalyzer(*(new Selector<ZZMCAnalyzer>(*this)),
		   configuration){}
  
  virtual ~ZZMCAnalyzer(){}

  void ZZplots(int id = -1);

  virtual void analyze();

  virtual void begin();

  virtual void end( TFile &);

  Long64_t nentries = 0;
  float m4L_gen = 0;
  int e =0;
 
 private:

  friend class Selector<ZZMCAnalyzer>;

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
