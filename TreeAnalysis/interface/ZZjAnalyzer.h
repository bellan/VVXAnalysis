#ifndef ZZjAnalyzer_h
#define ZZjAnalyzer_h

/** \class ZZjAnalyzer
 *  Concrete class for ZZj analysis
 *
 *  $Date: 2013/03/15 13:37:31 $
 *  $Revision: 1.4 $
 *  \author R. Bellan - UNITO <riccardo.bellan@cern.ch>
 */

#include "TRandom.h"
#include "EventAnalyzer.h"
#include "RegistrableAnalysis.h"
#include "VVXAnalysis/Commons/interface/Constants.h"
#include "VVXAnalysis/Commons/interface/LeptonScaleFactors.h"
class ZZjAnalyzer: public EventAnalyzer, RegistrableAnalysis<ZZjAnalyzer>{

public:
 ZZjAnalyzer(const AnalysisConfiguration& configuration)
   : EventAnalyzer(*(new Selector<ZZjAnalyzer>(*this)),
		   configuration){}
  
  virtual ~ZZjAnalyzer(){}

  void ZZplots(int id = -1);

  virtual void analyze();

  virtual void end( TFile &);

  // Jets obtained by gaussian JER smearing

  std::vector<phys::Jet> *  noJerJets;
  std::vector<double>    *  noJerPt;

 private:
  friend class Selector<ZZjAnalyzer>;
  template< class PAR >
    bool ZBosonDefinition(phys::Boson<PAR> cand) const{
    return fabs(cand.p4().M() - phys::ZMASS) < 20;
  }
  template< class PAR >
    bool WBosonDefinition(phys::Boson<PAR> cand) const{
    return fabs(cand.p4().M() - phys::WMASS) < 40;
  }

  std::vector<std::string> events2e2mu; 
  std::vector<std::string> events4e; 
  std::vector<std::string> events4mu;
  std::vector<std::string> events4l;

};
#endif

