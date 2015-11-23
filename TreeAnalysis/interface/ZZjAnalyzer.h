#ifndef ZZjAnalyzer_h
#define ZZjAnalyzer_h

/** \class ZZjAnalyzer
 *  Concrete class for ZZj analysis
 *
 *  $Date: 2013/03/15 13:37:31 $
 *  $Revision: 1.4 $
 *  \author R. Bellan - UNITO <riccardo.bellan@cern.ch>
 */


#include "EventAnalyzer.h"
#include "RegistrableAnalysis.h"
#include "VVXAnalysis/Commons/interface/Constants.h"
#include "VVXAnalysis/Commons/interface/LeptonScaleFactors.h"
class ZZjAnalyzer: public EventAnalyzer, RegistrableAnalysis<ZZjAnalyzer>{

public:
 ZZjAnalyzer(const AnalysisConfiguration& configuration)
   : EventAnalyzer(*(new Selector<ZZjAnalyzer>(*this)),
		   configuration),lepFR("../../ZZAnalysis/AnalysisStep/test/Macros/scale_factors_muons2012.root",
					"../../ZZAnalysis/AnalysisStep/test/Macros/scale_factors_ele2012.root",
					"../Commons/data/fakeRates.root",
					"../Commons/data/fakeRates.root"){}
  
  virtual ~ZZjAnalyzer(){}

  void ZZplots(int id = -1);

  virtual void analyze();

  virtual void end( TFile &);
 

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

 LeptonScaleFactors lepFR;

  std::vector<std::string> events2e2mu; 
  std::vector<std::string> events4e; 
  std::vector<std::string> events4mu;


};
#endif

