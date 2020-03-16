#ifndef ZZjGenAnalyzer_h
#define ZZjGenAnalyzer_h

/** \class ZZjGenAnalyzer
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
class ZZjGenAnalyzer: public EventAnalyzer, RegistrableAnalysis<ZZjGenAnalyzer>{

public:
 ZZjGenAnalyzer(const AnalysisConfiguration& configuration)
   : EventAnalyzer(*(new Selector<ZZjGenAnalyzer>(*this)),
		   configuration){}
  
  virtual ~ZZjGenAnalyzer(){}

  void ZZplots(std::string decay);
  virtual void analyze();
  virtual void end( TFile &);

  float m4L_gen;
  float mZ1_gen;
  float mZ2_gen;
  float mjj_gen;
  float deta_gen;
  float mjj_gen_cj;
  float deta_gen_cj;
  float ptjet1_gen;
  float ptjet2_gen;
  float etajet1_gen;
  float etajet2_gen; 
  float drzz_gen;
  int njets;
 

 private:
  friend class Selector<ZZjGenAnalyzer>;
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
  std::vector<std::string> eventsFull;

};
#endif

