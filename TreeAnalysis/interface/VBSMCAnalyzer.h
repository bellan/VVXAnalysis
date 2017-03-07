#ifndef VBSMCAnalyzer_h
#define VBSMCAnalyzer_h

//Analyzer for MC puropose. Generate plots for MC CrossSection and Acceptance

#include "EventAnalyzer.h"
#include "RegistrableAnalysis.h"
#include "VVXAnalysis/Commons/interface/Constants.h"
#include "VVXAnalysis/Commons/interface/LeptonScaleFactors.h"
#include <TTree.h>

class VBSMCAnalyzer: public EventAnalyzer, RegistrableAnalysis<VBSMCAnalyzer>{

public:
 VBSMCAnalyzer(const AnalysisConfiguration& configuration)
   : EventAnalyzer(*(new Selector<VBSMCAnalyzer>(*this)),
		   configuration){}
    
  virtual ~VBSMCAnalyzer(){}

  void ZZplots(std::string decay);
  virtual void analyze();
  virtual void begin();
  virtual void end( TFile &);

  float mjj_gen;
  float deta_gen;

  float m4L_gen;
  int   njets;
  
 private:

  friend class Selector<VBSMCAnalyzer>;

  std::vector<double> Xbins;

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
