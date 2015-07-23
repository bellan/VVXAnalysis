#ifndef ZZMCAnalyzer_h
#define ZZMCAnalyzer_h

//Analyzer for MC puropose. Generate plots for MC CrossSection and Acceptance

#include "EventAnalyzer.h"
#include "RegistrableAnalysis.h"
#include "VVXAnalysis/Commons/interface/Constants.h"
#include "VVXAnalysis/Commons/interface/LeptonScaleFactors.h"
#include <TTree.h>

class ZZMCAnalyzer: public EventAnalyzer, RegistrableAnalysis<ZZMCAnalyzer>{

public:
 ZZMCAnalyzer(const AnalysisConfiguration& configuration)
   : EventAnalyzer(*(new Selector<ZZMCAnalyzer>(*this)),
		   configuration),lepSF("../../ZZAnalysis/AnalysisStep/test/Macros/scale_factors_muons2012.root",
					"../../ZZAnalysis/AnalysisStep/test/Macros/scale_factors_ele2012.root",
					"../Commons/data/fakeRates_mu.root",
					"../Commons/data/fakeRates_el.root"){}
  
  
  virtual ~ZZMCAnalyzer(){}

  void ZZplots(std::string decay);

  virtual void analyze();

  virtual void begin();

  virtual void end( TFile &);

  Long64_t nentries;
  float m4L_gen;
  int njets;
  int PreCounter;

  TFile * UnfOverMC;
  TH1 * h_UnfOverMC_Mass; 
  TH1 * h_UnfOverMC_Jets;

 private:

  friend class Selector<ZZMCAnalyzer>;

  std::vector<double> Xbins; 
  std::vector<double> Ybins; 

  LeptonScaleFactors lepSF;

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
