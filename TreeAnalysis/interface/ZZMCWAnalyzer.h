#ifndef ZZMCWAnalyzer_h
#define ZZMCWAnalyzer_h

//Analyzer for MC puropose. Generate plots for MC CrossSection and Acceptance

#include "EventAnalyzer.h"
#include "RegistrableAnalysis.h"
#include "VVXAnalysis/Commons/interface/Constants.h"
#include "VVXAnalysis/Commons/interface/LeptonScaleFactors.h"
#include <TTree.h>

class ZZMCWAnalyzer: public EventAnalyzer, RegistrableAnalysis<ZZMCWAnalyzer>{

public:
 ZZMCWAnalyzer(const AnalysisConfiguration& configuration)
   : EventAnalyzer(*(new Selector<ZZMCWAnalyzer>(*this)),
		   configuration),lepSF("../../ZZAnalysis/AnalysisStep/test/Macros/scale_factors_muons2012.root",
  			   "../../ZZAnalysis/AnalysisStep/test/Macros/scale_factors_ele2012.root",
  			   "../data/fakeRates_mu.root",
  			   "../data/fakeRates_el.root"){}
  
  
  virtual ~ZZMCWAnalyzer(){}

  void ZZplots(std::string decay);

  virtual void analyze();

  virtual void begin();

  virtual void end( TFile &);

  Long64_t nentries;
  float m4L_gen;
  float mjj_gen;
  float deta_gen;
  float mjj_gen_cj;
  float deta_gen_cj; 
  float ptjet1_gen;
  float ptjet2_gen;
  float etajet1_gen;
  float etajet2_gen;
  int njets;
  int ncentraljets;
  int PreCounter;
 

  int bin_Jets;
  int bin_Mass;
  int bin_Mjj;
  int bin_Deta;
  int bin_CentralJets;
  int bin_CentralMjj;
  int bin_CentralDeta;  
  int bin_PtJet1;
  int bin_PtJet2;
  int bin_EtaJet1; 
  int bin_EtaJet2;
  float w_Jets;
  float w_Mass;
  float w_Mjj;
  float w_Deta;
  float w_CentralJets;
  float w_CentralMjj;
  float w_CentralDeta; 
  float w_PtJet1;
  float w_PtJet2;
  float w_EtaJet1;
  float w_EtaJet2;

  int bin_Jets_fr;
  int bin_Mass_fr;
  int bin_Mjj_fr;
  int bin_Deta_fr;
  int bin_CentralJets_fr;
  int bin_CentralMjj_fr;
  int bin_CentralDeta_fr;  
  int bin_PtJet1_fr;
  int bin_PtJet2_fr;
  int bin_EtaJet1_fr; 
  int bin_EtaJet2_fr;
  float w_Jets_fr;
  float w_Mass_fr;
  float w_Mjj_fr;
  float w_Deta_fr;
  float w_CentralJets_fr;
  float w_CentralMjj_fr;
  float w_CentralDeta_fr; 
  float w_PtJet1_fr;
  float w_PtJet2_fr;
  float w_EtaJet1_fr;
  float w_EtaJet2_fr;

  TFile * UnfOverMC;
  TFile * UnfOverMC_Pow;
  TFile * UnfOverMC_fr;
  TFile * UnfOverMC_fr_Pow;

  TH1 * h_UnfOverMC_Mass_4e; 
  TH1 * h_UnfOverMC_Jets_4e; 
  TH1 * h_UnfOverMC_Mjj_4e; 
  TH1 * h_UnfOverMC_Deta_4e; 
  TH1 * h_UnfOverMC_CentralJets_4e; 
  TH1 * h_UnfOverMC_CentralMjj_4e; 
  TH1 * h_UnfOverMC_CentralDeta_4e; 
  TH1 * h_UnfOverMC_PtJet1_4e; 
  TH1 * h_UnfOverMC_PtJet2_4e;
  TH1 * h_UnfOverMC_EtaJet1_4e; 
  TH1 * h_UnfOverMC_EtaJet2_4e;
  TH1 * h_UnfOverMC_Mass_2e2m; 
  TH1 * h_UnfOverMC_Jets_2e2m; 
  TH1 * h_UnfOverMC_Mjj_2e2m; 
  TH1 * h_UnfOverMC_Deta_2e2m; 
  TH1 * h_UnfOverMC_PtJet1_2e2m; 
  TH1 * h_UnfOverMC_PtJet2_2e2m;
  TH1 * h_UnfOverMC_EtaJet1_2e2m; 
  TH1 * h_UnfOverMC_EtaJet2_2e2m; 
  TH1 * h_UnfOverMC_CentralJets_2e2m; 
  TH1 * h_UnfOverMC_CentralMjj_2e2m; 
  TH1 * h_UnfOverMC_CentralDeta_2e2m; 
  TH1 * h_UnfOverMC_Mass_4m; 
  TH1 * h_UnfOverMC_Jets_4m; 
  TH1 * h_UnfOverMC_Mjj_4m; 
  TH1 * h_UnfOverMC_Deta_4m; 
  TH1 * h_UnfOverMC_PtJet1_4m; 
  TH1 * h_UnfOverMC_PtJet2_4m;
  TH1 * h_UnfOverMC_EtaJet1_4m; 
  TH1 * h_UnfOverMC_EtaJet2_4m; 
  TH1 * h_UnfOverMC_CentralJets_4m; 
  TH1 * h_UnfOverMC_CentralMjj_4m; 
  TH1 * h_UnfOverMC_CentralDeta_4m;
  TH1 * h_UnfOverMC_Mass; 
  TH1 * h_UnfOverMC_Jets; 
  TH1 * h_UnfOverMC_Mjj; 
  TH1 * h_UnfOverMC_Deta; 
  TH1 * h_UnfOverMC_CentralJets; 
  TH1 * h_UnfOverMC_CentralMjj; 
  TH1 * h_UnfOverMC_CentralDeta; 

  TH1 * h_UnfOverMC_fr_Mass_4e; 
  TH1 * h_UnfOverMC_fr_Jets_4e; 
  TH1 * h_UnfOverMC_fr_Mjj_4e; 
  TH1 * h_UnfOverMC_fr_Deta_4e; 
  TH1 * h_UnfOverMC_fr_CentralJets_4e; 
  TH1 * h_UnfOverMC_fr_CentralMjj_4e; 
  TH1 * h_UnfOverMC_fr_CentralDeta_4e; 
  TH1 * h_UnfOverMC_fr_PtJet1_4e; 
  TH1 * h_UnfOverMC_fr_PtJet2_4e;
  TH1 * h_UnfOverMC_fr_EtaJet1_4e; 
  TH1 * h_UnfOverMC_fr_EtaJet2_4e;
  TH1 * h_UnfOverMC_fr_Mass_2e2m; 
  TH1 * h_UnfOverMC_fr_Jets_2e2m; 
  TH1 * h_UnfOverMC_fr_Mjj_2e2m; 
  TH1 * h_UnfOverMC_fr_Deta_2e2m; 
  TH1 * h_UnfOverMC_fr_PtJet1_2e2m; 
  TH1 * h_UnfOverMC_fr_PtJet2_2e2m;
  TH1 * h_UnfOverMC_fr_EtaJet1_2e2m; 
  TH1 * h_UnfOverMC_fr_EtaJet2_2e2m; 
  TH1 * h_UnfOverMC_fr_CentralJets_2e2m; 
  TH1 * h_UnfOverMC_fr_CentralMjj_2e2m; 
  TH1 * h_UnfOverMC_fr_CentralDeta_2e2m; 
  TH1 * h_UnfOverMC_fr_Mass_4m; 
  TH1 * h_UnfOverMC_fr_Jets_4m; 
  TH1 * h_UnfOverMC_fr_Mjj_4m; 
  TH1 * h_UnfOverMC_fr_Deta_4m; 
  TH1 * h_UnfOverMC_fr_PtJet1_4m; 
  TH1 * h_UnfOverMC_fr_PtJet2_4m;
  TH1 * h_UnfOverMC_fr_EtaJet1_4m; 
  TH1 * h_UnfOverMC_fr_EtaJet2_4m; 
  TH1 * h_UnfOverMC_fr_CentralJets_4m; 
  TH1 * h_UnfOverMC_fr_CentralMjj_4m; 
  TH1 * h_UnfOverMC_fr_CentralDeta_4m;
  TH1 * h_UnfOverMC_fr_Mass; 
  TH1 * h_UnfOverMC_fr_Jets; 
  TH1 * h_UnfOverMC_fr_Mjj; 
  TH1 * h_UnfOverMC_fr_Deta; 
  TH1 * h_UnfOverMC_fr_CentralJets; 
  TH1 * h_UnfOverMC_fr_CentralMjj; 
  TH1 * h_UnfOverMC_fr_CentralDeta; 
 private:

  friend class Selector<ZZMCWAnalyzer>;

  std::vector<double> Xbins;
  std::vector<double> Xbins_mjj;  
  std::vector<double> Xbins_deta;
  std::vector<double> Ybins; 
  std::vector<double> Xbins_ptjet1;
  std::vector<double> Xbins_ptjet2;
  std::vector<double> Xbins_etajet1;
  std::vector<double> Xbins_etajet2;

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
