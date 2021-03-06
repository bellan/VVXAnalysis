#ifndef ZZMCAnalyzer_h
#define ZZMCAnalyzer_h

//Analyzer for MC puropose. Generate plots for MC CrossSection and Acceptance

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

  void FillHistosBase(std::string decay, float Wh, std::string type );

  virtual void analyze();

  virtual void begin();

  virtual void end( TFile &);

  Long64_t nentries;
  float m4L_gen;
  float pt4l_gen; 
  float mjj_gen;
  float deta_gen;
  float mjj_gen_cj;
  float deta_gen_cj;
  float ptjet1_gen;
  float ptjet2_gen;
  float etajet1_gen;
  float etajet2_gen; 
  float drzz_gen;
  float w_kf;
  float dphizz_gen;
  float ptzz_gen;
  Float_t scaleFacErrSq;
  Float_t scaleFacMuErrSq;
  Float_t scaleFacEleErrSq;


  int   njets;
  int   ncentraljets;
  int   PreCounter;
  int   inFiducialRegion;
  int   nEvent;

  std::string region;
  std::string sample;

  bool isTightFr;  
  std::vector<Float_t>  ScalVarVal;
  std::vector<std::string>  ScalVar;

 private:

  friend class Selector<ZZMCAnalyzer>;

  std::vector<double> Xbins_single;
  std::vector<double> Xbins;
  std::vector<double> Xbins_nJets;
  std::vector<double> Xbins_mjj;  
  std::vector<double> Xbins_deta;
  std::vector<double> Ybins; 
  std::vector<double> Xbins_ptjet1;
  std::vector<double> Xbins_ptjet2;
  std::vector<double> Xbins_etajet1;
  std::vector<double> Xbins_etajet2; 
  std::vector<double> Xbins_drzz;
  std::vector<double> Xbins_ptzz;
  std::vector<double> Xbins_dphizz;


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
