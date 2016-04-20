#ifndef ZZRecoAnalyzer_h
#define ZZRecoAnalyzer_h

/** \class ZZRecoAnalyzer
 *  Concrete class for ZZReco analysis
 *
 *  $Date: 2013/03/15 13:37:31 $
 *  $Revision: 1.4 $
 *  \author R. Bellan - UNITO <riccardo.bellan@cern.ch>
 */


#include "EventAnalyzer.h"
#include "RegistrableAnalysis.h"
#include "VVXAnalysis/Commons/interface/Constants.h"
#include "VVXAnalysis/Commons/interface/LeptonScaleFactors.h"

class ZZRecoAnalyzer: public EventAnalyzer, RegistrableAnalysis<ZZRecoAnalyzer>{
  
 public:
 ZZRecoAnalyzer(const AnalysisConfiguration& configuration)
   : EventAnalyzer(*(new Selector<ZZRecoAnalyzer>(*this)),
		   configuration),lepSF("../../ZZAnalysis/AnalysisStep/test/Macros/scale_factors_muons2012.root",
					"../../ZZAnalysis/AnalysisStep/test/Macros/scale_factors_ele2012.root",
					"../Commons/data/fakeRates_mu.root",
					"../Commons/data/fakeRates_el.root"){}
  
  virtual ~ZZRecoAnalyzer(){}
  
  void ZZplots(int id = -1, int e = 0);
  
  virtual void analyze();
  
  virtual void begin();
  
  virtual void end( TFile &);

  int e;

  Long64_t nentries;

  float m4L;
  float m4L_gen;
  float drzz;
  float drzz_gen;
  int ngenjets;
  float mjj_gen;
  float deta_gen;
  int ngencentraljets;
  float mjj_gen_cj;
  float deta_gen_cj;
  float ptjet1_gen;
  float ptjet2_gen; 
  float etajet1_gen;
  float etajet2_gen;
  float ptZZ;
  float min_dR;
  Int_t nUpJESDatajets;
  Int_t  nDownJESDatajets;
  Int_t nUpJESDatacentraljets;
  Int_t nDownJESDatacentraljets;
  Int_t nCentralJERjets;
  Int_t nUpJERjets;
  Int_t nDownJERjets; 
  Int_t nUpJESjets; 
  Int_t nDownJESjets;
  Int_t nCentralJERcentraljets;
  Int_t nUpJERcentraljets;
  Int_t nDownJERcentraljets; 
  Int_t nUpJEScentraljets; 
  Int_t nDownJEScentraljets;
  Int_t jet_count;
  Int_t inFiducialRegion;
  Int_t inGenAndReco;
  float w_kf;
  float dphizz;
  float dphizz_gen;
  float ptzz;
  float ptzz_gen;
  int nEvent;

  double JER_PtSmear(double pt, double width);
 // Jets obtained by gaussian JER smearing
  std::vector<phys::Jet> *CentralJER_jets;
  std::vector<phys::Jet> *UpJER_jets;
  std::vector<phys::Jet> *DownJER_jets;
  std::vector<phys::Jet> *CentralJER_centraljets;
  std::vector<phys::Jet> *UpJER_centraljets;
  std::vector<phys::Jet> *DownJER_centraljets;

  // Jets obtained correcting up and down for the JES uncertainty
  std::vector<phys::Jet> *UpJES_jets;
  std::vector<phys::Jet> *DownJES_jets;
  std::vector<phys::Jet> *UpJES_centraljets;
  std::vector<phys::Jet> *DownJES_centraljets;

  // Jets obtained correcting up and down for the JES uncertainty (data distributions = no JER correction applied)
  std::vector<phys::Jet> *UpJESData_jets;
  std::vector<phys::Jet> *DownJESData_jets;
  std::vector<phys::Jet> *UpJESData_centraljets;
  std::vector<phys::Jet> *DownJESData_centraljets;

  std::vector<double> *CentralJER_jetPt;
  std::vector<double> *UpJER_jetPt;
  std::vector<double> *DownJER_jetPt;
  std::vector<double> *UpJES_jetPt;
  std::vector<double> *DownJES_jetPt;
  std::vector<double> *UpJESData_jetPt;
  std::vector<double> *DownJESData_jetPt;


  TFile * UnfOverMC;
  TFile * UnfOverMC_Pow;
  TH1 * h_UnfOverMC_Mass; 
  TH1 * h_UnfOverMC_Jets; 
  TH1 * h_UnfOverMC_Mjj; 
  TH1 * h_UnfOverMC_Deta; 
  

 private:
  
  friend class Selector<ZZRecoAnalyzer>;
  
  std::vector<double> Xbins; 
  std::vector<double> Ybins; 
  std::vector<double> Xbins_deta;
  std::vector<double> Xbins_mjj;
  std::vector<double> Xbins_ptjet1;
  std::vector<double> Xbins_ptjet2;
  std::vector<double> Xbins_etajet1;
  std::vector<double> Xbins_etajet2;
  std::vector<double> Xbins_ptZZ;
  std::vector<double> Xbins_dphi; 
  std::vector<double> Xbins_drzz; 
  std::vector<double> Xbins_ptzz;
  std::vector<double> Xbins_dphizz;
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

