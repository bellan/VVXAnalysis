#ifndef ZZRecoAnalyzer_h
#define ZZRecoAnalyzer_h

/** \class ZZRecoAnalyzer
 *  Concrete class for ZZReco analysis
 *
 *  $Date: 2016/06/09 $
 *  $Revision: 1.4 $
 *  \author G. L. Pinna Angioni - UNITO <Gian.Luca.Pinna.Angioni@cern.ch>
 */


#include "EventAnalyzer.h"
#include "RegistrableAnalysis.h"
#include "VVXAnalysis/Commons/interface/Constants.h"
#include "VVXAnalysis/Commons/interface/LeptonScaleFactors.h"

class ZZRecoAnalyzer: public EventAnalyzer, RegistrableAnalysis<ZZRecoAnalyzer>{
  
 public:
 ZZRecoAnalyzer(const AnalysisConfiguration& configuration)
   : EventAnalyzer(*(new Selector<ZZRecoAnalyzer>(*this)),
		   configuration){}
  
  virtual ~ZZRecoAnalyzer(){}
  
  
  void FillHistosBase(std::string decay, float Wh, std::string type );
  
  void FillHistosJets(std::string decay, float Wh, std::vector<phys::Jet> *jetsVec ,std::string type);
  
  void FillMatrixHistosBase(std::string decay, float Wh,std::string type );  
  
  void FillMatrixHistosJets(std::string decay, float Wh,std::vector<phys::Jet> *jetsVec,std::vector<phys::Particle> *jetsGenVec,std::string type); 
  virtual void analyze();
  
  virtual void begin();
  
  virtual void end( TFile &);

  int e;

  Long64_t nentries;

  float m4L;
  float m4L_gen;
  float drzz;

  int njets;
  float ptJet1;
  float ptJet2;
  float ptJet3;
  float etaJet1;
  float etaJet2;
  float deta;
  float mjj;

  float drzz_gen;
  int njets_gen;
  float mjj_gen;
  float deta_gen;
  int ngencentraljets;
  float mjj_gen_cj;
  float deta_gen_cj;
  float ptJet1_gen;
  float ptJet2_gen; 
  float etaJet1_gen;
  float etaJet2_gen;
  float ptZZ;
  float min_dR;
  Int_t nUpJERjets;
  Int_t nDownJERjets; 
  Int_t nUpJERcentraljets;
  Int_t nDownJERcentraljets; 

  float dphi;
  float dphi_gen;

  Int_t inFiducialRegion;
  float w_kf;
  float dphizz;
  float dphizz_gen;
  float ptzz;
  float ptzz_gen;
  int nEvent;

  TStopwatch *st;


 // Jets obtained by gaussian JER smearing
  std::vector<phys::Jet> *UpJER_jets;
  std::vector<phys::Jet> *DownJER_jets;
  std::vector<phys::Jet> *UpJER_centraljets;
  std::vector<phys::Jet> *DownJER_centraljets;

  // Jets obtained correcting up and down for the JES uncertainty
  std::vector<phys::Jet> *UpJES_jets;
  std::vector<phys::Jet> *DownJES_jets;
  std::vector<phys::Jet> *UpJES_centraljets;
  std::vector<phys::Jet> *DownJES_centraljets;

  //  Jets obtained correcting up and down for the JES uncertainty (data distributions = no JER correction applied)
  std::vector<phys::Jet> *UpJESData_jets;
  std::vector<phys::Jet> *DownJESData_jets;
  std::vector<phys::Jet> *UpJESData_centraljets;
  std::vector<phys::Jet> *DownJESData_centraljets;


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
  std::vector<double> Xbins_ptJet1;
  std::vector<double> Xbins_ptJet2;
  std::vector<double> Xbins_ptJet3;
  std::vector<double> Xbins_etaJet1;
  std::vector<double> Xbins_etaJet2;
  std::vector<double> Xbins_ptZZ;
  std::vector<double> Xbins_dphi; 
  std::vector<double> Xbins_drzz; 
  std::vector<double> Xbins_ptzz;
  std::vector<double> Xbins_dphizz;
  std::vector<double> Xbins_nJets;
 
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

