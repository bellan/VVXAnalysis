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

class ZZRecoAnalyzer: public EventAnalyzer, RegistrableAnalysis<ZZRecoAnalyzer>{
  
 public:
 ZZRecoAnalyzer(const AnalysisConfiguration& configuration)
   : EventAnalyzer(*(new Selector<ZZRecoAnalyzer>(*this)),
		   configuration){}
  
  virtual ~ZZRecoAnalyzer(){}
  
  void ZZplots(int id = -1, int e = 0);
  
  virtual void analyze();
  
  virtual void begin();
  
  virtual void end( TFile &);

  int e;

  Long64_t nentries;

  float m4L;
  float m4L_gen;
  int ngenjets;
  float mjj_gen;
  float deta_gen;
  Int_t nCentralJERjets;
  Int_t nUpJERjets;
  Int_t nDownJERjets; 
  Int_t nUpJESjets; 
  Int_t nDownJESjets;

  double JER_PtSmear(double pt, double width);
  // Jets obtained by gaussian JER smearing
  std::vector<phys::Jet> *CentralJER_jets;
  std::vector<phys::Jet> *UpJER_jets;
  std::vector<phys::Jet> *DownJER_jets;
  
  // Jets obtained correcting up and down for the JES uncertainty
  std::vector<phys::Jet> *UpJES_jets;
  std::vector<phys::Jet> *DownJES_jets;

  // Jets obtained correcting up and down for the JES uncertainty (data distributions = no JER correction applied)
  std::vector<phys::Jet> *UpJESData_jets;
  std::vector<phys::Jet> *DownJESData_jets;

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

