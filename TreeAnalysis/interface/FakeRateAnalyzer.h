#ifndef FakeRateAnalyzer_h
#define FakeRateAnalyzer_h

/** \class FakeRateAnalyzer
 *  Concrete class for VVX analysis
 *
 *  $Date: 2013/03/15 13:37:31 $
 *  $Revision: 1.4 $
 *  \author R. Bellan - UNITO <riccardo.bellan@cern.ch>
 */


#include "EventAnalyzer.h"
#include "RegistrableAnalysis.h"
#include "VVXAnalysis/Commons/interface/Constants.h"

class FakeRateAnalyzer: public EventAnalyzer, RegistrableAnalysis<FakeRateAnalyzer>{

public:

  //, const std::string& filename, const double& lumi = 1., const double& externalXSection = -1., bool doBasicPlots = false

 FakeRateAnalyzer(const AnalysisConfiguration& configuration)
   : EventAnalyzer(*(new Selector<FakeRateAnalyzer>(*this)), 
		   configuration)
    // , lepSF("/home/bellan/Workspace/WZZ/NtupleTestBed/VVXAnalysis/TreeAnalysis/fakeRates.root","/home/bellan/Workspace/WZZ/NtupleTestBed/VVXAnalysis/TreeAnalysis/fakeRates.root",
   //	    "/home/bellan/Workspace/WZZ/NtupleTestBed/VVXAnalysis/TreeAnalysis/fakeRates.root","/home/bellan/Workspace/WZZ/NtupleTestBed/VVXAnalysis/TreeAnalysis/fakeRates.root")

{
    //theHistograms.profile(genCategory);
  }

  virtual ~FakeRateAnalyzer(){}

  virtual void analyze();
  
  virtual Int_t cut();

  void ZZplots(int id = -1);

  virtual void begin();
  virtual void end( TFile &);
  virtual void addOptions();


 private:
  friend class Selector<FakeRateAnalyzer>;
  template< class PAR >
    bool ZBosonDefinition(phys::Boson<PAR> cand) const{
    bool checkCharge = cand.daughter(0).charge() + cand.daughter(1).charge() == 0;
    return checkCharge && fabs(cand.p4().M() - phys::ZMASS) < 30;
  }
  /* bool ZBosonDefinition(phys::Boson<phys::Lepton> cand) const{ */
  /*   return fabs(cand.p4().M() - ZMASS) < 15; */
  /* } */
  
  /* bool ZBosonDefinition(phys::Boson<phys::Electron> cand) const{ */

  /*   return fabs(cand.p4().M() - ZMASS) < 15; */
  /* } */

  std::vector<double> xbins_ele;
  std::vector<double> ybins_ele;
  std::vector<double> xbins_mu;
  std::vector<double> ybins_mu;


  template< class PAR >
    bool WBosonDefinition(phys::Boson<PAR> cand) {

    bool gooddaughters = false;
    if(fabs(cand.daughter(0).eta()) < 2.5 && cand.daughter(0).pt() > 30 &&
       cand.daughter(0).passPUID() && cand.daughter(0).passLooseJetID() &&
       fabs(cand.daughter(1).eta()) < 2.5 && cand.daughter(1).pt() > 30 &&
       cand.daughter(1).passPUID() && cand.daughter(1).passLooseJetID())
      gooddaughters = true;

    if(fabs(cand.p4().M() - phys::WMASS) < 150 && gooddaughters)
      return true;
    return false;

  }


  std::vector<std::string> eventsStr; 
  std::vector<std::string> eventsN_eee; 
  std::vector<std::string> eventsN_eem; 
  std::vector<std::string> eventsN_mme;
  std::vector<std::string> eventsN_mmm;

  std::vector<std::string> eventsD_eee; 
  std::vector<std::string> eventsD_eem; 
  std::vector<std::string> eventsD_mme;
  std::vector<std::string> eventsD_mmm;


};
#endif

