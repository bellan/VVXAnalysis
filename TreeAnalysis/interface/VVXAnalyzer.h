#ifndef VVXAnalyzer_h
#define VVXAnalyzer_h

/** \class VVXAnalyzer
 *  Concrete class for VVX analysis
 *
 *  $Date: 2013/03/15 13:37:31 $
 *  $Revision: 1.4 $
 *  \author R. Bellan - UNITO <riccardo.bellan@cern.ch>
 */


#include "EventAnalyzer.h"
#include "RegistrableAnalysis.h"
#include "VVXAnalysis/Commons/interface/Constants.h"

class VVXAnalyzer: public EventAnalyzer, RegistrableAnalysis<VVXAnalyzer>{

public:

 VVXAnalyzer(const std::string& region, const std::string& filename, const double& lumi = 1., const double& externalXSection = -1., bool doBasicPlots = false)
   : EventAnalyzer(*(new Selector<VVXAnalyzer>(*this)), 
		   region,
		   filename, lumi, externalXSection, doBasicPlots){
    theHistograms.profile(genCategory);
  }

  virtual ~VVXAnalyzer(){}

  virtual void analyze();
  
  virtual Int_t cut();

  void ZZplots(int id = -1);
  
 private:
  friend class Selector<VVXAnalyzer>;
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
};
#endif

