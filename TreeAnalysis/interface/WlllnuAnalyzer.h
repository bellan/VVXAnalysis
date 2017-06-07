#ifndef WlllnuAnalyzer_h
#define WlllnuAnalyzer_h

/** \class WlllnuAnalyzer
 *  Concrete class for Wlllnu analysis
 *
 *  $Date: 2013/03/15 13:37:31 $
 *  $Revision: 1.4 $
 *  \author A. Corrado - UNITO <arianna.corrado@cern.ch>
 */


#include "EventAnalyzer.h"
#include "RegistrableAnalysis.h"
#include "VVXAnalysis/Commons/interface/Constants.h"
#include "VVXAnalysis/Commons/interface/LeptonScaleFactors.h"
#include "VVXAnalysis/DataFormats/interface/Boson.h"
#include "VVXAnalysis/DataFormats/interface/Particle.h"
#include "VVXAnalysis/Commons/interface/AriEle.h"

class WlllnuAnalyzer: public EventAnalyzer, RegistrableAnalysis<WlllnuAnalyzer>{

public:

  //, const std::string& filename, const double& lumi = 1., const double& externalXSection = -1., bool doBasicPlots = false
  
 WlllnuAnalyzer(const AnalysisConfiguration& configuration)
   : EventAnalyzer(*(new Selector<WlllnuAnalyzer>(*this)), 
		   configuration){
    //theHistograms.profile(genCategory);
  }

  virtual ~WlllnuAnalyzer(){}

  virtual void begin();
  
  virtual void analyze();
  
  virtual Int_t cut();
  
  
 
 private:
  Int_t nevents;
  
  friend class Selector<WlllnuAnalyzer>;
  template< class PAR >
    bool ZBosonDefinition(phys::Boson<PAR> cand) const{
    bool checkCharge = cand.daughter(0).charge() + cand.daughter(1).charge() == 0;
    return checkCharge && fabs(cand.p4().M() - phys::ZMASS) < 30;
  }

  
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

