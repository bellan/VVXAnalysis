#ifndef WZZAnalyzer_h
#define WZZAnalyzer_h

/** \class WZZAnalyzer
 *  Concrete class for WZZ analysis
 *
 *  $Date: 2013/03/15 13:37:31 $
 *  $Revision: 1.6 $
 *  \author R. Bellan - UNITO <riccardo.bellan@cern.ch>
 */


#include "EventAnalyzer.h"
#include "RegistrableAnalysis.h"
#include "VVXAnalysis/Commons/interface/Constants.h"

#include <TString.h>

class WZZAnalyzer: public EventAnalyzer, RegistrableAnalysis<WZZAnalyzer>{

public:

  //, const std::string& filename, const double& lumi = 1., const double& externalXSection = -1., bool doBasicPlots = false

 WZZAnalyzer(const AnalysisConfiguration& configuration)
   : EventAnalyzer(*(new Selector<WZZAnalyzer>(*this)), 
		   configuration){
    //theHistograms.profile(genCategory);
  }

  virtual ~WZZAnalyzer(){}

  virtual void analyze();

  virtual void genAnalyze();
  
  virtual Int_t signalCostraint();

  virtual Bool_t cut(Int_t, phys::Boson<phys::Jet>);

  virtual void Reconstruct(phys::Boson<phys::Jet>*);

  virtual void CompatibilityTest(phys::Boson<phys::Jet>, phys::Boson<phys::Particle>, std::string, std::string);

  virtual void printHistos(uint, std::string, phys::Boson<phys::Jet>);

 private:
  friend class Selector<WZZAnalyzer>;
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

