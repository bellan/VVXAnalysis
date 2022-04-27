#ifndef PPZZAnalyzer_h
#define PPZZAnalyzer_h

/** \class PPZZAnalyzer
 *  Class for analysis of Z boson pair production in diffraction events, with leading protons tagged by PPS
 *  $Date: 2022/4/20 $
 *  \author Giovanni Marozzo giovanni.marozzo@edu.unito.it
 */


#include "EventAnalyzer.h"
#include "RegistrableAnalysis.h"
#include "VVXAnalysis/Commons/interface/Constants.h"
#include "VVXAnalysis/DataFormats/interface/ProtonPair.h"

class PPZZAnalyzer: public EventAnalyzer, RegistrableAnalysis<PPZZAnalyzer>{

public:

  PPZZAnalyzer(const AnalysisConfiguration& configuration)
   : EventAnalyzer(*(new Selector<PPZZAnalyzer>(*this)), 
		   configuration){
   }

  virtual ~PPZZAnalyzer(){}

  virtual void analyze();
  
  virtual Int_t cut();


 private:
  friend class Selector<PPZZAnalyzer>;
  
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
