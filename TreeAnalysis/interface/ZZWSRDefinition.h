#ifndef ZZWSRDefinition_h
#define ZZWSRDefinition_h

/** \class ZZWSRDefinition
 *  Concrete class for QGC analysis
 *
 *  $Date: 2013/03/15 13:37:31 $
 *  $Revision: 1.4 $
 *  \author R. Bellan - UNITO <riccardo.bellan@cern.ch>
 */


#include "EventAnalyzer.h"
#include "RegistrableAnalysis.h"
#include "VVXAnalysis/Commons/interface/Constants.h"

class ZZWSRDefinition: public EventAnalyzer, RegistrableAnalysis<ZZWSRDefinition>{

public:

 ZZWSRDefinition(const AnalysisConfiguration& configuration)
   : EventAnalyzer(*(new Selector<ZZWSRDefinition>(*this)),
		   configuration){}

  virtual ~ZZWSRDefinition(){}

  virtual void analyze();

  virtual Int_t cut();

 private:

  friend class Selector<ZZWSRDefinition>; 
  template< class PAR >
    bool ZBosonDefinition(phys::Boson<PAR> cand) const {
    return fabs(cand.p4().M() - phys::ZMASS) < 20;
  }
  template< class PAR >
    bool WBosonDefinition(phys::Boson<PAR> cand) const {
    bool massRange = fabs(cand.p4().M() - phys::WMASS) < 40;
    bool jetPt = (cand.daughter(0).pt() > 25 && cand.daughter(1).pt() > 40) || (cand.daughter(0).pt() > 40 && cand.daughter(1).pt() > 25); 
    bool jetsID = cand.daughter(0).passPUID() && cand.daughter(0).passLooseJetID() && cand.daughter(1).passPUID() && cand.daughter(1).passLooseJetID();
    return massRange && jetPt && jetsID;
    //return jetPt && jetsID;
  }

  phys::Boson<phys::Lepton> myZ0;
  phys::Boson<phys::Lepton> myZ1;
  phys::Boson<phys::Jet> myW;


};
#endif

