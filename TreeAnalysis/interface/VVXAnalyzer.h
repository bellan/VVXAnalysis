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

  //, const std::string& filename, const double& lumi = 1., const double& externalXSection = -1., bool doBasicPlots = false

 VVXAnalyzer(const AnalysisConfiguration& configuration)
   : EventAnalyzer(*(new Selector<VVXAnalyzer>(*this)), 
		   configuration){
    //theHistograms.profile(genCategory);
  }

  virtual ~VVXAnalyzer(){}

  virtual void analyze();
  
  virtual Int_t cut();


 private:
  friend class Selector<VVXAnalyzer>;
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

  std::vector<phys::Particle> removeOverlaps(const std::vector<phys::Particle> &collectionX, const std::vector<phys::Boson<phys::Particle> >& collectionVB);

  
  template<class T1, class T2> 
  std::vector<phys::Boson<phys::Particle> > getVtoX(const std::vector<phys::Particle> & collectionX1,
						    const std::vector<phys::Particle> & collectionX2,
						    T1 idcondition, T2 masswindow, const double& referenceMass);

  
  class ZDaughtersIdCondition { 
  public: 
    // Comparator function 
    bool operator()(phys::Particle a, 
                    phys::Particle b){ 

      if(a.id() + b.id() == 0) return true; 
      else return false; 
    } 
  }; 

  class  WqqDaughtersIdCondition{ 
  public: 
    // Comparator function 
    bool operator()(phys::Particle a, 
                    phys::Particle b){ 

      if(abs((a.id() + b.id()))%2 == 1) return true; 
      else return false; 
    } 
  };

  class  WlnDaughtersIdCondition{ 
  public: 
    // Comparator function 
    bool operator()(phys::Particle l, 
                    phys::Particle n){ 

      if(abs(l.id() + n.id()) == 1 && abs(n.id()) > abs(l.id())) return true; 
      else return false; 
    } 
  };

  
 class ZMassWindow { 
  public: 
    // Comparator function 
    bool operator()(phys::Particle a){ 

      if(a.mass() <= 120 && a.mass() >= 60) return true; 
      else return false; 
    } 
  }; 

  class WMassWindow { 
  public: 
    // Comparator function 
    bool operator()(phys::Particle a){ 

      if(a.mass() <= 110 && a.mass() >= 50) return true; 
      else return false; 
    } 
  }; 

  


};
#endif

