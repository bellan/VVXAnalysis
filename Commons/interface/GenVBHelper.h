#ifndef GenVBHelper_h
#define GenVBHelper_h

/** \class GenVBHelper
 *  Concrete class for VVX analysis
 *
 *  $Date: 2013/03/15 13:37:31 $
 *  $Revision: 1.4 $
 *  \author R. Bellan - UNITO <riccardo.bellan@cern.ch>
 */


#include "VVXAnalysis/Commons/interface/Constants.h"
#include "VVXAnalysis/DataFormats/interface/Boson.h"

class GenVBHelper{

public:

  GenVBHelper(){}

  virtual ~GenVBHelper(){}

  virtual void analyze(const std::vector<phys::Particle>& genParticles, const std::vector<phys::Boson<phys::Particle> >& genVBParticles);

  std::vector<phys::Boson<phys::Particle> >  ZtoChLep()     const {return ZtoChLep_;}
  std::vector<phys::Boson<phys::Particle> >  ZtoNeutrinos() const {return ZtoNeutrinos_;}
  std::vector<phys::Boson<phys::Particle> >  ZtoQ()         const {return ZtoQ_;}
  std::vector<phys::Boson<phys::Particle> >  WtoQ()         const {return WtoQ_;}
  std::vector<phys::Boson<phys::Particle> >  WtoLep()       const {return WtoLep_;}


  

 private:

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

  // New data memebers
  std::vector<phys::Boson<phys::Particle> > ZtoChLep_;
  std::vector<phys::Boson<phys::Particle> > ZtoNeutrinos_;
  std::vector<phys::Boson<phys::Particle> > ZtoQ_;
  std::vector<phys::Boson<phys::Particle> > WtoQ_;
  std::vector<phys::Boson<phys::Particle> > WtoLep_;
};
#endif

