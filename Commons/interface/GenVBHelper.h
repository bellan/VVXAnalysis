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
#include "VVXAnalysis/DataFormats/interface/GenParticle.h"
#include "VVXAnalysis/DataFormats/interface/Boson.h"

class GenVBHelper{

public:

  GenVBHelper(){}

  virtual ~GenVBHelper(){}

  virtual void analyze(const std::vector<phys::GenParticle>& genGenParticles, const std::vector<phys::Boson<phys::GenParticle> >& genVBGenParticles);

  std::vector<phys::Boson<phys::GenParticle> >  ZtoChLep()     const {return ZtoChLep_;}
  std::vector<phys::Boson<phys::GenParticle> >  ZtoNeutrinos() const {return ZtoNeutrinos_;}
  std::vector<phys::Boson<phys::GenParticle> >  ZtoQ()         const {return ZtoQ_;}
  std::vector<phys::Boson<phys::GenParticle> >  WtoQ()         const {return WtoQ_;}
  std::vector<phys::Boson<phys::GenParticle> >  WtoLep()       const {return WtoLep_;}


  

 private:

  std::vector<phys::GenParticle> removeOverlaps(const std::vector<phys::GenParticle> &collectionX, const std::vector<phys::Boson<phys::GenParticle> >& collectionVB);

  
  template<class T1, class T2> 
  std::vector<phys::Boson<phys::GenParticle> > getVtoX(const std::vector<phys::GenParticle> & collectionX1,
						    const std::vector<phys::GenParticle> & collectionX2,
						    T1 idcondition, T2 masswindow, const double& referenceMass);

  
  class ZDaughtersIdCondition { 
  public: 
    // Comparator function 
    bool operator()(phys::GenParticle a, 
                    phys::GenParticle b){ 

      if(a.id() + b.id() == 0) return true; 
      else return false; 
    } 
  }; 

  class  WqqDaughtersIdCondition{ 
  public: 
    // Comparator function 
    bool operator()(phys::GenParticle a, 
                    phys::GenParticle b){ 

      if(abs((a.id() + b.id()))%2 == 1) return true; 
      else return false; 
    } 
  };

  class  WlnDaughtersIdCondition{ 
  public: 
    // Comparator function 
    bool operator()(phys::GenParticle l, 
                    phys::GenParticle n){ 

      if(abs(l.id() + n.id()) == 1 && abs(n.id()) > abs(l.id())) return true; 
      else return false; 
    } 
  };

  
 class ZMassWindow { 
  public: 
    // Comparator function 
    bool operator()(phys::GenParticle a){ 

      if(a.mass() <= 120 && a.mass() >= 60) return true; 
      else return false; 
    } 
  }; 

  class WMassWindow { 
  public: 
    // Comparator function 
    bool operator()(phys::GenParticle a){ 

      if(a.mass() <= 110 && a.mass() >= 50) return true; 
      else return false; 
    } 
  }; 

  // New data memebers
  std::vector<phys::Boson<phys::GenParticle> > ZtoChLep_;
  std::vector<phys::Boson<phys::GenParticle> > ZtoNeutrinos_;
  std::vector<phys::Boson<phys::GenParticle> > ZtoQ_;
  std::vector<phys::Boson<phys::GenParticle> > WtoQ_;
  std::vector<phys::Boson<phys::GenParticle> > WtoLep_;
};
#endif

