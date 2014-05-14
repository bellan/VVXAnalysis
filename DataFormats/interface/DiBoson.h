#ifndef ZZWAnalysis_DataFormats_DiBoson_H
#define ZZWAnalysis_DataFormats_DiBoson_H

/** \class DiBoson
 *  No description available.
 *
 *  $Date: 2013/03/15 13:37:31 $
 *  $Revision: 1.3 $
 *  \author R. Bellan - UNITO <riccardo.bellan@cern.ch>
 */

#include "Particle.h"
#include "Boson.h"

namespace phys {
  
  template <typename P1, typename P2>
    class DiBoson: public Particle {
    
  public:
    /// Constructor
  DiBoson(): Particle()
      , regionWord_(0)
      , isBestCand_(false)
      , passFullSel_(false)
      {}

  DiBoson(const Boson<P1>& vb1, const Boson<P2>& vb2)
    : Particle(vb1.p4()+vb2.p4(),0,0)
      , regionWord_(0)
      , isBestCand_(false)
      , passFullSel_(false)
      {}

    template<typename T>
      Boson<T> daughter(int i) const{
      if(i == 0) return daughter0_;
      else if(i == 1) return daughter1_;
      else { std::cout << "*** DiBoson's daughter not found! ***" << " " << i << std::endl; abort();}
    }

    // Best candidate in the Control/Search region
    bool isBestCandidate() const {return isBestCand_;}

    // True if pass all requirements on di-boson quantities and on
    // its daughters and grand daughters and...
    bool passFullSelection() const {return passFullSel_;}

    // Type of search/control region
    short region() const {return regionWord_;}

  private:

    Boson<P1> daughter0_;
    Boson<P2> daughter1_;
    
    Short_t regionWord_;
    Bool_t isBestCand_;
    Bool_t passFullSel_;

    ClassDef(DiBoson, 1) //
  };
}
#endif
