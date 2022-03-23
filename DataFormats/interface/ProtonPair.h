#ifndef ZZWAnalysis_DataFormats_ProtonPair_H
#define ZZWAnalysis_DataFormats_ProtonPair_H

/** \class ProtonPair
 *
 *  $Date: 2022/03/23  $
 *  $Revision: ??? $
 *  \author G. Marozzo - UNITO <giovanni.marozzo@edu.unito.it>
 */

#include "Proton.h"

namespace phys {
  
  template <typename P1, typename P2>
    class ProtonPair: {
    
    friend class ::TreePlanter;
  public:
    /// Constructor
    ProtonPair():
      , passFullSel_(false)
      {}
      
  ProtonPair(const Proton<P1>& p1, const Proton<P2>& p2):
      , daughter0_(p1)
      , daughter1_(p2)
      , passFullSel_(false)
      {}
    

    Proton<P1> first()  const {return daughter0_;}
    Proton<P2> second() const {return daughter1_;}

    Proton<P1> *firstPtr()  {return &daughter0_;}
    Proton<P2> *secondPtr() {return &daughter1_;}


    bool passFullSelection() const {return passFullSel_;}
    void setPassFullSel(Bool_t  fs) {passFullSel_ = fs;}

    //methods to extract the expected mass and rapidity value of the central system
    Double_t mpp() const {return (daughter0_.Ebeam()+daughter1_.Ebeam())*daughter0_.xi()*daughter1_.xi();}
    Double_t ypp() const {return -0.5*log(daughter0_.xi()/daughter1_.xi());}
    

  private:
  
    Proton daughter0_;
    Proton daughter1_;
    
    Bool_t  passFullSel_;

    ClassDef(ProtonPair, 1) //
  }
}
#endif
