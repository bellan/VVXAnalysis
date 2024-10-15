#ifndef VVXAnalysis_DataFormats_RecoParticle_H
#define VVXAnalysis_DataFormats_RecoParticle_H

/** \class Particle
 *  Base class for reconstructed particles
 *
 *  $Date: 2024/10/15 12:12:30 $
 *  $Revision: 1.0 $
 *  \author A. Mecca - Torino <alberto.mecca@cern.ch>
 */

#include "Particle.h"

namespace phys {
  class RecoParticle : public Particle {

    friend class ::TreePlanter;

  public:
    RecoParticle(const TLorentzVector& pin = TLorentzVector(0.,0.,0.,0.), float q = 0, int id = 0)
      : Particle(pin, q, id)
      , efficiencySF_(1.)
      , efficiencySFUnc_(0.)
      , fakeRateSF_(1.)
      , fakeRateSFUnc_(0.)
      {}

    virtual ~RecoParticle();

    virtual Float_t efficiencySF()   const {return efficiencySF_;}
    virtual Float_t efficiencySFUnc()const {return efficiencySFUnc_;}
    virtual Float_t fakeRateSF()     const {return fakeRateSF_;}
    virtual Float_t fakeRateSFUnc()  const {return fakeRateSFUnc_;} 
    virtual Float_t fakeRateSFVar()  const {return fakeRateSFUnc()*fakeRateSFUnc();}
    virtual Float_t energyScaleUp()  const {return scale_total_up_;}
    virtual Float_t energyScaleDn()  const {return scale_total_dn_;}
    virtual Float_t energySigmaUp()  const {return sigma_total_up_;}
    virtual Float_t energySigmaDn()  const {return sigma_total_dn_;}

    Bool_t passFullSel() const {return true;}
    void setEfficenySFUnc(float effSfUnc ) {efficiencySFUnc_ = effSfUnc;}

  protected:
    Float_t efficiencySF_, efficiencySFUnc_;
    Float_t fakeRateSF_, fakeRateSFUnc_;
    Float_t scale_total_up_, scale_total_dn_; // Momentum/energy scale uncertainty
    Float_t sigma_total_up_, sigma_total_dn_; // Resolution uncertainty
    
    ClassDef(RecoParticle, 1)
  };
}
#endif
