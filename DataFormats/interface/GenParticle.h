#ifndef VVXAnalysis_DataFormats_GenParticle_H
#define VVXAnalysis_DataFormats_GenParticle_H

/** \class Particle
 *  Base class for reconstructed particles
 *
 *  $Date: 2024/10/15 12:12:30 $
 *  $Revision: 1.0 $
 *  \author A. Mecca - Torino <alberto.mecca@cern.ch>
 */

#include "Particle.h"

namespace phys {
  class GenParticle : public Particle {

    friend class ::TreePlanter;

  public:
    GenParticle(const TLorentzVector& mom = TLorentzVector(0.,0.,0.,0.), float q = 0, int id = 0, std::bitset<15> flags = 0)
      : Particle(mom, q, id)
      , genStatusFlags_(flags)
      , motherId_(-99)
      {}

    virtual ~GenParticle();

    // Gen info
    std::bitset<15> genStatusFlags() const { return genStatusFlags_; }
    void setGenStatusBit(GenStatusBit bit, int val = 1) {genStatusFlags_.set(bit, val);}

    void setMotherId(int pid) {motherId_ = pid;}
    int motherId() const {return motherId_;}

    void printStatusBits() const {
      std::cout << "isPrompt\t"  << genStatusFlags().test(phys::GenStatusBit::isPrompt) << std::endl
		<< "isDecayedLeptonHadron\t" << genStatusFlags().test(phys::GenStatusBit::isDecayedLeptonHadron)<< std::endl
		<< "isTauDecayProduct\t" << genStatusFlags().test(phys::GenStatusBit::isTauDecayProduct)<< std::endl
		<< "isPromptTauDecayProduct\t" << genStatusFlags().test(phys::GenStatusBit::isPromptTauDecayProduct)<< std::endl
		<< "isDirectTauDecayProduct\t" << genStatusFlags().test(phys::GenStatusBit::isDirectTauDecayProduct)<< std::endl
		<< "isDirectPromptTauDecayProduct\t" << genStatusFlags().test(phys::GenStatusBit::isDirectPromptTauDecayProduct)<< std::endl
		<< "isDirectHadronDecayProduct\t" << genStatusFlags().test(phys::GenStatusBit::isDirectHadronDecayProduct)<< std::endl
		<< "isHardProcess\t" << genStatusFlags().test(phys::GenStatusBit::isHardProcess)<< std::endl
		<< "fromHardProcess\t" << genStatusFlags().test(phys::GenStatusBit::fromHardProcess)<< std::endl
		<< "isHardProcessTauDecayProduct\t" << genStatusFlags().test(phys::GenStatusBit::isHardProcessTauDecayProduct)<< std::endl
		<< "isDirectHardProcessTauDecayProduct\t" << genStatusFlags().test(phys::GenStatusBit::isDirectHardProcessTauDecayProduct)<< std::endl
		<< "fromHardProcessBeforeFSR\t" << genStatusFlags().test(phys::GenStatusBit::fromHardProcessBeforeFSR)<< std::endl
		<< "isFirstCopy\t" << genStatusFlags().test(phys::GenStatusBit::isFirstCopy)<< std::endl
		<< "isLastCopy\t" << genStatusFlags().test(phys::GenStatusBit::isLastCopy)<< std::endl
		<< "isLastCopyBeforeFSR\t" << genStatusFlags().test(phys::GenStatusBit::isLastCopyBeforeFSR)<< std::endl;
    }

  protected:
    std::bitset<15> genStatusFlags_;
    Int_t motherId_;
    
    ClassDef(GenParticle, 1)
  };
}
#endif
