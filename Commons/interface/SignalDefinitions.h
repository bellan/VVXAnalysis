//-----------FUNCTION: definition of the two ZZ bosons from leptons-------



#ifndef VVXAnalysis_Commons_SignalDefinitions_H
#define VVXAnalysis_Commons_SignalDefinitions_H

#include "VVXAnalysis/DataFormats/interface/GenParticle.h"
#include "VVXAnalysis/DataFormats/interface/Boson.h"
#include "VVXAnalysis/DataFormats/interface/Lepton.h"
#include "VVXAnalysis/DataFormats/interface/Electron.h"
#include "VVXAnalysis/DataFormats/interface/GenStatusBit.h"
#include <tuple>
#include <bitset>

namespace vvx{

int makeVBosonsFromIds(int j0Id, int j1Id);


}
namespace zzw{

  typedef std::tuple<int,           // GenCategory
    phys::Boson<phys::GenParticle>,    // Z0 --> ll
    phys::Boson<phys::GenParticle>,    // Z1 --> ll
    phys::Boson<phys::GenParticle>,    // Z2 --> jj
    phys::Boson<phys::GenParticle>     //  W --> jj
    > GenTopology;
  

  std::pair<phys::Boson<phys::GenParticle>, phys::Boson<phys::GenParticle> > makeZBosonsFromLeptons(const std::vector<phys::GenParticle>& lm, const std::vector<phys::GenParticle>& lp, int leptonCode, float mZ);
 
 GenTopology getGenTopology(int signalDefinition,
			    const std::vector<phys::GenParticle> &theGenl, const std::vector<phys::GenParticle> &theGenj, 
			    const std::vector<phys::GenParticle> &theGenZ, const std::vector<phys::GenParticle> &theGenW); 
}

namespace zz{
  
  std::tuple<bool, phys::Boson<phys::GenParticle>, phys::Boson<phys::GenParticle> > getZZ(const std::vector<phys::Boson<phys::GenParticle> >  &Zll);



  typedef std::tuple<int,           // SignCategory
    phys::Boson<phys::GenParticle>,    // Z0 --> ll
    phys::Boson<phys::GenParticle>,    // Z1 --> ll
    phys::Boson<phys::GenParticle>,    // Z2 --> ll
    phys::Boson<phys::GenParticle>,    // Z3 --> jj
    phys::Boson<phys::GenParticle>,    // W0 --> lv
    phys::Boson<phys::GenParticle>,    // W1 --> lv
    phys::Boson<phys::GenParticle>     // W2 --> jj
    > SignalTopology;


  SignalTopology getSignalTopology       (const std::vector<phys::GenParticle> &theGenl, std::vector<phys::GenParticle> &theGenj, std::vector<phys::GenParticle> &theGenjAK8);
  
  bool inTriggerPlateau(const std::vector<phys::GenParticle>& leptons);

  bool inTightFiducialRegion(const zz::SignalTopology &topology);

  bool inHiggsFiducialRegion(const zz::SignalTopology &topology);


}

namespace wz{
  std::tuple<bool, phys::Boson<phys::GenParticle>, phys::Boson<phys::GenParticle> > getWZ(const std::vector<phys::GenParticle>& lepMinus,
										    const std::vector<phys::GenParticle>& lepPlus,
										    const std::vector<phys::GenParticle>& neutrinos);
}



namespace vv{

  bool checkLeptonAcceptance(const phys::GenParticle& lepton, 
			     const double& pt_e, const double&eta_e, 
			     const double& pt_mu, const double& eta_mu);

  bool inLeptonAcceptance(const std::vector<phys::GenParticle>& leptons,
			  const double& pt_e, const double&eta_e, 
			  const double& pt_mu, const double& eta_mu);


  std::vector<phys::Boson<phys::GenParticle> > categorizeHadronicPartOftheEvent(std::vector<phys::GenParticle> &theGenj,
									     std::vector<phys::GenParticle> &theGenjAK8,
									     const std::vector<phys::Boson<phys::GenParticle> >& bosonsToLeptons, 
									     std::bitset<16>& topology);
  


}




#endif
