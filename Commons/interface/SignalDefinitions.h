//-----------FUNCTION: definition of the two ZZ bosons from leptons-------



#ifndef VVXAnalysis_Commons_SignalDefinitions_H
#define VVXAnalysis_Commons_SignalDefinitions_H

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
    phys::Boson<phys::Particle>,    // Z0 --> ll
    phys::Boson<phys::Particle>,    // Z1 --> ll
    phys::Boson<phys::Particle>,    // Z2 --> jj
    phys::Boson<phys::Particle>     //  W --> jj
    > GenTopology;
  

 std::pair<phys::Boson<phys::Particle>, phys::Boson<phys::Particle> > makeZBosonsFromLeptons(const std::vector<phys::Particle>& lm, const std::vector<phys::Particle>& lp, int leptonCode, float mZ);
 
 GenTopology getGenTopology(int signalDefinition,
			    const std::vector<phys::Particle> &theGenl, const std::vector<phys::Particle> &theGenj, 
			    const std::vector<phys::Particle> &theGenZ, const std::vector<phys::Particle> &theGenW); 
}

namespace zz{
  
  std::tuple<bool, phys::Boson<phys::Particle>, phys::Boson<phys::Particle> > getZZ(const std::vector<phys::Boson<phys::Particle> >  &Zll);



  typedef std::tuple<int,           // SignCategory
    phys::Boson<phys::Particle>,    // Z0 --> ll
    phys::Boson<phys::Particle>,    // Z1 --> ll
    phys::Boson<phys::Particle>,    // Z2 --> ll
    phys::Boson<phys::Particle>,    // Z3 --> jj
    phys::Boson<phys::Particle>,    // W0 --> lv
    phys::Boson<phys::Particle>,    // W1 --> lv
    phys::Boson<phys::Particle>     // W2 --> jj
    > SignalTopology;


  SignalTopology getSignalTopology       (const std::vector<phys::Particle> &theGenl, std::vector<phys::Particle> &theGenj, std::vector<phys::Particle> &theGenjAK8);
  
  bool inTriggerPlateau(const std::vector<phys::Particle>& leptons);

  bool inTightFiducialRegion(const zz::SignalTopology &topology);

  bool inHiggsFiducialRegion(const zz::SignalTopology &topology);


}

namespace wz{
  std::tuple<bool, phys::Boson<phys::Particle>, phys::Boson<phys::Particle> > getWZ(const std::vector<phys::Particle>& lepMinus,
										    const std::vector<phys::Particle>& lepPlus,
										    const std::vector<phys::Particle>& neutrinos);
}



namespace vv{

  bool checkLeptonAcceptance(const phys::Particle& lepton, 
			     const double& pt_e, const double&eta_e, 
			     const double& pt_mu, const double& eta_mu);

  bool inLeptonAcceptance(const std::vector<phys::Particle>& leptons,
			  const double& pt_e, const double&eta_e, 
			  const double& pt_mu, const double& eta_mu);


  std::vector<phys::Boson<phys::Particle> > categorizeHadronicPartOftheEvent(std::vector<phys::Particle> &theGenj,
									     std::vector<phys::Particle> &theGenjAK8,
									     const std::vector<phys::Boson<phys::Particle> >& bosonsToLeptons, 
									     std::bitset<16>& topology);
  


}




#endif
