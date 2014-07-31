//-----------FUNCTION: definition of the two ZZ bosons from leptons-------

#ifndef VVXAnalysis_Commons_SignalDefinitions_H
#define VVXAnalysis_Commons_SignalDefinitions_H

#include "VVXAnalysis/DataFormats/interface/Boson.h"
#include "VVXAnalysis/DataFormats/interface/Lepton.h"
#include "VVXAnalysis/DataFormats/interface/Electron.h"
#include <ZZAnalysis/AnalysisStep/interface/bitops.h>
#include <tuple>

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
  std::tuple<bool, phys::Boson<phys::Lepton>, phys::Boson<phys::Lepton> >
    zz4l(const std::vector<phys::Boson<phys::Lepton> >  &Zmm,
	 const std::vector<phys::Boson<phys::Electron> > &Zee);
  
  std::tuple<bool, phys::Boson<phys::Particle>, phys::Boson<phys::Particle> > getZZ(const std::vector<phys::Boson<phys::Particle> >  &Zll);



  typedef std::tuple<int,           // SignCategory
    phys::Boson<phys::Particle>,    // Z0 --> ll
    phys::Boson<phys::Particle>,    // Z1 --> ll
    phys::Boson<phys::Particle>,    // Z2 --> ll
    phys::Boson<phys::Particle>,    // Z3 --> jj
    phys::Boson<phys::Particle>,    // W0 --> lv
    phys::Boson<phys::Particle>     // W1 --> jj
    > SignalTopology;


  SignalTopology getSignalTopology( const std::vector<phys::Particle> &theGenl, const std::vector<phys::Particle> &theGenj);			    
}



#endif
