//-----------FUNCTION: definition of the two ZZ bosons from leptons-------

#include <utility>
#include "DataFormats/Candidate/interface/Candidate.h"
#include "VVXAnalysis/DataFormats/interface/Boson.h"
#include <tuple>


namespace zzw{

  typedef std::tuple<int,           // GenCategory
    phys::Boson<phys::Particle>,    // Z0 --> ll
    phys::Boson<phys::Particle>,    // Z1 --> ll
    phys::Boson<phys::Particle>,    // Z2 --> jj
    phys::Boson<phys::Particle>     //  W --> jj
    > GenTopology;

  
  std::pair<phys::Boson<phys::Particle>, phys::Boson<phys::Particle> > makeZBosonsFromLeptons(const std::vector<const reco::Candidate *>& lm, const std::vector<const reco::Candidate *>& lp, int leptonCode, float mZ);
  
  std::pair<phys::Boson<phys::Particle>, phys::Boson<phys::Particle> > makeZBosonsFromLeptons(const std::vector<phys::Particle>& lm, const std::vector<phys::Particle>& lp, int leptonCode, float mZ);

  int makeVBosonsFromIds(int j0Id, int j1Id);
  
  GenTopology getGenTopology(int signalDefinition,
			     const std::vector<phys::Particle> &theGenl, const std::vector<phys::Particle> &theGenj, 
			     const std::vector<phys::Particle> &theGenZ, const std::vector<phys::Particle> &theGenW);
}
