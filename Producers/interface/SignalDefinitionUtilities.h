//-----------FUNCTION: definition of the two ZZ bosons from leptons-------

#include <utility>
#include "DataFormats/Candidate/interface/Candidate.h"
#include "VVXAnalysis/DataFormats/interface/Boson.h"

std::pair<phys::Boson<phys::Particle>, phys::Boson<phys::Particle> > makeZBosonsFromLeptons(const std::vector<const reco::Candidate *>& lm, const std::vector<const reco::Candidate *>& lp, int leptonCode, float mZ);

std::pair<phys::Boson<phys::Particle>, phys::Boson<phys::Particle> > makeZBosonsFromLeptons(const std::vector<phys::Particle>& lm, const std::vector<phys::Particle>& lp, int leptonCode, float mZ);

phys::Particle convert(const reco::Candidate&);

int makeVBosonsFromIds(int j0Id, int j1Id);
