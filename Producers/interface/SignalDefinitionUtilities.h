//-----------FUNCTION: definition of the two ZZ bosons from leptons-------

#include <utility>
#include "DataFormats/Candidate/interface/Candidate.h"
#include "VVXAnalysis/DataFormats/interface/Boson.h"

std::pair<phys::Boson<phys::Particle>, phys::Boson<phys::Particle> > makeZbosonsFromLeptons(const std::vector<const reco::Candidate *>& lm, const std::vector<const reco::Candidate *>& lp, int leptonCode, float mZ);


