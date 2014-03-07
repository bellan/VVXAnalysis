//-----------FUNCTION: definition of the two ZZ bosons from leptons-------

#include <utility>
#include "Boson.h"
#include "DataFormats/Candidate/interface/Candidate.h"


std::pair<Boson*,Boson*> makeZbosonsFromLeptons(const std::vector<const reco::Candidate *>& lm, const std::vector<const reco::Candidate *>& lp, int leptonCode, float mZ);


