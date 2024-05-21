#include <VVXAnalysis/DataFormats/interface/Particle.h>
#include <VVXAnalysis/DataFormats/interface/Lepton.h>
#include <VVXAnalysis/DataFormats/interface/Photon.h>
#include <VVXAnalysis/DataFormats/interface/Boson.h>
#include <VVXAnalysis/DataFormats/interface/DiBoson.h>

phys::Boson<phys::Particle> dummyBosonParticle;
phys::Boson<phys::Lepton>   dummyBosonLepton;
phys::Boson<phys::Photon>   dummyBosonPhoton;

phys::DiBoson<phys::Particle, phys::Particle> dummyDiBosonParticleParticle;
phys::DiBoson<phys::Lepton  , phys::Lepton>   dummyDiBosonLeptonLepton;
phys::DiBoson<phys::Photon  , phys::Photon>   dummyDiBosonPhotonPhoton;
