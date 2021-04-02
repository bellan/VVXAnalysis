#include "Electron.h"
#include "Lepton.h"
#include "Jet.h"
#include "Boson.h"
#include "DiBoson.h"

namespace phys{
  typedef phys::DiBoson<phys::Lepton  , phys::Lepton>   ZZ4MU;
  typedef phys::DiBoson<phys::Electron, phys::Electron> ZZ4E;
  typedef phys::DiBoson<phys::Electron, phys::Lepton>   ZZ2E2MU;
  typedef phys::DiBoson<phys::Particle, phys::Particle> DiBosonParticle;
  typedef phys::DiBoson<phys::Lepton, phys::Lepton>     DiBosonLepton;
  
  typedef phys::Boson<phys::Particle>  BosonParticle;
  typedef phys::Boson<phys::Lepton>    BosonLepton;
  
  typedef std::pair<phys::Particle, phys::Particle>              pairParticle;
  typedef std::pair<phys::Boson<phys::Particle>, phys::Particle> pairBosonParticle;
}
