#include "Electron.h"
#include "Lepton.h"
#include "Jet.h"
#include "Boson.h"
#include "DiBoson.h"

namespace phys{
  typedef phys::DiBoson<phys::Lepton  , phys::Lepton>   ZZ4MU;
  typedef phys::DiBoson<phys::Electron, phys::Electron> ZZ4E;
  typedef phys::DiBoson<phys::Electron, phys::Lepton>   ZZ2E2MU;
  typedef phys::Boson<phys::Particle>  BosonParticle;
}
