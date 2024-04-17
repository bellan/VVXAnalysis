#include <VVXAnalysis/DataFormats/interface/Particle.h>
#include <VVXAnalysis/DataFormats/interface/Lepton.h>
#include <VVXAnalysis/DataFormats/interface/Photon.h>
#include <VVXAnalysis/DataFormats/interface/Jet.h>
#include <VVXAnalysis/DataFormats/interface/Electron.h>
#include <VVXAnalysis/DataFormats/interface/Boson.h>
#include <VVXAnalysis/DataFormats/interface/DiBoson.h>

phys::Boson<phys::Particle> dummyBosonParticle;
phys::Boson<phys::Lepton>   dummyBosonLepton;
phys::Boson<phys::Photon>   dummyBosonPhoton;
phys::Boson<phys::Jet>      dummyBosonJet;
phys::Boson<phys::Electron> dummyBosonElectron;

phys::DiBoson<phys::Particle, phys::Particle> dummyDiBosonParticleParticle;
phys::DiBoson<phys::Lepton  , phys::Lepton>   dummyDiBosonLeptonLepton;
phys::DiBoson<phys::Photon  , phys::Photon>   dummyDiBosonPhotonPhoton;
phys::DiBoson<phys::Jet     , phys::Jet>      dummyDiBosonJetJet;
phys::DiBoson<phys::Electron, phys::Electron> dummyDiBosonElectronElectron;

phys::DiBoson<phys::Lepton  , phys::Jet     > dummyDiBosonLeptonJet;
phys::DiBoson<phys::Jet     , phys::Lepton  > dummyDiBosonJetLepton;
phys::DiBoson<phys::Lepton  , phys::Photon  > dummyDiBosonLeptonPhoton;
phys::DiBoson<phys::Photon  , phys::Lepton  > dummyDiBosonPhotonLepton;
phys::DiBoson<phys::Photon  , phys::Jet     > dummyDiBosonPhotonJet;
phys::DiBoson<phys::Jet     , phys::Photon  > dummyDiBosonJetPhoton;
