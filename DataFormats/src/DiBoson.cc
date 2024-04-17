#include "VVXAnalysis/DataFormats/interface/Boson.h"
#include "VVXAnalysis/DataFormats/interface/DiBoson.h"
#include "VVXAnalysis/DataFormats/interface/Lepton.h"
#include "VVXAnalysis/DataFormats/interface/Photon.h"
#include "VVXAnalysis/DataFormats/interface/Jet.h"
#include "VVXAnalysis/DataFormats/interface/Electron.h"

template <class P1, class P2>
phys::DiBoson<P1,P2>::~DiBoson(){}

template<> phys::DiBoson<phys::Particle, phys::Particle >::~DiBoson(){}
template<> phys::DiBoson<phys::Lepton  , phys::Lepton   >::~DiBoson(){}
template<> phys::DiBoson<phys::Photon  , phys::Photon   >::~DiBoson(){}
template<> phys::DiBoson<phys::Jet     , phys::Jet      >::~DiBoson(){}
template<> phys::DiBoson<phys::Electron, phys::Electron >::~DiBoson(){}

template<> phys::DiBoson<phys::Lepton  , phys::Jet      >::~DiBoson(){}
template<> phys::DiBoson<phys::Jet     , phys::Lepton   >::~DiBoson(){}
template<> phys::DiBoson<phys::Lepton  , phys::Photon   >::~DiBoson(){}
template<> phys::DiBoson<phys::Photon  , phys::Lepton   >::~DiBoson(){}
template<> phys::DiBoson<phys::Photon  , phys::Jet      >::~DiBoson(){}
template<> phys::DiBoson<phys::Jet     , phys::Photon   >::~DiBoson(){}
