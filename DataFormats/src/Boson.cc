#include "VVXAnalysis/DataFormats/interface/Boson.h"
#include "VVXAnalysis/DataFormats/interface/Lepton.h"
#include "VVXAnalysis/DataFormats/interface/Photon.h"
#include "VVXAnalysis/DataFormats/interface/Jet.h"
#include "VVXAnalysis/DataFormats/interface/Electron.h"

template <class T>
phys::Boson<T>::~Boson(){}

template<> phys::Boson<phys::Particle>::~Boson(){}
template<> phys::Boson<phys::Lepton  >::~Boson(){}
template<> phys::Boson<phys::Photon  >::~Boson(){}
template<> phys::Boson<phys::Jet     >::~Boson(){}
template<> phys::Boson<phys::Electron>::~Boson(){}
