#include "VVXAnalysis/DataFormats/interface/Boson.h"
#include "VVXAnalysis/DataFormats/interface/Lepton.h"
#include "VVXAnalysis/DataFormats/interface/Photon.h"

template <class T>
phys::Boson<T>::~Boson(){}

template<> phys::Boson<phys::Particle>::~Boson(){}
template<> phys::Boson<phys::Lepton  >::~Boson(){}
template<> phys::Boson<phys::Photon  >::~Boson(){}
