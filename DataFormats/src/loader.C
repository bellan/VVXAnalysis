#include "../interface/Particle.h"
#include "../interface/GenStatusBit.h"
// #include "../DataFormats/interface/Lepton.h"
// #include "../DataFormats/interface/Jet.h"
// #include "../DataFormats/interface/Electron.h"
#include "../interface/Boson.h"
// #include "../DataFormats/interface/DiBoson.h"
#include "../interface/Proton.h"
#include "../interface/ProtonPair.h"

#ifdef __CINT__

#pragma link C++ class  phys::Proton+;
#pragma link C++ class  phys::ProtonPair+;
#pragma link C++ class  phys::Particle+;
// #pragma link C++ class  phys::Lepton+;
// #pragma link C++ class  phys::Jet+;
// #pragma link C++ class  phys::Electron+;
#pragma link C++ class  phys::Boson<phys::Particle>+;
// #pragma link C++ class  phys::Boson<phys::Lepton>+;
// #pragma link C++ class  phys::Boson<phys::Electron>+;
// #pragma link C++ class  phys::Boson<phys::Jet>+;
// #pragma link C++ class  phys::DiBoson<phys::Lepton  , phys::Lepton >+;
// #pragma link C++ class  phys::DiBoson<phys::Electron, phys::Lepton >+;
// #pragma link C++ class  phys::DiBoson<phys::Lepton  , phys::Electron >+;
// #pragma link C++ class  phys::DiBoson<phys::Electron, phys::Electron >+;

#pragma link C++ enum GenStatusBits;

#pragma link C++ class  std::vector<phys::Proton>;
#pragma link C++ class  std::vector<phys::ProtonPair>;
#pragma link C++ class  std::vector<phys::Particle>;
// #pragma link C++ class  std::vector<phys::Lepton>;
// #pragma link C++ class  std::vector<phys::Jet>;
// #pragma link C++ class  std::vector<phys::Electron>;
#pragma link C++ class  std::vector<phys::Boson<phys::Particle> >;
// #pragma link C++ class  std::vector<phys::Boson<phys::Lepton> >;
// #pragma link C++ class  std::vector<phys::Boson<phys::Electron> >;
// #pragma link C++ class  std::vector<phys::Boson<phys::Jet> >;
// #pragma link C++ class  std::vector<phys::DiBoson<phys::Lepton  , phys::Lepton > >;
// #pragma link C++ class  std::vector<phys::DiBoson<phys::Electron, phys::Lepton > >;
// #pragma link C++ class  std::vector<phys::DiBoson<phys::Lepton  , phys::Electron > >;
// #pragma link C++ class  std::vector<phys::DiBoson<phys::Electron, phys::Electron > >;

#endif
