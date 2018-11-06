
#include "../interface/Particle.h"
#include "../interface/Lepton.h"
#include "../interface/Jet.h"
#include "../interface/Electron.h"
#include "../interface/Boson.h"
#include "../interface/DiBoson.h"
#include "../interface/GenEventInfo.h"
#include "../interface/MELA.h"

#ifdef __CINT__

#pragma link C++ class  phys::GenEventInfo+;
#pragma link C++ class  phys::MELA+;

#pragma link C++ class  phys::Particle+;
#pragma link C++ class  phys::Lepton+;
#pragma link C++ class  phys::Jet+;
#pragma link C++ class  phys::Electron+;
#pragma link C++ class  phys::Boson<phys::Particle>+;
#pragma link C++ class  phys::Boson<phys::Lepton>+;
#pragma link C++ class  phys::Boson<phys::Electron>+;
#pragma link C++ class  phys::Boson<phys::Jet>+;
#pragma link C++ class  phys::DiBoson<phys::Particle, phys::Particle >+;
#pragma link C++ class  phys::DiBoson<phys::Lepton  , phys::Lepton >+;
#pragma link C++ class  phys::DiBoson<phys::Electron, phys::Lepton >+;
#pragma link C++ class  phys::DiBoson<phys::Lepton  , phys::Electron >+;
#pragma link C++ class  phys::DiBoson<phys::Electron, phys::Electron >+;

#pragma link C++ class  std::vector<phys::Particle>;
#pragma link C++ class  std::vector<phys::Lepton>;
#pragma link C++ class  std::vector<phys::Jet>;
#pragma link C++ class  std::vector<phys::Electron>;
#pragma link C++ class  std::vector<phys::Boson<phys::Particle> >;
#pragma link C++ class  std::vector<phys::Boson<phys::Lepton> >;
#pragma link C++ class  std::vector<phys::Boson<phys::Electron> >;
#pragma link C++ class  std::vector<phys::Boson<phys::Jet> >;
#pragma link C++ class  std::vector<phys::DiBoson<phys::Particle, phys::Particle > >;
#pragma link C++ class  std::vector<phys::DiBoson<phys::Lepton  , phys::Lepton > >;
#pragma link C++ class  std::vector<phys::DiBoson<phys::Electron, phys::Lepton > >;
#pragma link C++ class  std::vector<phys::DiBoson<phys::Lepton  , phys::Electron > >;
#pragma link C++ class  std::vector<phys::DiBoson<phys::Electron, phys::Electron > >;
#pragma link C++ class  std::pair<phys::Boson<phys::Lepton>, phys::Lepton>+;
#pragma link C++ class  std::vector<std::pair<phys::Boson<phys::Lepton>, phys::Lepton> >;


#endif
