
#include "../interface/Particle.h"
#include "../interface/Lepton.h"
#include "../interface/Jet.h"
#include "../interface/Electron.h"
#include "../interface/Boson.h"

#ifdef __CINT__

#pragma link C++ class phys::Particle+;
#pragma link C++ class phys::Lepton+;
#pragma link C++ class phys::Jet+;
#pragma link C++ class phys::Electron+;
#pragma link C++ class phys::Boson<phys::Particle>+;
#pragma link C++ class phys::Boson<phys::Lepton>+;
#pragma link C++ class phys::Boson<phys::Electron>+;
#pragma link C++ class phys::Boson<phys::Jet>+;

#pragma link C++ class  std::vector<std::vector<phys::Particle> >;
#pragma link C++ class  std::vector<std::vector<phys::Lepton> >;
#pragma link C++ class  std::vector<std::vector<phys::Jet> >;
#pragma link C++ class  std::vector<std::vector<phys::Electron> >;
#pragma link C++ class  std::vector<std::vector<phys::Boson<phys::Particle> > >;
#pragma link C++ class  std::vector<std::vector<phys::Boson<phys::Lepton> > >;
#pragma link C++ class  std::vector<std::vector<phys::Boson<phys::Electron> > >;
#pragma link C++ class  std::vector<std::vector<phys::Boson<phys::Jet> > >;

#endif
