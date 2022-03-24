#include "VVXAnalysis/DataFormats/interface/Particle.h"
#include "VVXAnalysis/DataFormats/interface/Lepton.h"
#include "VVXAnalysis/DataFormats/interface/Jet.h"
#include "VVXAnalysis/DataFormats/interface/Electron.h"
#include "VVXAnalysis/DataFormats/interface/Photon.h"
#include "VVXAnalysis/DataFormats/interface/Boson.h"
#include "VVXAnalysis/DataFormats/interface/DiBoson.h"
#include "VVXAnalysis/DataFormats/interface/GenEventWeights.h"
#include "VVXAnalysis/DataFormats/interface/MELA.h"
#include "VVXAnalysis/DataFormats/interface/RegionsCounter.h"
#include "VVXAnalysis/DataFormats/interface/RegionTypes.h"

#ifdef __CINT__

#pragma link C++ class  phys::GenEventWeights+;
#pragma link C++ class  phys::MELA+;
#pragma link C++ class  phys::RegionTypes+;
#pragma link C++ class  phys::RegionsCounter+;
#pragma link C++ class  std::map<phys::RegionTypes,Int_t>+;
#pragma link C++ class  phys::Jet::JetScores+;


#pragma link C++ class  phys::Particle+;
#pragma link C++ class  phys::Lepton+;
#pragma link C++ class  phys::Jet+;
#pragma link C++ class  phys::Electron+;
#pragma link C++ class  phys::Photon+;
#pragma link C++ class  phys::Proton+;
#pragma link C++ class  phys::Boson<phys::Particle>+;
#pragma link C++ class  phys::Boson<phys::Lepton>+;
#pragma link C++ class  phys::Boson<phys::Electron>+;
#pragma link C++ class  phys::Boson<phys::Jet>+;
#pragma link C++ class  phys::Boson<phys::Photon>+;
#pragma link C++ class  phys::DiBoson<phys::Particle, phys::Particle >+;
#pragma link C++ class  phys::DiBoson<phys::Lepton  , phys::Lepton >+;
#pragma link C++ class  phys::DiBoson<phys::Electron, phys::Lepton >+;
#pragma link C++ class  phys::DiBoson<phys::Lepton  , phys::Electron >+;
#pragma link C++ class  phys::DiBoson<phys::Electron, phys::Electron >+;

#pragma link C++ class  std::vector<phys::Particle>;
#pragma link C++ class  std::vector<phys::Lepton>;
#pragma link C++ class  std::vector<phys::Jet>;
#pragma link C++ class  std::vector<phys::Electron>;
#pragma link C++ class  std::vector<phys::Photon>;
#pragma link C++ class  std::vector<phys::Proton>;
#pragma link C++ class  std::vector<phys::Boson<phys::Particle> >;
#pragma link C++ class  std::vector<phys::Boson<phys::Lepton> >;
#pragma link C++ class  std::vector<phys::Boson<phys::Electron> >;
#pragma link C++ class  std::vector<phys::Boson<phys::Jet> >;
#pragma link C++ class  std::vector<phys::Boson<phys::Photon> >;
#pragma link C++ class  std::vector<phys::DiBoson<phys::Particle, phys::Particle > >;
#pragma link C++ class  std::vector<phys::DiBoson<phys::Lepton  , phys::Lepton > >;
#pragma link C++ class  std::vector<phys::DiBoson<phys::Electron, phys::Lepton > >;
#pragma link C++ class  std::vector<phys::DiBoson<phys::Lepton  , phys::Electron > >;
#pragma link C++ class  std::vector<phys::DiBoson<phys::Electron, phys::Electron > >;
#pragma link C++ class  std::pair<phys::Boson<phys::Lepton>, phys::Lepton>+;
#pragma link C++ class  std::vector<std::pair<phys::Boson<phys::Lepton>, phys::Lepton> >;


#endif
