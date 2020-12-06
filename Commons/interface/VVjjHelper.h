#ifndef VVjjHelper_h
#define VVjjHelper_h

/** \class VVjjHelper
 *  Helper class for VVjj analysis
 *
 *  $Date: 2020/11/18 $
 *  $Revision: 0.5 $
 *  \author E. Racca - UNITO <eleonora.racca@cern.ch>
 */


#include "VVXAnalysis/Commons/interface/Constants.h"
#include "VVXAnalysis/DataFormats/interface/DiBoson.h"
#include "VVXAnalysis/DataFormats/interface/Particle.h"
#include "VVXAnalysis/DataFormats/interface/TypeDefs.h"

using namespace std;
using namespace phys;

class VVjjHelper{

public:

  VVjjHelper(){}

  virtual ~VVjjHelper(){}
  
  static void test();

  static bool FindDiBoson(vector<Particle> &genparticles, VVtype &VV, string eventtype);

  static bool FindLeadingJets(vector<Particle> &genjets, vector<Particle> &jets);
  

 private:

  // Private member functions
 

  // Data memebers
  const float range = 30.;
  
  int counter1;
  int counter2;
  int counter3;
  int counter4;
  
  void LeptonSearch(vector<Particle> &genparticles, string eventtype);
  VVtype BuildVV(string eventtype);
  unsigned int GetAllLeptonsNumber();
  unsigned int GetNeutrinosNumber();
  
  //vector<Particle> electrons_;
  //vector<Particle> muons_;
  vector<Particle> neutrinos_;
  vector<Particle> leptons_;
  
};
#endif
