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

#include "VVXAnalysis/TreeAnalysis/interface/Histogrammer.h"

using namespace std;
using namespace phys;

class VVjjHelper{

public:

  VVjjHelper(Histogrammer *histopointer){
    histo_ = histopointer;
  }

  virtual ~VVjjHelper(){}
  
  static void test();
    
  void LeptonSearch(const vector<Particle> &genparticles, string eventkind);
  void FindLeadingJets(vector<Particle> &jetcollection, vector<Particle> &particlecollection, Particle &Jet0, Particle &Jet1);
  DiBosonParticle BuildVV(string eventkind);
  unsigned int GetAllLeptonsNumber();
  unsigned int GetNeutrinosNumber();

 private:

  // Private member functions
 

  // Data memebers
  const float rangeVmass = 30.;
  Histogrammer *histo_;
    
  vector<Particle> neutrinos_;
  vector<Particle> leptons_;
  
};
#endif
