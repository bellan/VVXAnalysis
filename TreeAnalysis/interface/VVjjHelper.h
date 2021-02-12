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
#include "VVXAnalysis/Commons/interface/Utilities.h"
#include <time.h>

using namespace std;
using namespace phys;
using namespace physmath;

class VVjjHelper{

public:

  VVjjHelper(Histogrammer *histopointer){
    histo_ = histopointer;
  }

  virtual ~VVjjHelper(){}
  
  static void test(int number);
    
  
 private:

  friend class WZAnalyzer;
  friend class ZZjjAnalyzer;

  // Getter functions
  unsigned int GetLeptonsNumber();
  unsigned int GetNeutrinosNumber();

  
  // Data memebers
  const float rangeVmass = 30.;
  Histogrammer *histo_;
    
  vector<Particle> neutrinos_;
  vector<Particle> leptons_;
  
  // Private member functions
  void printTime(float btime, float etime){
    cout << "\nExecution time: " << (int)((etime - btime)/3600) << " h " << (((int)(etime - btime)%3600)/60) << " m " << etime - btime - (int)((etime - btime)/3600)*3600 - (((int)(etime - btime)%3600)/60)*60 << " s." << endl;
  }
  
  void LeptonSearch(const vector<Particle> &genparticles, string eventkind, vector<Particle> &lepm, vector<Particle> &lepp, vector<Particle> &neutrino);
  void FindLeadingJets(const vector<Particle> *jetcollection, Particle &Jet0, Particle &Jet1, const vector<Particle> *particlecollection);
  void FindLeadingJets(const vector<Jet> *jetcollection, Particle &Jet0, Particle &Jet1);

  
  // Histogram functions
  void PlotParticle(const Particle &particle, string name, float weight);
  void PlotJets(const Particle &Jet0, const Particle &Jet1, string prename, float weight);

  template <class BOS>
  void PlotBoson(const BOS &particle, string name, float weight);
  
  template <class DiBOS>
  void PlotDiBoson(const DiBOS& particle, string name, float weight);

};

#endif
