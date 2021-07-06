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

#include "TSpline.h"
#include "TFile.h"

using namespace std;
using namespace phys;
using namespace physmath;

class VVjjHelper{

public:

  VVjjHelper(Histogrammer *histopointer){
  	// theHistograms
    histo_ = histopointer;
    
    // MELA spline
    TFile* f = TFile::Open("../../ZZAnalysis/AnalysisStep/data/cconstants/SmoothKDConstant_m4l_DjjVBF13TeV.root");
  	MELASpline = *((TSpline3*)(f->Get("sp_gr_varReco_Constant_Smooth")->Clone()));
  	f->Close();
  	delete f;    
  }

  virtual ~VVjjHelper(){
    delete histo_;
  }
  
  static void test(int number);

  
 private:

  friend class WZAnalyzer;
  friend class ZZjjAnalyzer;

  
  // Data memebers
  const float rangeVmass = 30.;
  Histogrammer *histo_;
  TSpline3 MELASpline;
  
  
  // Private member functions
  void printTime(float btime, float etime);
  
  float getSpline(double ZZmass) const;
  
  void FindLeadingJets(const vector<Particle> *jetcollection, Particle &Jet0, Particle &Jet1, const vector<Particle> *particlecollection);
  void FindLeadingJets(const vector<Jet> *jetcollection, Particle &Jet0, Particle &Jet1);

  
  // Histogram functions
  void PlotParticle(const Particle &particle, string name, const float weight, string suffix);
  void PlotJets(const Particle &Jet0, const Particle &Jet1, string prename, const float weight, string suffix);

  template <class BOS>
  void PlotBoson(const BOS &particle, string name, const float weight, string suffix);
  
  template <class DiBOS>
  void PlotDiBoson(const DiBOS& particle, string name, const float weight, string suffix);

};

#endif
