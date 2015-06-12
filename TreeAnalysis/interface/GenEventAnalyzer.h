#ifndef VVXAnalysis_TreeAnalysis_GenEventAnalyzer_H
#define VVXAnalysis_TreeAnalysis_GenEventAnalyzer_H

/** \class GenEventAnalyzer
 *  Base class for event analyzers for MadGraph test events. Analyzers have to inherit from this class
 *  and implement the pure virtual function analyze(), called each event.
 *
 */

#include <vector>
#include <string>
#include <utility>
#include <algorithm>
#include <map>

#include "VVXAnalysis/DataFormats/interface/Particle.h"
#include "VVXAnalysis/DataFormats/interface/Boson.h"
#include "VVXAnalysis/TreeAnalysis/interface/Histogrammer.h"
//#include "VVXAnalysis/TreeAnalysis/interface/GenMCInfo.h"
#include "VVXAnalysis/Commons/interface/SignalDefinitions.h"

class TFile;
class TTree;
class TBranch;
class TH1;

class GenEventAnalyzer {
public:

  GenEventAnalyzer(std::string filename, double lumi = 1., float xsec = 1.);

  virtual ~GenEventAnalyzer();

  // To steer the loop over all events. User is not supposed to change this.
  virtual void     loop(const std::string outputfile);

 protected:
  // Functions to be overloaded in the concrete instance of the GenEventAnalyzer class.
  virtual void  begin() {};
  virtual Int_t cut();
  virtual void  analyze();
  virtual void  end(TFile &);
  
 private:
  // Structural functions to access the tree branches. User is not supposed to change these.
  virtual Int_t    GetEntry(Long64_t entry);
  virtual Long64_t LoadTree(Long64_t entry);
  virtual void     Init(TTree *tree);

  void makeBasicPlots(const std::string &selection, const zz::SignalTopology& zzSignalTopology);

 private:
  TTree *theTree;
  int fCurrent;
  double lumi_;
  float xsec_;
  float theWeight;
  int numberOfAnalyzedEvents_;
  int numberOfInputEvents_;
  
 protected:
  /*   static const double ZMASS; */
  /*   static const double WMASS; */
  /*   static const double HMASS; */
  
  // Histograms helper class 
  Histogrammer theHistograms;
  
  // MC helper class
  // GenMCInfo theGenMCInfo;
  
/*   double theWeight; */
/*   double theSampleWeight; */
/*   double theCutCounter; */
/*   double theInputWeightedEvents; */
  
  // Access to the branches

  //std::vector<phys::Particle> *genParticlesIn; TBranch * b_genParticlesIn;  
  std::vector<phys::Particle>               *genParticles;   TBranch *b_genParticles;
  std::vector<phys::Boson<phys::Particle> > *genVBParticles; TBranch *b_genVBParticles;
  std::vector<phys::Particle>               *pgenJets;       TBranch *b_pgenJets;

  // Jets with pT > 30 GeV and |eta| < 4.7 (not in the tree)
  std::vector<phys::Particle> *genJets;
  
  // Central jets (not in the tree)
  std::vector<phys::Particle> *centralGenJets;

  int passSignal;
  int passSignalTightFiducialRegion;
  int passHiggsFiducialRegion;

};

#endif
