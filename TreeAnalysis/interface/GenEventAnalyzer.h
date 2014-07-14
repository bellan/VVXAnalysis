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
#include "VVXAnalysis/TreeAnalysis/interface/Histogrammer.h"
//#include "VVXAnalysis/TreeAnalysis/interface/GenMCInfo.h"


class TFile;
class TTree;
class TBranch;
class TH1;

class GenEventAnalyzer {
public:

  GenEventAnalyzer(std::string filename, double lumi = 1.);

  virtual ~GenEventAnalyzer();

  // To steer the loop over all events. User is not supposed to change this.
  virtual void     loop(const std::string outputfile);

 protected:
  // Functions to be overloaded in the concrete instance of the GenEventAnalyzer class.
  virtual void  begin() {};
  virtual Int_t cut();
  virtual void  analyze();
  virtual void  end(TFile &) {};
  
 private:
  // Structural functions to access the tree branches. User is not supposed to change these.
  virtual Int_t    GetEntry(Long64_t entry);
  virtual Long64_t LoadTree(Long64_t entry);
  virtual void     Init(TTree *tree);
  
 private:
  TTree *theTree;
  int fCurrent;
  
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
  
  std::vector<phys::Particle> *genParticles; TBranch *b_genParticles;
  std::vector<phys::Particle> *genParticlesIn; TBranch * b_genParticlesIn;
  
};

#endif
