#ifndef VVXAnalysis_TreeAnalysis_EventAnalyzer_H
#define VVXAnalysis_TreeAnalysis_EventAnalyzer_H

/** \class EventAnalyzer
 *  Base class for event analyzers. Analyzers have to inherit from this class 
 *  and implement the pure virtual function analyze(), called each event.
 *
 *  $Date: 2013/03/15 13:37:31 $
 *  $Revision: 1.4 $
 *  \author R. Bellan - UNITO <riccardo.bellan@cern.ch>
 */

#include <vector>
#include <string>
#include <utility>
#include <map>

#include "Electron.h"
#include "Lepton.h"
#include "Jet.h"

#include "Histogrammer.h"
class TFile;
class TTree;
class TBranch;
class TH1;

class EventAnalyzer {
public:
  enum METType {Std,NoMu,NoEl};

  EventAnalyzer(std::string filename, double lumi = 1., double externalXSection = -1., bool doBasicPlots = false);
  virtual ~EventAnalyzer();

  struct PtComparator{
    template<typename LEP>
    bool operator()( const LEP & a , 
		     const LEP & b) const{ 
      return a.p4().Pt() > b.p4().Pt(); 
    }
  };
  
  // To steer the loop over all events. User is not supposed to change this.
  virtual void     loop(const std::string outputfile);

 protected:
  // Functions to be overloaded in the concrete instance of the EventAnalyzer class.
  virtual void  begin() {}
  virtual Int_t cut();
  virtual void  analyze() = 0;
  virtual void  end(TFile &) {};  
  
 private:
  // Structural functions to access the tree branches. User is not supposed to change these.
  virtual Int_t    GetEntry(Long64_t entry);
  virtual Long64_t LoadTree(Long64_t entry);
  virtual void     Init(TTree *tree);
  
  // Some basic plots. User may want to change these, thou they should be used only for very basic plots.
  virtual void     fillBasicPlots();
  void fillParticlePlots         (const std::string &type, const phys::Particle &particle);
  void fillLeptonPlots           (const std::string &type, const phys::Lepton   &lepton);
  void fillJetPlots              (const std::string &type, const phys::Jet      &jet);
  void fillElectronPlots         (const std::string &type, const phys::Electron &electron);
  void fillExtraPlotsForElectrons(const std::string &type, const phys::Electron &electron);

 private:
  TTree *theTree;
  int fCurrent; 
  bool doBasicPlots_;

 protected:
  // Histograms helper class
  Histogrammer histograms;

  double theWeight;
  double theSampleWeight;
  int    theCutCounter;

  // Access to the branches
  Int_t    event    ; TBranch *b_event;
  Int_t    run      ; TBranch *b_run;
  Int_t    lumiBlock; TBranch *b_lumiBlock;
  Int_t    nvtx     ; TBranch *b_nvtx;
  Double_t rho      ; TBranch *b_rho;
  Double_t weight   ; TBranch *b_weight;
  Double_t puweight ; TBranch *b_puweight;
  Double_t xsec     ; TBranch *b_xsec;
  Int_t genCategory ; TBranch *b_genCategory;
  Int_t totalEvents ; TBranch *b_totalEvents;

  //MET
  phys::Particle *met   ; TBranch *b_met;

  //Muons
  std::vector<phys::Lepton> *muons; TBranch *b_muons;    

  //Electrons
  std::vector<phys::Electron> *electrons; TBranch *b_electrons;

  // Jets  
  std::vector<phys::Jet> *jets; TBranch *b_jets;

  // GenParticle 
  std::vector<phys::Particle> *genParticles; TBranch *b_genParticles;
  
};
#endif


