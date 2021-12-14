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
#include <algorithm>
#include <map>
#include <bitset>

#include "AnalysisConfiguration.h"

#include "VVXAnalysis/DataFormats/interface/Lepton.h"
#include "VVXAnalysis/DataFormats/interface/Jet.h"
#include "VVXAnalysis/DataFormats/interface/Photon.h"
#include "VVXAnalysis/DataFormats/interface/Boson.h"
#include "VVXAnalysis/DataFormats/interface/DiBoson.h"
#include "VVXAnalysis/DataFormats/interface/TypeDefs.h"
#include "VVXAnalysis/DataFormats/interface/MELA.h"
#include "VVXAnalysis/DataFormats/interface/GenEventWeights.h"

#include "VVXAnalysis/TreeAnalysis/interface/Histogrammer.h"
#include "VVXAnalysis/TreeAnalysis/interface/MCInfo.h"

#include "VVXAnalysis/DataFormats/interface/RegionTypes.h"
//#include "VVXAnalysis/Commons/interface/LeptonScaleFactors.h"
#include "VVXAnalysis/Commons/interface/Utilities.h"
#include "VVXAnalysis/Commons/interface/GenVBHelper.h"
#include <ZZAnalysis/AnalysisStep/interface/bitops.h>



class TFile;
class TTree;
class TBranch;
class TH1;

class SelectorBase {
public:
  virtual bool operator()(const phys::Boson<phys::Lepton>&)   const = 0;
  virtual bool operator()(const phys::Boson<phys::Jet>&)      const = 0;
  virtual ~SelectorBase() = 0;
};


template <typename Analysis>
class Selector : public SelectorBase {
public:
  Analysis& analysis;
  Selector(Analysis& ananalysis) : analysis(ananalysis) { };
  virtual ~Selector(){};
  virtual bool operator()(const phys::Boson<phys::Lepton>&   boson) const { return analysis.ZBosonDefinition(boson); }
  virtual bool operator()(const phys::Boson<phys::Jet>&      boson) const { return analysis.WBosonDefinition(boson); }
};


class EventAnalyzer {
public:
  
  enum METType {Std,NoMu,NoEl};
  typedef std::pair<phys::Boson<phys::Lepton>, phys::Lepton> ZLCompositeCandidate;
  typedef std::vector<ZLCompositeCandidate> ZLCompositeCandidates;

  EventAnalyzer(SelectorBase& aSelector, const AnalysisConfiguration& configuration);

  virtual ~EventAnalyzer();

  
  std::string fileName;
  
    
  // Class for specific selections
  SelectorBase& select;

  // To steer the loop over all events. User is not supposed to change this.
  virtual void     loop(const std::string outputfile);

  inline void printEvent(int opt = 0) const {
    if(opt == 0)      std::cout << "Run: " << run << " luminosity block: " << lumiBlock << " event: " << event << std::endl; 
    else if(opt == 1) std::cout << run << ":" << lumiBlock << ":" << event << std::endl; 
    else              std::cout << "Run: " << run << " luminosity block: " << lumiBlock << " event: " << event << std::endl; 
  }

 protected:
  // Functions to be overloaded in the concrete instance of the EventAnalyzer class.
  virtual void  begin() {}
  virtual Int_t cut();
  virtual void  analyze() = 0;
  virtual void  end(TFile &) {};  
  virtual void  addOptions(){};
  TTree * tree() const {return theTree;}
  
 private:
  // Structural functions to access the tree branches. User is not supposed to change these.
  virtual Int_t    GetEntry(Long64_t entry);
  virtual Long64_t LoadTree(Long64_t entry);
  virtual void     Init(TTree *tree);
  
  //void applyLeptonScaleFactors();
  
  // Some basic plots. User may want to change these, thou they should be used only for very basic plots.
  virtual void     fillBasicPlots();
  void fillParticlePlots         (const std::string &type, const phys::Particle &particle);
  void fillLeptonPlots           (const std::string &type, const phys::Lepton   &lepton);
  void fillJetPlots              (const std::string &type, const phys::Jet      &jet);

 private:
  TTree *theTree;
  //LeptonScaleFactors leptonScaleFactors_;

  int fCurrent;
  int  maxNumEvents_;
  bool doBasicPlots_;
  bool doSF;

 protected:

  // Region
   phys::RegionTypes region_;

  // Histograms helper class
  Histogrammer theHistograms;

  // MC helper class
  MCInfo theMCInfo;

  double theWeight;
  double theSampleWeight;
  double theCutCounter;
  double theInputWeightedEvents;

  // Counters about SR and CRs, with trigger requirements too. The numbers are unweighted.
  int year;
  
  // Access to the branches
  Long64_t event     ; TBranch *b_event;
  Int_t    run       ; TBranch *b_run;
  Int_t    lumiBlock ; TBranch *b_lumiBlock;
  Int_t    nvtx      ; TBranch *b_nvtx;
  
  TBranch *b_genEventWeights;

  //TBranch *b_mcprocweight;
  Int_t genCategory  ; TBranch *b_genCategory;

  std::bitset<16> topology;

  Bool_t  passTrigger; TBranch *b_passTrigger;
  Bool_t  passSkim   ; TBranch *b_passSkim;
  Short_t triggerWord; TBranch *b_triggerWord;
  Int_t   pregionWord;  TBranch *b_regionWord;

  phys::MELA *mela;
  TBranch *b_mela;

  std::bitset<128> regionWord;

  //MET
  phys::Particle *met   ; TBranch *b_met;

  //Muons
  std::vector<phys::Lepton> *muons; TBranch *b_muons;    

  //Electrons
  std::vector<phys::Lepton> *electrons; TBranch *b_electrons;

  // Persistent Jets (no eta cut, pT > 20 GeV)  
  std::vector<phys::Jet> *pjets;    TBranch *b_pjets;
  std::vector<phys::Jet> *pjetsAK8; TBranch *b_pjetsAK8;

  // Jets with pT > 30 GeV and |eta| < 4.7 (not in the tree)
  std::vector<phys::Jet> *jets;
  std::vector<phys::Jet> *jetsAK8;
  
  // Central jets (not in the tree)
  std::vector<phys::Jet> *centralJets;
  
  // Photons
  std::vector<phys::Photon> *photons; TBranch *b_photons;

  // Bosons Candidate
  std::vector<phys::Boson<phys::Jet> >      *VhadCand; TBranch *b_VhadCand;
  
  // Bosons (Not in the tree)
  std::vector<phys::Boson<phys::Jet> >      *Vhad;

  // DiBoson, if in SR, or Z+ll if in CR
  phys::DiBoson<phys::Lepton  , phys::Lepton> *ZZ; TBranch *b_ZZ;

  // DiBoson, WZ
  phys::DiBoson<phys::Lepton  , phys::Lepton> *ZW; TBranch *b_ZW;

  // Z->ll
  phys::Boson<phys::Lepton>                    *Z; TBranch *b_Z;

  // Z+L 
  ZLCompositeCandidate *ZL; TBranch *b_ZLCand;

  // GenParticle 
  std::vector<phys::Particle>               *genParticles;   TBranch *b_genParticles;
  std::vector<phys::Particle>               *genTaus;        TBranch *b_genTaus;
  std::vector<phys::Boson<phys::Particle> > *genVBParticles; TBranch *b_genVBParticles;
  
  // GenJets
  std::vector<phys::Particle>               *pgenJets;    TBranch *b_pgenJets;
  std::vector<phys::Particle>               *pgenJetsAK8; TBranch *b_pgenJetsAK8;

  // Jets with pT > 30 GeV and |eta| < 4.7 (not in the tree)
  std::vector<phys::Particle> *genJets;
  std::vector<phys::Particle> *genJetsAK8;

  
  // Central jets (not in the tree)
  std::vector<phys::Particle> *centralGenJets;

  // Helper for gen VB analysis
  GenVBHelper genVBHelper_;

  
};

#endif
