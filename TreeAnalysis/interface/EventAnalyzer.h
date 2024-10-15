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
//#include "TTree.h"


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
//#include "VVXAnalysis/TreeAnalysis/interface/FeatSelector.h"
//CT: helper for feature tree filling

#include "VVXAnalysis/TreeAnalysis/interface/SampleInfo.h"

#include "VVXAnalysis/DataFormats/interface/RegionTypes.h"
//#include "VVXAnalysis/Commons/interface/LeptonScaleFactors.h"
#include "VVXAnalysis/Commons/interface/Utilities.h"
#include "VVXAnalysis/Commons/interface/GenVBHelper.h"
#include <ZZAnalysis/AnalysisStep/interface/bitops.h>



class TFile;
class TTree;
class TBranch;
class TH1;
class TString;
struct FeatList{
  double f_weight, f_mll, f_ptGamma, f_ptl1, f_ptl2, f_ptJ0, f_ptJ1, f_etaL0, f_etaL1, f_etaJ0, f_etaJ1, f_dPhiL0G, f_dPhiL1G, f_dPhiLL, f_dPhiJ0G, f_dPhiJ1G, f_dPhiJJ, f_dPhiL0J0, f_dPhiL1J0, f_dPhiL0J1, f_dPhiL1J1, f_recoVMass, f_FWMT0, f_FWMT1, f_FWMT2, f_FWMT3, f_FWMT4;//, f_FWMT5, f_FWMT6;
  int f_nbOfCutsPassed;
};


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
  virtual void  fillFeatTree(FeatList&, bool&) {};
  virtual void  end(TFile &) {};  
  virtual void  finish() {};
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

  void InitOut(FeatList&, TTree*);
    
 private:
  TTree *theTree;
  //LeptonScaleFactors leptonScaleFactors_;

  int fCurrent;
  int  maxNumEvents_;
  bool doBasicPlots_;
  bool doSF;

  
 protected:

  // Region
  std::vector<phys::RegionTypes> regions_;
  phys::RegionTypes region_;  // The current region

  // Histograms helper class
  Histogrammer* theHistograms;  // Points to the Histogrammer for the current region
  std::map<phys::RegionTypes, Histogrammer> mapRegionHisto_;  // Maps every region to a corresponding Histogrammer

  //CT: Feature Tree helper class
  TTree* theFeatTree;  // Points to the FeatSelector for the current region
  //  std::map<phys::RegionTypes, TTree> mapRegionTree_;  // Maps every region to a corresponding FeatSelector
  //  std::map<TString, Double_t>* featList_;//AKA typedef featMap featMap_
  //  std::map<phys::RegionTypes, std::map<TString, Double_t>*> mapRegionFeatList_;
  bool doFeats_;

  
  // MC helper class
  SampleInfo theSampleInfo;

  double theWeight;
  double theSampleWeight;
  double theCutCounter;
  double theInputWeightedEvents;

  // Counters about SR and CRs, with trigger requirements too. The numbers are unweighted.
  int year;
  std::string subEra_;
  
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

  std::bitset<32> regionWord;

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
