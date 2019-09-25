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
#include "VVXAnalysis/DataFormats/interface/Boson.h"
#include "VVXAnalysis/DataFormats/interface/DiBoson.h"
#include "VVXAnalysis/DataFormats/interface/TypeDefs.h"
#include "VVXAnalysis/DataFormats/interface/MELA.h"
#include "VVXAnalysis/DataFormats/interface/GenEventWeights.h"

#include "VVXAnalysis/TreeAnalysis/interface/Histogrammer.h"
#include "VVXAnalysis/TreeAnalysis/interface/MCInfo.h"

#include "VVXAnalysis/Commons/interface/RegionTypes.h"
#include "VVXAnalysis/Commons/interface/LeptonScaleFactors.h"
//#include "VVXAnalysis/Commons/interface/Utilities.h"
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

  static double deltaR (double rapidity1, double phi1, double rapidity2, double phi2 ) {
    
    double Drapidity = fabs(rapidity1 - rapidity2);
    double Dphi      = fabs(phi1 - phi2);
    if (Dphi>TMath::Pi()) Dphi -= TMath::TwoPi();

    double dR = sqrt(Drapidity*Drapidity  + Dphi*Dphi);
    
    return dR;
    
  }
  template <typename T>
    static double deltaR (const T&a, const T&b ) {
    
    double Drapidity = a.Eta() - b.Eta();
    double Dphi      = a.Phi() - b.Phi(); 
    if (Dphi>TMath::Pi()) Dphi -= TMath::TwoPi();
   
    double dR = sqrt(Drapidity*Drapidity  + Dphi*Dphi);
    
    return dR;
  }
  
  std::string fileName;
  
    static double deltaPhi (const TLorentzVector &a, const TLorentzVector &b) {
    
    double phi1 = a.Phi();
    double phi2 = b.Phi();
    
    double DPhi = std::abs(phi2 - phi1);
    
    if(DPhi > TMath::Pi()) DPhi = (TMath::TwoPi()) - DPhi;
    
    return DPhi;   
  }
    
  // Class for specific selections
  SelectorBase& select;

  // To steer the loop over all events. User is not supposed to change this.
  virtual void     loop(const std::string outputfile);

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
  
  void applyLeptonScaleFactors();
  
  // Some basic plots. User may want to change these, thou they should be used only for very basic plots.
  virtual void     fillBasicPlots();
  void fillParticlePlots         (const std::string &type, const phys::Particle &particle);
  void fillLeptonPlots           (const std::string &type, const phys::Lepton   &lepton);
  void fillJetPlots              (const std::string &type, const phys::Jet      &jet);

 private:
  TTree *theTree;
  LeptonScaleFactors leptonScaleFactors_;

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
  int unweightedEventsInSR;
  int unweightedEventsIn2P2FCR;
  int unweightedEventsIn3P1FCR;

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

  // Bosons Candidate
  std::vector<phys::Boson<phys::Jet> >      *VhadCand; TBranch *b_VhadCand;
  
  // Bosons (Not in the tree)
  std::vector<phys::Boson<phys::Jet> >      *Vhad;

  // DiBoson, if in SR, or Z+ll if in CR
  phys::DiBoson<phys::Lepton  , phys::Lepton> *ZZ; TBranch *b_ZZ;

  // Z+L 
  ZLCompositeCandidates *ZLCand; TBranch *b_ZLCand;

  // Z+L  (Not in the tree)
  ZLCompositeCandidates *ZL;


  // GenParticle 
  std::vector<phys::Particle>               *genParticles;   TBranch *b_genParticles;
  std::vector<phys::Boson<phys::Particle> > *genVBParticles; TBranch *b_genVBParticles;
  
  // GenJets
  std::vector<phys::Particle>               *pgenJets;    TBranch *b_pgenJets;
  std::vector<phys::Particle>               *pgenJetsAK8; TBranch *b_pgenJetsAK8;

  // Jets with pT > 30 GeV and |eta| < 4.7 (not in the tree)
  std::vector<phys::Particle> *genJets;
  std::vector<phys::Particle> *genJetsAK8;

  
  // Central jets (not in the tree)
  std::vector<phys::Particle> *centralGenJets;

};

#endif
