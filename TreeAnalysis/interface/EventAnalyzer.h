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

#include "VVXAnalysis/DataFormats/interface/Electron.h"
#include "VVXAnalysis/DataFormats/interface/Lepton.h"
#include "VVXAnalysis/DataFormats/interface/Jet.h"
#include "VVXAnalysis/DataFormats/interface/Boson.h"
#include "VVXAnalysis/DataFormats/interface/DiBoson.h"
#include "VVXAnalysis/DataFormats/interface/TypeDefs.h"

#include "VVXAnalysis/TreeAnalysis/interface/Histogrammer.h"
#include "VVXAnalysis/TreeAnalysis/interface/MCInfo.h"


class TFile;
class TTree;
class TBranch;
class TH1;

class SelectorBase {
public:
  virtual bool operator()(const phys::Boson<phys::Lepton>&)   const = 0;
  virtual bool operator()(const phys::Boson<phys::Electron>&) const = 0;
  virtual bool operator()(const phys::Boson<phys::Jet>&)      const = 0;
};

template <typename Analysis>
class Selector : public SelectorBase {
public:
  Analysis& analysis;
  Selector(Analysis& ananalysis) : analysis(ananalysis) { };
  virtual bool operator()(const phys::Boson<phys::Lepton>&   boson) const { return analysis.ZBosonDefinition(boson); }
  virtual bool operator()(const phys::Boson<phys::Electron>& boson) const { return analysis.ZBosonDefinition(boson); }
  virtual bool operator()(const phys::Boson<phys::Jet>&      boson) const { return analysis.WBosonDefinition(boson); }
};


class EventAnalyzer {
public:
  
  enum METType {Std,NoMu,NoEl};

  EventAnalyzer(SelectorBase& aSelector, std::string filename, double lumi = 1., double externalXSection = -1., bool doBasicPlots = true);

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
  int  triggerChecks_;

 protected:

  // Histograms helper class
  Histogrammer theHistograms;

  // MC helper class
  MCInfo theMCInfo;

  double theWeight;
  double theSampleWeight;
  double theCutCounter;
  double theInputWeightedEvents;

  // Access to the branches
  Int_t    event    ; TBranch *b_event;
  Int_t    run      ; TBranch *b_run;
  Int_t    lumiBlock; TBranch *b_lumiBlock;
  Int_t    nvtx     ; TBranch *b_nvtx;
  Double_t rho      ; TBranch *b_rho;
  
		      TBranch *b_puweight;
  		      TBranch *b_mcprocweight;
  Int_t genCategory ; TBranch *b_genCategory;

  Bool_t  passTrigger; TBranch *b_passTrigger;
  Bool_t  passSkim   ; TBranch *b_passSkim;
  Short_t triggerWord; TBranch *b_triggerWord;
  

  //MET
  phys::Particle *met   ; TBranch *b_met;

  //Muons
  std::vector<phys::Lepton> *muons; TBranch *b_muons;    

  //Electrons
  std::vector<phys::Electron> *electrons; TBranch *b_electrons;

  // Persistent Jets (no eta cut, pT > 20 GeV)  
  std::vector<phys::Jet> *pjets; TBranch *b_pjets;

  // Jets with pT > 30 GeV and |eta| < 4.7 (not in the tree)
  std::vector<phys::Jet> *jets;
  
  // Central jets (not in the tree)
  std::vector<phys::Jet> *centralJets;

  // Bosons Candidate
  std::vector<phys::Boson<phys::Lepton> >   *ZmmCand; TBranch *b_ZmmCand;
  std::vector<phys::Boson<phys::Electron> > *ZeeCand; TBranch *b_ZeeCand;
  std::vector<phys::Boson<phys::Jet> >      *WjjCand; TBranch *b_WjjCand;
  
  // Bosons (Not in the tree)
  std::vector<phys::Boson<phys::Lepton> >   *Zmm;
  std::vector<phys::Boson<phys::Electron> > *Zee;
  std::vector<phys::Boson<phys::Jet> >      *Wjj;

  // DiBosons
  std::vector<phys::DiBoson<phys::Lepton  , phys::Lepton>   >  *ZZ4m  ; TBranch *b_ZZ4m; 
  std::vector<phys::DiBoson<phys::Electron, phys::Electron> >  *ZZ4e  ; TBranch *b_ZZ4e;  
  std::vector<phys::DiBoson<phys::Electron, phys::Lepton>   >  *ZZ2e2m; TBranch *b_ZZ2e2m;

  // DiBoson (Not in the tree)
  phys::DiBoson<phys::Lepton  , phys::Lepton> *ZZ;

  // GenParticle 
  std::vector<phys::Particle>               *genParticles;   TBranch *b_genParticles;
  std::vector<phys::Boson<phys::Particle> > *genVBParticles; TBranch *b_genVBParticles;
};

#endif
