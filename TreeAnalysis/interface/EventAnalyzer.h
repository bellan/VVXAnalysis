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
  typedef std::pair<const phys::Particle*, phys::Boson<phys::Lepton> > ZParticlePair;
  typedef std::pair<const phys::Particle*, phys::Boson<phys::Jet> >    WParticlePair;

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
  


  struct PtComparator{
    template<typename LEP>
    bool operator()( const LEP & a , 
		     const LEP & b) const{ 
      return a.p4().Pt() > b.p4().Pt(); 
    }
  };

 struct WPtComparator{
    bool operator()( const phys::Boson<phys::Jet> & w1 ,
		     const phys::Boson<phys::Jet> & w2) const{
      return w1.pt() > w2.pt();
    }
  };

 struct WJetPtComparator{
    bool operator()( const phys::Boson<phys::Jet> & w1 ,
		     const phys::Boson<phys::Jet> & w2) const{
      double w1pt1 = w1.daughter(0).pt();
      double w1pt2 = w1.daughter(1).pt();
      if (w1pt2>w1pt1) std::swap(w1pt1,w1pt2);
      double w2pt1 = w2.daughter(0).pt();
      double w2pt2 = w2.daughter(1).pt();
      if (w2pt2>w2pt1) std::swap(w2pt1,w2pt2);

      if (w1pt2==w2pt2) {
	return w1pt1>w2pt1;
      } else return (w1pt2>w2pt2);
    }
  };


  
  struct MassComparator{
    MassComparator(const double& ref): ref_(ref){}
    template<typename PAR>
    bool operator()(const PAR & a , 
		    const PAR & b) const{ 
      return fabs(a.p4().M()-ref_) < fabs(b.p4().M()-ref_); 
    }
    template<typename PAR>
    bool operator()(const PAR * a , 
		    const PAR * b) const{ 
      return fabs(a->p4().M()-ref_) < fabs(b->p4().M()-ref_); 
    }

    double ref_;
  };
 
  // Class for specific selections
  SelectorBase& select;

  struct ZdeltaRComparator{
    bool operator()(const ZParticlePair & a ,
                    const ZParticlePair & b) const{
      return deltaR(a.first->p4().Rapidity(), a.first->p4().Phi(), a.second.p4().Rapidity(), a.second.p4().Phi()) < deltaR(b.first->p4().Rapidity(), b.first->p4().Phi(), b.second.p4().Rapidity(), b.second.p4().Phi());
    }
  };

  struct WdeltaRComparator{
    bool operator()(const WParticlePair & a ,
                    const WParticlePair & b) const{
      return deltaR(a.first->p4().Rapidity(), a.first->p4().Phi(), a.second.p4().Rapidity(), a.second.p4().Phi()) < deltaR(b.first->p4().Rapidity(), b.first->p4().Phi(), b.second.p4().Rapidity(), b.second.p4().Phi());
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
  static const double ZMASS;
  static const double WMASS;
  static const double HMASS;

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


  //MET
  phys::Particle *met   ; TBranch *b_met;

  //Muons
  std::vector<phys::Lepton> *muons; TBranch *b_muons;    

  //Electrons
  std::vector<phys::Electron> *electrons; TBranch *b_electrons;

  // Jets  
  std::vector<phys::Jet> *jets; TBranch *b_jets;

  // Bosons Candidate
  std::vector<phys::Boson<phys::Lepton> >   *ZmmCand; TBranch *b_ZmmCand;
  std::vector<phys::Boson<phys::Electron> > *ZeeCand; TBranch *b_ZeeCand;
  std::vector<phys::Boson<phys::Jet> >      *WjjCand; TBranch *b_WjjCand;
  
  // Bosons (Not in the tree)
  std::vector<phys::Boson<phys::Lepton> >   *Zmm;
  std::vector<phys::Boson<phys::Electron> > *Zee;
  std::vector<phys::Boson<phys::Jet> >      *Wjj;

  // GenParticle 
  std::vector<phys::Particle>               *genParticles;   TBranch *b_genParticles;
  std::vector<phys::Boson<phys::Particle> > *genVBParticles; TBranch *b_genVBParticles;
};

#endif
