#ifndef ZZWAnalysis_Producers_TreePlanter_H
#define ZZWAnalysis_Producers_TreePlanter_H

/** \class TreePlanter
 *  No description available.
 *
 *  $Date: $
 *  $Revision: $
 *  \author R. Bellan - UNITO <riccardo.bellan@cern.ch>
 */

#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"

#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/CompositeCandidate.h"

#include "VVXAnalysis/DataFormats/interface/Lepton.h"
#include "VVXAnalysis/DataFormats/interface/Electron.h"
#include "VVXAnalysis/DataFormats/interface/Jet.h"
#include "VVXAnalysis/DataFormats/interface/Boson.h"

#include "ZZAnalysis/AnalysisStep/interface/PUReweight.h"

class TTree;
namespace cmg{class PFJet;}

class TreePlanter: public edm::EDAnalyzer {
  
 public:
  
  /// Constructor
  TreePlanter(const edm::ParameterSet &);
  
  /// Destructor
  virtual ~TreePlanter(){};
  
  // Operations
  virtual void beginJob();
  virtual void endLuminosityBlock(edm::LuminosityBlock const& lumi, edm::EventSetup const& setup);
  virtual void analyze(const edm::Event& event, const edm::EventSetup& setup);
  virtual void endRun(const edm::Run& run, const edm::EventSetup& setup);
  virtual void endJob();
  void initTree();
  
  void fillEventInfo(const edm::Event& event);

  template<typename LEP>
    phys::Lepton fillLepton(const LEP& particle) const;

   phys::Electron fillElectron(const pat::Electron &electron) const;

   phys::Jet fillJet(const cmg::PFJet &jet) const;

   template<typename PAR>
     std::vector<phys::Boson<PAR> > fillBosons(const edm::Handle<edm::View<pat::CompositeCandidate> > & edmBosons,  const std::vector<PAR> & physDaughtersCand, int type = 23) const;

 private:
  struct MinPairComparator{
    bool operator()( const std::pair<int,double> & a , 
		     const std::pair<int,double> & b) const{ 
      return a.second < b.second; 
    }
  };

  
 private:

  TTree *theTree;

  PUReweight       PUWeighter_;
  

  // ------------------- Event info in the tree ------------------- //
  Int_t event_;
  Int_t run_;
  Int_t lumiBlock_;
  

  Double_t mcprocweight_;
  Double_t puweight_;
  Double_t summcprocweights_;
  Double_t sumpuweights_;
  Double_t sumpumcprocweights_;

  Double_t xsec_;
  Int_t genCategory_;
  Int_t nobservedPUInt_; 
  Int_t ntruePUInt_;

  phys::Particle  met_;
  Int_t           nvtx_;
  Double_t        rho_;
  
  // ------------------- Objects in the tree ------------------- //
  std::vector<phys::Lepton>                 muons_;
  std::vector<phys::Electron>               electrons_;
  std::vector<phys::Jet>                    jets_;
  std::vector<phys::Boson<phys::Lepton> >   Zmm_;
  std::vector<phys::Boson<phys::Electron> > Zee_;
  std::vector<phys::Boson<phys::Jet> >      Wjj_;
  std::vector<phys::Particle>               genParticles_;

  // ------------------- Input Labels ------------------- //
  edm::InputTag theMuonLabel;
  edm::InputTag theElectronLabel;
  edm::InputTag theJetLabel;
  edm::InputTag theZmmLabel;
  edm::InputTag theZeeLabel;
  edm::InputTag theWLabel;
  edm::InputTag theMETLabel;
  edm::InputTag theVertexLabel;
  edm::InputTag thePUInfoLabel;
  edm::InputTag theGenCategoryLabel;  
  edm::InputTag	theGenCollectionLabel;

  // --------------------------------------------------------- //

  // Ordinary data members
  bool isMC_;
  int  sampleType_;
  int  setup_;
  std::vector<double> theXSections;
  double externalCrossSection_;
  Int_t theNumberOfEvents;
  Int_t theNumberOfAnalyzedEvents;
};
#endif

