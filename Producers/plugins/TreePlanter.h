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
#include "VVXAnalysis/DataFormats/interface/DiBoson.h"

#include "VVXAnalysis/Producers/interface/FilterController.h"

#include "ZZAnalysis/AnalysisStep/interface/PUReweight.h"

class TTree;
namespace cmg{class PFJet;}
class MCHistoryTools;


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
  
  bool fillEventInfo(const edm::Event& event);
  
  template<typename LEP>
    phys::Lepton fillLepton(const LEP& particle) const;
  
  phys::Lepton fill(const pat::Muon &muon) const;

  phys::Electron fill(const pat::Electron &electron) const;
  
  phys::Jet fill(const cmg::PFJet &jet) const;
  
  template<typename T, typename PAR>
    std::vector<phys::Boson<PAR> > fillBosons(const edm::Handle<edm::View<pat::CompositeCandidate> > & edmBosons, int type = 23) const;

  template<typename T, typename PAR>
    phys::Boson<PAR> fillBoson(const pat::CompositeCandidate & v, int type = 23) const;

  template<typename T1, typename PAR1, typename T2, typename PAR2>
    std::vector<phys::DiBoson<PAR1,PAR2> > fillDiBosons(const edm::Handle<edm::View<pat::CompositeCandidate> > & edmDiBosons) const;

  int computeCRFlag(const pat::CompositeCandidate & vv) const;

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
  FilterController filterController_;
  MCHistoryTools   *mcHistoryTools_;
  // ------------------- Event info in the tree ------------------- //
  Int_t event_;
  Int_t run_;
  Int_t lumiBlock_;
  
  Bool_t  passTrigger_;
  Bool_t  passSkim_;
  Short_t triggerWord_;
  

  Int_t preSkimCounter_;
  Int_t postSkimCounter_;

  Double_t mcprocweight_;
  Double_t puweight_;

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

  std::vector<phys::Boson<phys::Lepton>   > Zmm_;
  std::vector<phys::Boson<phys::Electron> > Zee_;
  std::vector<phys::Boson<phys::Jet>      > Wjj_;

  std::vector<phys::DiBoson<phys::Lepton  , phys::Lepton>   > ZZ4m_;
  std::vector<phys::DiBoson<phys::Electron, phys::Electron> > ZZ4e_;
  std::vector<phys::DiBoson<phys::Electron, phys::Lepton>   > ZZ2e2m_;

  std::vector<phys::Particle>               genParticles_;
  std::vector<phys::Boson<phys::Particle> > genVBParticles_;

  // ------------------- Input Labels ------------------- //
  edm::InputTag theMuonLabel;
  edm::InputTag theElectronLabel;
  edm::InputTag theJetLabel;
  edm::InputTag theZmmLabel;
  edm::InputTag theZeeLabel;
  edm::InputTag theWLabel;
  edm::InputTag theZZ4mLabel;
  edm::InputTag theZZ4eLabel;
  edm::InputTag theZZ2e2mLabel;
  edm::InputTag theMETLabel;
  edm::InputTag theVertexLabel;
  edm::InputTag thePUInfoLabel;
  edm::InputTag theGenCategoryLabel;
  edm::InputTag theGenVBCollectionLabel;
  edm::InputTag	theGenCollectionLabel;

  // --------------------------------------------------------- //

  // Ordinary data members
  bool isMC_;
  int  sampleType_;
  int  setup_;
  bool applyTrigger_;
  bool applySkim_;
  bool applyMCSel_;

  std::vector<double> theXSections;
  double externalCrossSection_;
  Double_t summcprocweights_;
  Double_t sumpuweights_;
  Double_t sumpumcprocweights_;
  Int_t theNumberOfEvents;
  Int_t theNumberOfAnalyzedEvents;


  std::vector<std::string> skimPaths_;
};
#endif

