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
#include "DataFormats/Common/interface/MergeableCounter.h"

#include <SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h>
#include "SimDataFormats/GeneratorProducts/interface/GenRunInfoProduct.h"

#include "VVXAnalysis/DataFormats/interface/Lepton.h"
#include "VVXAnalysis/DataFormats/interface/Jet.h"
#include "VVXAnalysis/DataFormats/interface/Boson.h"
#include "VVXAnalysis/DataFormats/interface/DiBoson.h"
#include "VVXAnalysis/DataFormats/interface/GenEventWeights.h"
#include "VVXAnalysis/DataFormats/interface/MELA.h"


#include "VVXAnalysis/Producers/interface/FilterController.h"
#include "VVXAnalysis/Commons/interface/LeptonScaleFactors.h"
#include "ZZAnalysis/AnalysisStep/interface/PileUpWeight.h"
#include "DataFormats/PatCandidates/interface/MET.h"
#include <JetMETCorrections/Modules/interface/JetResolution.h>

#include <CommonLHETools/LHEHandler/interface/LHEHandler.h>
#include <JHUGenMELA/MELA/interface/Mela.h>

class TTree;
namespace pat{class Jet;}

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
  
  bool fillGenInfo(const edm::Event& event);
  bool fillEventInfo(const edm::Event& event);
  
  template<typename LEP>
  phys::Lepton fillLepton(const LEP& particle) const;

  phys::Lepton fill(const pat::Muon &muon) const;

  phys::Lepton fill(const pat::Electron &electron) const;

  phys::Lepton fill(const reco::RecoCandidate& lep) const;

  phys::Jet fill(const pat::Jet &jet) const;

  std::vector<phys::Boson<phys::Lepton> > fillLepBosons(const edm::Handle<edm::View<pat::CompositeCandidate> > & edmBosons, int type) const;
  std::vector<phys::Boson<phys::Jet> >    fillHadBosons(const edm::Handle<edm::View<pat::CompositeCandidate> > & edmBosons, int type) const;


  template<typename T, typename PAR>
    phys::Boson<PAR> fillBoson(const pat::CompositeCandidate & v, int type, bool requireQualityCriteria) const;


  template<typename T1, typename T2>
    phys::DiBoson<phys::Lepton,phys::Lepton> fillDiBoson(const pat::CompositeCandidate& edmDiBosons) const;

  std::vector<phys::DiBoson<phys::Lepton,phys::Lepton> > fillDiBosons(const edm::Handle<edm::View<pat::CompositeCandidate> > & edmDiBosons) const;

  std::vector<std::pair<phys::Boson<phys::Lepton>, phys::Lepton> > fillZLCandidates(const edm::Handle<edm::View<pat::CompositeCandidate> > & edmZLs) const;

  phys::DiBoson<phys::Lepton,phys::Lepton> fillZWCandidate(const edm::Handle<edm::View<pat::CompositeCandidate> > & edmZWs) const;

  

  int computeRegionFlag(const pat::CompositeCandidate & vv) const;

 private:
  struct MinPairComparator{
    bool operator()( const std::pair<int,double> & a , 
		     const std::pair<int,double> & b) const{ 
      return a.second < b.second; 
    }
  };

  
 private:

  TTree *theTree;

  Int_t            setup_;
  PileUpWeight     PUWeighter_; 
  FilterController filterController_;
  // To get Lepton efficiency scale factors. Temporary here!
  LeptonScaleFactors leptonScaleFactors_;
  Int_t signalDefinition_;

  // ------------------- Event info in the tree ------------------- //
  Long64_t event_;
  Int_t run_;
  Int_t lumiBlock_;
  
  Bool_t  passTrigger_;
  Bool_t  passSkim_;
  Short_t triggerWord_;
  
  Int_t genCategory_;

  Int_t preSkimCounter_;
  Int_t postSkimCounter_;
  Int_t postSkimSignalCounter_;
  Int_t signalCounter_;
  Int_t postSkimSignalEvents_;

  phys::GenEventWeights genEventWeights_;
  phys::MELA MELA_;


  phys::Particle  met_;

  Int_t           nvtx_;
  
  // ------------------- Objects in the tree ------------------- //
  // all good isolated muons BUT the ones coming from ZZ decay
  std::vector<phys::Lepton>                 muons_;
  // all good isolated electrons BUT the ones coming from ZZ decay
  std::vector<phys::Lepton>                 electrons_;
  // jets which do not contains leptons from ZZ or other good isolated leptons
  std::vector<phys::Jet>                    jets_;

  // fat jets a-kT R = 0.8
  std::vector<phys::Jet>                    jetsAK8_;

  // V --> jj, with V = W,Z
  std::vector<phys::Boson<phys::Jet>      > Vhad_;

  // ZZ in the SR, Z + ll in the CR
  phys::DiBoson<phys::Lepton  , phys::Lepton> ZZ_;
  
  // Z + 1 loose lepton, for fake rate studies
  std::vector<std::pair<phys::Boson<phys::Lepton>, phys::Lepton> > ZL_;

  // ZW in the SR
  phys::DiBoson<phys::Lepton  , phys::Lepton> ZW_;


  std::vector<phys::Particle>               genParticles_;
  std::vector<phys::Particle>               genTaus_;
  std::vector<phys::Boson<phys::Particle> > genVBParticles_;
  std::vector<phys::Particle>               genJets_;
  std::vector<phys::Particle>               genJetsAK8_;

  // ------------------- Input Labels ------------------- //
  edm::EDGetTokenT<pat::MuonCollection>                 theMuonToken;
  edm::EDGetTokenT<pat::ElectronCollection>             theElectronToken;
  edm::EDGetTokenT<std::vector<pat::Jet> >              theJetToken;
  edm::EDGetTokenT<std::vector<pat::Jet> >              theJetAK8Token;
  edm::EDGetTokenT<edm::View<pat::CompositeCandidate> > theVhadToken;
  edm::EDGetTokenT<edm::View<pat::CompositeCandidate> > theZZToken;
  edm::EDGetTokenT<edm::View<pat::CompositeCandidate> > theZLToken;
  edm::EDGetTokenT<edm::View<pat::CompositeCandidate> > theZWToken;
  edm::EDGetTokenT<pat::METCollection>                  theMETToken;
  edm::EDGetTokenT<std::vector<reco::Vertex> >          theVertexToken;
  edm::EDGetTokenT<double>                              theRhoToken;
  edm::EDGetTokenT<float>                               thekfactorToken_ggZZ    ;
  edm::EDGetTokenT<float>                               thekfactorToken_qqZZM   ;
  edm::EDGetTokenT<float>                               thekfactorToken_qqZZPt  ;
  edm::EDGetTokenT<float>                               thekfactorToken_qqZZdPhi;
  edm::EDGetTokenT<float>                               thekfactorToken_EWKqqZZ ;

  //edm::EDGetTokenT<pat::METCollection>                  theMETNoHFToken;

  // thePUInfoLabel;
  edm::EDGetTokenT<int>                         theGenCategoryToken;
  edm::EDGetTokenT<edm::View<reco::Candidate> > theGenVBCollectionToken;
  edm::EDGetTokenT<edm::View<reco::Candidate> >	theGenCollectionToken;
  edm::EDGetTokenT<edm::View<reco::Candidate> >	theGenTauCollectionToken;
  edm::EDGetTokenT<edm::View<reco::Candidate> > theGenJetCollectionToken;
  edm::EDGetTokenT<edm::View<reco::Candidate> > theGenJetAK8CollectionToken;
  edm::EDGetTokenT<GenEventInfoProduct>         theGenInfoToken;
  edm::EDGetTokenT<GenRunInfoProduct>           theGenInfoTokenInRun;

  // Tokens for counters
  edm::EDGetTokenT<edm::MergeableCounter> thePreSkimCounterToken;
  edm::EDGetTokenT<edm::MergeableCounter> prePreselectionCounterToken_;
  edm::EDGetTokenT<edm::MergeableCounter> postPreselectionCounterToken_;
  edm::EDGetTokenT<edm::MergeableCounter> signalCounterToken_;
  edm::EDGetTokenT<edm::MergeableCounter> postSkimSignalCounterToken_;
  edm::EDGetTokenT<edm::MergeableCounter> srCounterToken_;
  edm::EDGetTokenT<edm::MergeableCounter> cr2P2FCounterToken_;
  edm::EDGetTokenT<edm::MergeableCounter> cr3P1FCounterToken_;          

  // --------------------------------------------------------- //

  // Ordinary data members
  std::string sampleName_;

  bool isMC_;
  bool isSignal_;
  int  sampleType_;
  bool applyTrigger_;
  bool applySkim_;
  bool applyMCSel_;
  bool addLHEKinematics_;

  std::vector<double> theXSections;
  double rho_;


  
  LHEHandler* theLHEHandler; 


  double externalCrossSection_; 
  Double_t summcprocweights_; 
  Double_t sumpuweights_;
  Double_t sumpumcprocweights_;
  Int_t theNumberOfEvents;
  Int_t theNumberOfAnalyzedEvents;
  Int_t eventsInEtaAcceptance_;
  Int_t eventsInEtaPtAcceptance_;
  Int_t eventsInSR_;
  Int_t eventsIn2P2FCR_;
  Int_t eventsIn3P1FCR_;

  std::vector<std::string> skimPaths_;
};
#endif

