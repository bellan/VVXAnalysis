#include "VVXAnalysis/Producers/plugins/TreePlanter.h"
#include "FWCore/Framework/interface/ConsumesCollector.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/Framework/interface/LuminosityBlock.h"
#include "FWCore/Framework/interface/Run.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/Common/interface/View.h"
#include "DataFormats/Common/interface/MergeableCounter.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/Common/interface/Ptr.h"
#include "DataFormats/JetReco/interface/GenJet.h"
#include <JetMETCorrections/Modules/interface/JetResolution.h>
#include <DataFormats/PatCandidates/interface/PFParticle.h>

#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"

#include "SimDataFormats/GeneratorProducts/interface/LHERunInfoProduct.h"
#include <SimDataFormats/GeneratorProducts/interface/LHEEventProduct.h>

#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/MET.h"
#include "VVXAnalysis/Commons/interface/Utilities.h"
#include "VVXAnalysis/Commons/interface/Constants.h"

#include "ZZAnalysis/AnalysisStep/interface/bitops.h"
#include "VVXAnalysis/DataFormats/interface/GenStatusBit.h"

#include "TTree.h"
#include "TVectorD.h"

#include <boost/foreach.hpp>
#define foreach BOOST_FOREACH

#include <boost/assign/std/vector.hpp>

#include <typeinfo>


using namespace boost::assign;

using std::cout;
using std::endl;


TreePlanter::TreePlanter(const edm::ParameterSet &config)
  : setup_(config.getParameter<int>("setup"))
  , PUWeighter_      (setup_,setup_)
  , filterController_(config,consumesCollector())
  , leptonScaleFactors_(setup_,
			// FIXME need to be updated for Run2Legacy
			edm::FileInPath("VVXAnalysis/Commons/data/fakeRate_20feb2017.root").fullPath(),
  			edm::FileInPath("VVXAnalysis/Commons/data/fakeRate_20feb2017.root").fullPath())

  , signalDefinition_(config.getParameter<int>("signalDefinition"   ))
  , passTrigger_(false)
  , passSkim_(false)
  , triggerWord_(0)
  , preSkimCounter_  (0)
  , postSkimCounter_ (0)
  , postSkimSignalCounter_(0)
  , signalCounter_(0)
  , postSkimSignalEvents_(0)
  , theMuonToken     (consumes<pat::MuonCollection>                (config.getParameter<edm::InputTag>("muons"    )))
  , theElectronToken (consumes<pat::ElectronCollection>            (config.getParameter<edm::InputTag>("electrons")))
  , theJetToken      (consumes<std::vector<pat::Jet> >             (config.getParameter<edm::InputTag>("jets"     )))
  , theJetAK8Token   (consumes<std::vector<pat::Jet> >             (config.getParameter<edm::InputTag>("jetsAK8"  )))
  , theVhadToken     (consumes<edm::View<pat::CompositeCandidate> >(config.getParameter<edm::InputTag>("Vhad"     )))
  , theZZToken       (consumes<edm::View<pat::CompositeCandidate> >(config.getParameter<edm::InputTag>("ZZ"       )))
  , theZLToken       (consumes<edm::View<pat::CompositeCandidate> >(config.getParameter<edm::InputTag>("ZL"       )))
  , theZWToken       (consumes<edm::View<pat::CompositeCandidate> >(config.getParameter<edm::InputTag>("ZW"       )))
  , theMETToken      (consumes<pat::METCollection>                 (config.getParameter<edm::InputTag>("MET"      )))
  , theVertexToken   (consumes<std::vector<reco::Vertex> >         (config.getParameter<edm::InputTag>("Vertices" )))
  , theRhoToken      (consumes<double>                             (edm::InputTag("fixedGridRhoFastjetAll","")))
  , thekfactorToken_ggZZ     (consumes<float>                      (edm::InputTag("kFactor","ggZZ"           )))
  , thekfactorToken_qqZZM    (consumes<float>                      (edm::InputTag("kFactor","qqZZM"          )))
  , thekfactorToken_qqZZPt   (consumes<float>                      (edm::InputTag("kFactor","qqZZPt"         )))
  , thekfactorToken_qqZZdPhi (consumes<float>                      (edm::InputTag("kFactor","qqZZdPhi"       )))
  , thekfactorToken_EWKqqZZ  (consumes<float>                      (edm::InputTag("kFactor","EWKqqZZ"        )))
  , thePreSkimCounterToken       (consumes<edm::MergeableCounter,edm::InLumi>(edm::InputTag("preSkimCounter"         )))
  , prePreselectionCounterToken_ (consumes<edm::MergeableCounter,edm::InLumi>(edm::InputTag("prePreselectionCounter" )))
  , postPreselectionCounterToken_(consumes<edm::MergeableCounter,edm::InLumi>(edm::InputTag("postPreselectionCounter")))
  , signalCounterToken_          (consumes<edm::MergeableCounter,edm::InLumi>(edm::InputTag("signalCounter"          )))
  , postSkimSignalCounterToken_  (consumes<edm::MergeableCounter,edm::InLumi>(edm::InputTag("postSkimSignalCounter"  )))
  , srCounterToken_              (consumes<edm::MergeableCounter,edm::InLumi>(edm::InputTag("srCounter"              ))) 
  , cr2P2FCounterToken_          (consumes<edm::MergeableCounter,edm::InLumi>(edm::InputTag("cr2P2FCounter"          )))
  , cr3P1FCounterToken_          (consumes<edm::MergeableCounter,edm::InLumi>(edm::InputTag("cr3P1FCounter"          )))
  , sampleName_      (config.getParameter<std::string>("sampleName"))
  , isMC_            (config.getUntrackedParameter<bool>("isMC",false))
  , sampleType_      (config.getParameter<int>("sampleType"))
  , applyTrigger_    (config.getUntrackedParameter<bool>("TriggerRequired", false)) 
  , applySkim_       (config.getUntrackedParameter<bool>("SkimRequired"   , true)) 
  , applyMCSel_      (config.getUntrackedParameter<bool>("DoMCSelection"  , false)) 
  , addLHEKinematics_(config.getParameter<bool>("AddLHEKinematics"))
  , externalCrossSection_(-1.)
  , summcprocweights_    (0.)
  , sumpuweights_        (0.) 
  , sumpumcprocweights_  (0.)
  , theNumberOfEvents(0)
  , theNumberOfAnalyzedEvents(0)
  , eventsInEtaAcceptance_(0)
  , eventsInEtaPtAcceptance_(0)
  , eventsInSR_(0)
  , eventsIn2P2FCR_(0)
  , eventsIn3P1FCR_(0){
 
  edm::Service<TFileService> fs;
  theTree = fs->make<TTree>("ElderTree","ElderTree");

  if(isMC_){
    consumesMany<std::vector< PileupSummaryInfo > >();
    consumesMany<LHEEventProduct>();
    theGenCategoryToken      = consumes<int>                        (config.getUntrackedParameter<edm::InputTag>("GenCategory"    , edm::InputTag("genCategory")));
    theGenCollectionToken    = consumes<edm::View<reco::Candidate> >(config.getUntrackedParameter<edm::InputTag>("GenCollection"  , edm::InputTag("genParticlesFromHardProcess")));
    theGenTauCollectionToken = consumes<edm::View<reco::Candidate> >(config.getUntrackedParameter<edm::InputTag>("GenTaus"  , edm::InputTag("genTaus")));
    //theGenJetCollectionToken = consumes<edm::View<reco::Candidate> >(config.getUntrackedParameter<edm::InputTag>("GenJets"        , edm::InputTag("selectedGenJets")));
    theGenJetCollectionToken    = consumes<edm::View<reco::Candidate> >(config.getUntrackedParameter<edm::InputTag>("GenJets"        , edm::InputTag("genCategory","genJets")));
    theGenJetAK8CollectionToken = consumes<edm::View<reco::Candidate> >(config.getUntrackedParameter<edm::InputTag>("GenJetsAK8"      , edm::InputTag("genCategory","genJetsAK8")));
    theGenVBCollectionToken  = consumes<edm::View<reco::Candidate> >(config.getUntrackedParameter<edm::InputTag>("GenVBCollection", edm::InputTag("genCategory","vectorBosons")));
    theGenInfoToken          = consumes<GenEventInfoProduct>          (edm::InputTag("generator"));
    theGenInfoTokenInRun     = consumes<GenRunInfoProduct,edm::InRun>(edm::InputTag("generator"));
    externalCrossSection_    = config.getUntrackedParameter<double>("XSection",-1);
    theLHEHandler            = new LHEHandler(((MELAEvent::CandidateVVMode)(config.getParameter<int>("VVMode")+1))
					      , config.getParameter<int>("VVDecayMode")
					      , (addLHEKinematics_ ? LHEHandler::doHiggsKinematics : LHEHandler::noKinematics)
					      , setup_ // means year
					      , LHEHandler::tryNNPDF30
					      , LHEHandler::tryNLO);
    
  }
   
  initTree();
}


void TreePlanter::beginJob(){
  theTree->Branch("event"     , &event_); 
  theTree->Branch("run"       , &run_); 
  theTree->Branch("lumiBlock" , &lumiBlock_); 

  theTree->Branch("passTrigger" , &passTrigger_); 
  theTree->Branch("passSkim"    , &passSkim_); 
  theTree->Branch("triggerWord" , &triggerWord_); 

  theTree->Branch("genCategory" , &genCategory_);

  theTree->Branch("genEventWeights", &genEventWeights_);
  theTree->Branch("MELA"        , &MELA_);


  theTree->Branch("met"   , &met_);
  theTree->Branch("nvtxs" , &nvtx_);

  theTree->Branch("muons"     , &muons_);
  theTree->Branch("electrons" , &electrons_);
  theTree->Branch("jets"      , &jets_); 
  theTree->Branch("jetsAK8"   , &jetsAK8_); 
  theTree->Branch("VhadCand"  , &Vhad_);
  theTree->Branch("ZZCand"    , &ZZ_); 
  theTree->Branch("ZWCand"    , &ZW_); 
  theTree->Branch("ZLCand"    , &ZL_); 

  theTree->Branch("genParticles"  , &genParticles_);
  theTree->Branch("genTaus"       , &genTaus_);
  theTree->Branch("genVBParticles", &genVBParticles_);
  theTree->Branch("genJets"       , &genJets_);
  theTree->Branch("genJetsAK8"    , &genJetsAK8_);
}

void TreePlanter::endLuminosityBlock(edm::LuminosityBlock const& lumi, edm::EventSetup const& setup)
{
  // Beware: preSkimCounter for H->ZZ->4l means a skim done at path level

  Float_t Nevt_preskim = -1.;
  edm::Handle<edm::MergeableCounter> preSkimCounter;
  if (lumi.getByToken(thePreSkimCounterToken, preSkimCounter)) { // Counter before skim. Does not exist for non-skimmed samples.
    Nevt_preskim = preSkimCounter->value;
  }  
  
  edm::Handle<edm::MergeableCounter> prePathCounter;
  lumi.getByToken(prePreselectionCounterToken_, prePathCounter);       // Counter of input events from the input pattuple

  // Nevt_gen: this is the number before any skim
  if (Nevt_preskim>=0.) {
    theNumberOfEvents += Nevt_preskim; 
  } else {
    theNumberOfEvents += prePathCounter->value;    
  }

  // Beware: pre/post Skim here means before/after the preselection path at Tree building level!
  edm::Handle<edm::MergeableCounter> counter;
  bool found = lumi.getByToken(prePreselectionCounterToken_, counter);
  if(found) preSkimCounter_ += counter->value;
  
  found = lumi.getByToken(postPreselectionCounterToken_, counter);
  if(found) postSkimCounter_ += counter->value;

  found = lumi.getByToken(signalCounterToken_, counter);
  if(found) signalCounter_ += counter->value;

  found = lumi.getByToken(postSkimSignalCounterToken_, counter);
  if(found) postSkimSignalCounter_ += counter->value;

  found = lumi.getByToken(srCounterToken_, counter);
  if(found) eventsInSR_ += counter->value;

  found = lumi.getByToken(cr2P2FCounterToken_, counter);
  if(found) eventsIn2P2FCR_ += counter->value;

  found = lumi.getByToken(cr3P1FCounterToken_, counter);
  if(found) eventsIn3P1FCR_ += counter->value;
}


void TreePlanter::endRun(const edm::Run& run, const edm::EventSetup& setup){
  if(isMC_){
    edm::Handle<GenRunInfoProduct> genRunInfo;
    run.getByToken(theGenInfoTokenInRun, genRunInfo);
    theXSections.push_back(genRunInfo->crossSection());
  }
}

void TreePlanter::endJob(){

  edm::Service<TFileService> fs;
  TTree *countTree = fs->make<TTree>("HollyTree","HollyTree");
  
  countTree->Branch("setup"         , &setup_); 
  countTree->Branch("analyzedEvents", &theNumberOfAnalyzedEvents);
  countTree->Branch("eventsInSR"    , &eventsInSR_);
  countTree->Branch("eventsIn2P2FCR", &eventsIn2P2FCR_);
  countTree->Branch("eventsIn3P1FCR", &eventsIn3P1FCR_);

  if(isMC_){

    Double_t internalCrossSection = 0.;
    foreach(const double &xsec, theXSections) internalCrossSection += xsec;
    
    internalCrossSection = internalCrossSection/theXSections.size();
    
    cout << "================================================"                << endl
	 << "Total number of generated events: " << theNumberOfEvents         << endl
	 << "Total number of analyzed events: "  << theNumberOfAnalyzedEvents << endl
	 <<" Internal cross section: "           << internalCrossSection      << endl
	 <<" External cross section: "           << externalCrossSection_     << endl
	 << "================================================"                << endl;
    
    countTree->Branch("signalDefinition"      , &signalDefinition_);
    countTree->Branch("genEvents"             , &theNumberOfEvents);
    countTree->Branch("internalCrossSection"  , &internalCrossSection);
    countTree->Branch("externalCrossSection"  , &externalCrossSection_);
    countTree->Branch("summcprocweight"       , &summcprocweights_);
    countTree->Branch("sumpuweight"           , &sumpuweights_);
    countTree->Branch("sumpumcprocweight"     , &sumpumcprocweights_);
    countTree->Branch("preSkimCounter"        , &preSkimCounter_);
    countTree->Branch("postSkimCounter"       , &postSkimCounter_);
    countTree->Branch("postSkimSignalCounter" , &postSkimSignalCounter_);
    countTree->Branch("signalCounter"         , &signalCounter_);
    countTree->Branch("postSkimSignalEvents"  , &postSkimSignalEvents_);
    countTree->Branch("eventsInEtaAcceptance"   , &eventsInEtaAcceptance_);
    countTree->Branch("eventsInEtaPtAcceptance" , &eventsInEtaPtAcceptance_);
  }
  
  countTree->Fill();
}


void TreePlanter::initTree(){
  event_     = -1;
  run_       = -1;
  lumiBlock_ = -1;

  passTrigger_ = false;
  passSkim_    = false;
  triggerWord_ = 0;

  genCategory_ = -1;

  genEventWeights_ = phys::GenEventWeights();
  MELA_ = phys::MELA();

  met_    = phys::Particle();
  nvtx_   = -1;
  
  muons_     = std::vector<phys::Lepton>();
  electrons_ = std::vector<phys::Lepton>();
  jets_      = std::vector<phys::Jet>();
  jetsAK8_   = std::vector<phys::Jet>();
  ZZ_        = phys::DiBoson<phys::Lepton  , phys::Lepton>();
  Vhad_      = std::vector<phys::Boson<phys::Jet>      >();   

  ZL_        = std::vector<std::pair<phys::Boson<phys::Lepton>, phys::Lepton> >();
  ZW_        = phys::DiBoson<phys::Lepton  , phys::Lepton>();

  genParticles_ = std::vector<phys::Particle>();
  genTaus_ = std::vector<phys::Particle>();
  genVBParticles_ = std::vector<phys::Boson<phys::Particle> >();
  genJets_ = std::vector<phys::Particle>();
  genJetsAK8_ = std::vector<phys::Particle>();
}

bool TreePlanter::fillEventInfo(const edm::Event& event){
    
  // Check trigger request. Actually, it is a very very loose request, not the actual one, that instead should be
  // asked to the specific final state

  passTrigger_ = filterController_.passTrigger(NONE, event, triggerWord_);
  if (applyTrigger_ && !passTrigger_) return false;


  // Check Skim requests
  passSkim_ = filterController_.passSkim(event, triggerWord_);

  run_       = event.id().run();
  event_     = event.id().event(); 
  lumiBlock_ = event.luminosityBlock();

  edm::Handle<pat::METCollection> met;      event.getByToken(theMETToken    , met);
  met_ = phys::Particle(phys::Particle::convert(met->front().p4()));


  edm::Handle<std::vector<reco::Vertex> > vertices; event.getByToken(theVertexToken, vertices);
  nvtx_ = vertices->size();
    
  edm::Handle<double> rhoHandle;   
  event.getByToken(theRhoToken, rhoHandle);
  
  rho_ = *rhoHandle;

 
  return true;
}



bool TreePlanter::fillGenInfo(const edm::Event& event){
  // Fill some info abut acceptance before cutting away events. Beware: if the signal is defined a-posteriori, we will have a problem. For that case, we need to
  // explicitly check here that we are counting signal and not irreducible background.

  // Apply MC filter
  if (!filterController_.passMCFilter(event)) return false;
    
  
  edm::Handle<edm::View<reco::Candidate> > genParticles; event.getByToken(theGenCollectionToken, genParticles);
  edm::Handle<edm::View<reco::Candidate> > genTaus; event.getByToken(theGenTauCollectionToken, genTaus);
  edm::Handle<GenEventInfoProduct>         genInfo     ; event.getByToken(theGenInfoToken, genInfo);


  // Fill Pile-up info 
  std::vector<edm::Handle<std::vector< PileupSummaryInfo > > >  PupInfos; 
  event.getManyByType(PupInfos);
  edm::Handle<std::vector< PileupSummaryInfo > > puInfo = PupInfos.front(); 
  
  foreach(const PileupSummaryInfo& pui, *puInfo)
    if(pui.getBunchCrossing() == 0) { 
      int ntruePUInt            = pui.getTrueNumInteractions();
      genEventWeights_.puweight_   = PUWeighter_.weight(ntruePUInt);
      genEventWeights_.puweightUp_ = PUWeighter_.weight(ntruePUInt, PileUpWeight::PUvar::VARUP);
      genEventWeights_.puweightDn_ = PUWeighter_.weight(ntruePUInt, PileUpWeight::PUvar::VARDOWN);
      break;
    }
  

  
  // Info about the MC weight
  genEventWeights_.mcprocweight_ = genInfo->weight();
  
  // Sum of weight, particularly imprtant for MCs that return also negative weights
  // or, in general, weighted events
  sumpuweights_       += genEventWeights_.puweight_;
  summcprocweights_   += genEventWeights_.mcprocweight_;
  sumpumcprocweights_ += genEventWeights_.puweight_*genEventWeights_.mcprocweight_;
  
    
  // NLO k-factors
  edm::Handle<float> kFactorHandle_ggZZ     ;   
  edm::Handle<float> kFactorHandle_qqZZM    ;   
  edm::Handle<float> kFactorHandle_qqZZPt   ;   
  edm::Handle<float> kFactorHandle_qqZZdPhi ;   
  edm::Handle<float> kFactorHandle_EWKqqZZ  ;   
  
  if(event.getByToken( thekfactorToken_ggZZ    ,  kFactorHandle_ggZZ    ))  genEventWeights_.kFactor_ggZZ_     = *kFactorHandle_ggZZ    ;
  if(event.getByToken( thekfactorToken_qqZZM   ,  kFactorHandle_qqZZM   ))  genEventWeights_.kFactor_qqZZM_    = *kFactorHandle_qqZZM   ;
  if(event.getByToken( thekfactorToken_qqZZPt  ,  kFactorHandle_qqZZPt  ))  genEventWeights_.kFactor_qqZZPt_   = *kFactorHandle_qqZZPt  ;
  if(event.getByToken( thekfactorToken_qqZZdPhi,  kFactorHandle_qqZZdPhi))  genEventWeights_.kFactor_qqZZdPhi_ = *kFactorHandle_qqZZdPhi;
  if(event.getByToken( thekfactorToken_EWKqqZZ ,  kFactorHandle_EWKqqZZ ))  genEventWeights_.kFactor_EWKqqZZ_  = *kFactorHandle_EWKqqZZ ;

  
  edm::Handle<int> genCategory;
  event.getByToken(theGenCategoryToken, genCategory);
  genCategory_ = *genCategory;
  isSignal_ = ((genCategory_ ^ signalDefinition_) & signalDefinition_) == 0; 
  if(isSignal_) ++postSkimSignalEvents_;
  

  // load only gen leptons and gen photon status 1, from hard process
  event.getByToken(theGenCollectionToken,  genParticles);
  for (edm::View<reco::Candidate>::const_iterator p = genParticles->begin(); p != genParticles->end(); ++p){
    const reco::GenParticle* gp = dynamic_cast<const reco::GenParticle*>(&(*p));
    genParticles_.push_back(phys::Particle(gp->p4(), phys::Particle::computeCharge(gp->pdgId()), gp->pdgId()));
  }

  event.getByToken(theGenTauCollectionToken,  genTaus);
  for (edm::View<reco::Candidate>::const_iterator p = genTaus->begin(); p != genTaus->end(); ++p){
    const reco::GenParticle* gp = dynamic_cast<const reco::GenParticle*>(&(*p));
    genTaus_.push_back(phys::Particle(gp->p4(), phys::Particle::computeCharge(gp->pdgId()), gp->pdgId()));
  }



  
  event.getByToken(theGenVBCollectionToken,  genParticles);
  for(edm::View<reco::Candidate>::const_iterator p = genParticles->begin(); p != genParticles->end(); ++p)
    { 
      if(fabs(p->pdgId()) == 24 || p->pdgId() == 23)
	genVBParticles_.push_back(phys::Boson<phys::Particle>(phys::Particle(p->daughter(0)->p4(), phys::Particle::computeCharge(p->daughter(0)->pdgId()), p->daughter(0)->pdgId()),
							      phys::Particle(p->daughter(1)->p4(), phys::Particle::computeCharge(p->daughter(1)->pdgId()), p->daughter(1)->pdgId()),
							      p->pdgId()));
    }

  // Get the gen jet collection
  edm::Handle<edm::View<reco::Candidate> > genJets;
  event.getByToken(theGenJetCollectionToken,  genJets);
    
  // Still need to clean the genjets. In SignalDefinition, each category cleans the jet collection when it set the bits.
  // However, for the final... REMOVE this comment
  for(edm::View<reco::Candidate>::const_iterator jet = genJets->begin(); jet != genJets->end(); ++jet)
    genJets_.push_back(phys::Particle(jet->p4(), phys::Particle::computeCharge(jet->pdgId()), jet->pdgId()));
 

  // Get the gen jet collection
  edm::Handle<edm::View<reco::Candidate> > genJetsAK8;
  event.getByToken(theGenJetAK8CollectionToken,  genJetsAK8);
  
  // Still need to clean the genjets. In SignalDefinition, each category cleans the jet collection when it set the bits.
  // However, for the final... REMOVE this comment
  for(edm::View<reco::Candidate>::const_iterator jet = genJetsAK8->begin(); jet != genJetsAK8->end(); ++jet)
    genJetsAK8_.push_back(phys::Particle(jet->p4(), phys::Particle::computeCharge(jet->pdgId()), jet->pdgId()));
 





  // LHE information
  edm::Handle<LHEEventProduct> lhe_evt;
  std::vector<edm::Handle<LHEEventProduct> > lhe_handles;
  event.getManyByType(lhe_handles);

  if (lhe_handles.size()>0){
    lhe_evt = lhe_handles.front();
    theLHEHandler->setHandle(&lhe_evt);
    theLHEHandler->extract();
    //    FillLHECandidate();
    
    genEventWeights_.LHEPDFScale_                      = theLHEHandler->getPDFScale();
    genEventWeights_.LHEweight_QCDscale_muR1_muF1_     = theLHEHandler->getLHEWeight(0, 1.);
    genEventWeights_.LHEweight_QCDscale_muR1_muF2_     = theLHEHandler->getLHEWeight(1, 1.);
    genEventWeights_.LHEweight_QCDscale_muR1_muF0p5_   = theLHEHandler->getLHEWeight(2, 1.);
    genEventWeights_.LHEweight_QCDscale_muR2_muF1_     = theLHEHandler->getLHEWeight(3, 1.);
    genEventWeights_.LHEweight_QCDscale_muR2_muF2_     = theLHEHandler->getLHEWeight(4, 1.);
    genEventWeights_.LHEweight_QCDscale_muR2_muF0p5_   = theLHEHandler->getLHEWeight(5, 1.);
    genEventWeights_.LHEweight_QCDscale_muR0p5_muF1_   = theLHEHandler->getLHEWeight(6, 1.);
    genEventWeights_.LHEweight_QCDscale_muR0p5_muF2_   = theLHEHandler->getLHEWeight(7, 1.);
    genEventWeights_.LHEweight_QCDscale_muR0p5_muF0p5_ = theLHEHandler->getLHEWeight(8, 1.);
    genEventWeights_.LHEweight_PDFVariation_Up_        = theLHEHandler->getLHEWeight_PDFVariationUpDn( 1, 1.);
    genEventWeights_.LHEweight_PDFVariation_Dn_        = theLHEHandler->getLHEWeight_PDFVariationUpDn(-1, 1.);
    genEventWeights_.LHEweight_AsMZ_Up_                = theLHEHandler->getLHEWeigh_AsMZUpDn( 1, 1.);
    genEventWeights_.LHEweight_AsMZ_Dn_                = theLHEHandler->getLHEWeigh_AsMZUpDn(-1, 1.);
    
    theLHEHandler->clear();
  }
  
  
  return true;
}

void TreePlanter::analyze(const edm::Event& event, const edm::EventSetup& setup){

  initTree();

  bool goodEvent = isMC_ ? (fillGenInfo(event) && fillEventInfo(event)) : fillEventInfo(event);
  if(!goodEvent) return;
  ++theNumberOfAnalyzedEvents;

  // Load a bunch of objects from the event
  edm::Handle<pat::MuonCollection>       muons          ; event.getByToken(theMuonToken    ,     muons);
  edm::Handle<pat::ElectronCollection>   electrons      ; event.getByToken(theElectronToken, electrons);
  edm::Handle<std::vector<pat::Jet> >    jets           ; event.getByToken(theJetToken     ,      jets);
  edm::Handle<std::vector<pat::Jet> >    jetsAK8        ; event.getByToken(theJetAK8Token  ,   jetsAK8);
  //  edm::Handle<edm::View<pat::CompositeCandidate> > Vhad ; event.getByToken(theVhadToken    ,      Vhad);
  edm::Handle<edm::View<pat::CompositeCandidate> > ZZ   ; event.getByToken(theZZToken      ,        ZZ);
  edm::Handle<edm::View<pat::CompositeCandidate> > ZL   ; event.getByToken(theZLToken      ,        ZL);
  edm::Handle<edm::View<pat::CompositeCandidate> > ZW   ; event.getByToken(theZWToken      ,        ZW);

  // No further selection on muons, electrons, or jets. All is made in the .py file
  foreach(const pat::Muon&     muon    , *muons    ) muons_.push_back(fill(muon));
  foreach(const pat::Electron& electron, *electrons) electrons_.push_back(fill(electron));
  foreach(const pat::Jet&      jet     , *jets     ) jets_.push_back(fill(jet));
  foreach(const pat::Jet&      jet     , *jetsAK8  ) jetsAK8_.push_back(fill(jet)); // FIXME: need jet class extension
  

  // The bosons are selected requiring that their daughters pass the quality criteria to be good daughters
  // Vhad_ = fillHadBosons(Vhad, 24);


  // The bosons have NOT any requirement on the quality of their daughters, only the flag is set (because of the same code is usd for CR too)
  std::vector<phys::DiBoson<phys::Lepton,phys::Lepton> > ZZs = fillDiBosons(ZZ);

  // Fill Z+l pairs for fake rate measurements
  ZL_ = fillZLCandidates(ZL);
  ZW_ = fillZWCandidate(ZW);

  bool oneZZInSR = false;

  if(ZZ->size() > 1) {

    // Debug in case of more than 1 ZZ candidate
    cout << "----------------------------------------------------" << endl;
    cout << "More than one ZZ candidate!! " << ZZ->size() <<  " candidates in event: " << event_ << endl;
    typedef phys::DiBoson<phys::Lepton,phys::Lepton> ZZlep;
    foreach(const ZZlep& zz , ZZs){
      cout << "....................." << endl;
      cout << zz << " SR? " << test_bit(zz.regionWord_,Channel::ZZ) << " CR2P2F? " << test_bit(zz.regionWord_,CRZLLos_2P2F) << " CR3P1F? " << test_bit(zz.regionWord_,CRZLLos_3P1F) << " ZZOnShell " <<test_bit(zz.regionWord_,Channel::ZZOnShell)<<endl;
      cout << "daughter 0: "   << zz.first() << endl;
      cout << "daughter 0.1: " << zz.first().daughter(0) << " is good? " <<  zz.first().daughter(0).isGood() << " pass full sel? " << zz.first().daughter(0).passFullSel() <<  endl;
      cout << "daughter 0.1: " << zz.first().daughter(1) << " is good? " <<  zz.first().daughter(1).isGood() << " pass full sel? " << zz.first().daughter(1).passFullSel() <<  endl;
      cout << "daughter 1: "   << zz.second() << endl;
      cout << "daughter 1.1: " << zz.second().daughter(0) << " is good? " <<  zz.second().daughter(0).isGood() << " pass full sel? " << zz.second().daughter(0).passFullSel() <<  endl;
      cout << "daughter 1.1: " << zz.second().daughter(1) << " is good? " <<  zz.second().daughter(1).isGood() << " pass full sel? " << zz.second().daughter(1).passFullSel() <<  endl;
      cout << "....................." << endl;
      
      if(!zz.passTrigger()) continue;
      
      // If more than a candidate is found, then give precedence to SR type ZZ
      if(test_bit(zz.regionWord_,Channel::ZZ) || test_bit(zz.regionWord_,Channel::ZZOnShell)){
	ZZ_ = zz;   
	oneZZInSR = true;
      }
      // Otherwise, select the ZZ accordingly to the same logic as the ZZ is chosen
      if(!oneZZInSR){
	if(abs(zz.first().mass() - phys::ZMASS) < abs(ZZ_.first().mass() - phys::ZMASS)) 
	  ZZ_ = zz;
	if((zz.first().mass() - ZZ_.first().mass()) < 0.001 && (zz.second().daughter(0).pt()+zz.second().daughter(1).pt()) > (ZZ_.second().daughter(0).pt() + ZZ_.second().daughter(1).pt()) )
	  ZZ_ = zz;
      }
    }
    
    cout << "----------------------------------------------------" << endl;
  }

  else if(ZZs.size() == 1 && ZZs.front().passTrigger()) ZZ_ = ZZs.front();
  // case ZZ == 1 !trigger is missing
  else if(isMC_  && ZL_.empty() && !ZW_.isValid() && !isSignal_ && applySkim_) return;
  else if(!isMC_ && ZL_.empty() && !ZW_.isValid() &&               applySkim_) return;
  
  theTree->Fill();

}

template<typename LEP>
phys::Lepton TreePlanter::fillLepton(const LEP& lepton) const{

  phys::Lepton output(phys::Particle::convert(lepton.p4()),lepton.charge(),lepton.pdgId());
 

  output.dxy_                 = lepton.userFloat("dxy"              );               
  output.dz_                  = lepton.userFloat("dz"               );                
  output.sip_                 = lepton.userFloat("SIP"              );
  output.pfCombRelIso_        = lepton.userFloat("combRelIsoPF"     );  
  output.pfCombRelIsoFSRCorr_ = lepton.userFloat("combRelIsoPFFSRCorr");
  output.matchHLT_            = lepton.userFloat("HLTMatch"         );
  output.isGood_              = lepton.userFloat("isGood"           ) && lepton.userFloat("passCombRelIsoPFFSRCorr");
  if(abs(lepton.pdgId())      == 11){ 
    output.isInCracks_        = lepton.userFloat("isCrack"          );
    output.scEta_             = lepton.userFloat("SCeta"            );
    // Deal with very rare cases when SCeta is out of 2.5 bonds
    if ( output.eta() <= 2.5 && output.scEta_ >= 2.5) output.scEta_ = 2.49;
    else if ( output.eta() >= -2.5 && output.scEta_ <= -2.5) output.scEta_ = -2.49;
    }
  std::pair<double,double> effSF = leptonScaleFactors_.efficiencyScaleFactor(output); 
  output.efficiencySF_    = effSF.first;
  output.efficiencySFUnc_ = effSF.second;
  output.setFakeRateSF(leptonScaleFactors_.fakeRateScaleFactor(output));

  return output;
}

phys::Lepton TreePlanter::fill(const pat::Electron &electron) const{

  return fillLepton(electron);

}


phys::Lepton TreePlanter::fill(const pat::Muon& mu) const{
  return fillLepton(mu);
}


phys::Jet TreePlanter::fill(const pat::Jet &jet) const{
  

  phys::Jet output(phys::Particle::convert(jet.p4()),jet.charge(),1);
  
 
  output.csvtagger_      = jet.hasUserFloat("bTagger")                   ? jet.userFloat("bTagger")                : -999;
  
  output.deepAK8_TvsQCD_      =  jet.bDiscriminator("pfDeepBoostedDiscriminatorsJetTags:TvsQCD");
  output.deepAK8_WvsQCD_      = jet.bDiscriminator("pfDeepBoostedDiscriminatorsJetTags:WvsQCD");
  output.deepAK8MD_TvsQCD_    = jet.bDiscriminator("pfMassDecorrelatedDeepBoostedDiscriminatorsJetTags:TvsQCD");
  output.deepAK8MD_WvsQCD_    = jet.bDiscriminator("pfMassDecorrelatedDeepBoostedDiscriminatorsJetTags:WvsQCD");
  output.deepAK8MD_ZHbbvsQCD_ = jet.bDiscriminator("pfMassDecorrelatedDeepBoostedDiscriminatorsJetTags:ZHbbvsQCD");
  output.deepAK8MD_ZHccvsQCD_ = jet.bDiscriminator("pfMassDecorrelatedDeepBoostedDiscriminatorsJetTags:ZHccvsQCD");
  
  //output.bTaggers = jet.getPairDiscri(); //TEST get the whole list of tags
  
  
  output.qgLikelihood_   = jet.hasUserFloat("qgLikelihood")              ? jet.userFloat("qgLikelihood")           : -999;
  output.fullPuId_       = jet.hasUserInt  ("pileupJetIdUpdated:fullId") ? jet.userInt("pileupJetIdUpdated:fullId"): -999;
  output.passLooseId_    = jet.hasUserFloat("looseJetID")                ? jet.userFloat("looseJetID")             : -999;

  //girth - See http://arxiv.org/abs/1106.3076v2 - eqn (2)
  double gsumcharged = 0.0;
  double gsum = 0.0;
  float sumpt2 = 0;
  float sumpt = 0;

  for(unsigned i=0; i<jet.numberOfDaughters();++i){
    sumpt2 += jet.daughter(i)->pt() * jet.daughter(i)->pt();
    sumpt += jet.daughter(i)->pt();


    const double dr = reco::deltaR(*jet.daughter(i),jet);
    gsum += (dr*( jet.daughter(i)->pt()/jet.pt() ) );
    if(jet.daughter(i)->charge() == 0) continue;
    else
      gsumcharged += (dr*( jet.daughter(i)->pt()/jet.pt() ) );
  }
  output.girth_        = gsum;
  output.girth_charged_ = gsumcharged;
  //end girth
  output.ptd_ = sqrt( sumpt2 ) / sumpt; 
  
  output.jetArea_    = jet.jetArea();
  output.secvtxMass_ = jet.hasUserFloat("vtxMass") ? jet.userFloat("vtxMass") : -999;
  
  output.mcPartonFlavour_ = jet.partonFlavour();
  
  // JEC
  output.rawFactor_ = jet.jecFactor(0);
  // JES
  output.jecUnc_    = jet.hasUserFloat("jec_unc") ? jet.userFloat("jec_unc") : -999;

  // JER
  if(isMC_){                                                                                                         
    output.pt_nojer_    = jet.hasUserFloat("pt_nojer") ? jet.userFloat("pt_nojer") : -999;
    output.pt_jerup_    = jet.hasUserFloat("pt_jerup") ? jet.userFloat("pt_jerup") : -999;
    output.pt_jerdn_    = jet.hasUserFloat("pt_jerdn") ? jet.userFloat("pt_jerdn") : -999;
  }  

  // Variables for AK8 jets

  if(setup_ == 2016 || setup_ == 2017){
    
     output.tau1_           = jet.hasUserFloat("NjettinessAK8:tau1") ? jet.userFloat("NjettinessAK8:tau1") : -999;
     output.tau2_           = jet.hasUserFloat("NjettinessAK8:tau2") ? jet.userFloat("NjettinessAK8:tau2") : -999;
     output.tau3_           = jet.hasUserFloat("NjettinessAK8:tau3") ? jet.userFloat("NjettinessAK8:tau3") : -999;
     output.corrPrunedMass_ = jet.hasUserFloat("ak8PFJetsCHSCorrPrunedMass") ? jet.userFloat("ak8PFJetsCHSCorrPrunedMass") : -999;
  
     output.prunedMass_     = jet.hasUserFloat("ak8PFJetsCHSPrunedMass")     ? jet.userFloat("ak8PFJetsCHSPrunedMass") : -999;
     output.softDropMass_   = jet.hasUserFloat("ak8PFJetsCHSSoftDropMass")   ? jet.userFloat("ak8PFJetsCHSSoftDropMass") : -999;
  
     output.puppiTau1_      = jet.hasUserFloat("ak8PFJetsPuppiValueMap:NjettinessAK8PuppiTau1") ? jet.userFloat("ak8PFJetsPuppiValueMap:NjettinessAK8PuppiTau1") : -999;
     output.puppiTau2_      = jet.hasUserFloat("ak8PFJetsPuppiValueMap:NjettinessAK8PuppiTau1") ? jet.userFloat("ak8PFJetsPuppiValueMap:NjettinessAK8PuppiTau1") : -999;
     output.puppiTau3_      = jet.hasUserFloat("ak8PFJetsPuppiValueMap:NjettinessAK8PuppiTau1") ? jet.userFloat("ak8PFJetsPuppiValueMap:NjettinessAK8PuppiTau1") : -999;
     output.puppiMass_      = jet.hasUserFloat("ak8PFJetsPuppiValueMap:mass") ? jet.userFloat("ak8PFJetsPuppiValueMap:mass") : -999;

     // Try with new schema
     output.tau1_           = jet.hasUserFloat("ak8PFJetsCHSValueMap:NjettinessAK8CHSTau1") ? jet.userFloat("ak8PFJetsCHSValueMap:NjettinessAK8CHSTau1") : -999;
     output.tau2_           = jet.hasUserFloat("ak8PFJetsCHSValueMap:NjettinessAK8CHSTau2") ? jet.userFloat("ak8PFJetsCHSValueMap:NjettinessAK8CHSTau2") : -999;
     output.tau3_           = jet.hasUserFloat("ak8PFJetsCHSValueMap:NjettinessAK8CHSTau3") ? jet.userFloat("ak8PFJetsCHSValueMap:NjettinessAK8CHSTau3") : -999;
     
     output.prunedMass_     = jet.hasUserFloat("ak8PFJetsCHSValueMap:ak8PFJetsCHSPrunedMass")   ? jet.userFloat("ak8PFJetsCHSValueMap:ak8PFJetsCHSPrunedMass")   : -999;
     output.softDropMass_   = jet.hasUserFloat("ak8PFJetsCHSValueMap:ak8PFJetsCHSSoftDropMass") ? jet.userFloat("ak8PFJetsCHSValueMap:ak8PFJetsCHSSoftDropMass") : -999;
     
     output.puppiTau1_      = jet.hasUserFloat("NjettinessAK8Puppi:tau1") ? jet.userFloat("NjettinessAK8Puppi:tau1") : -999;
     output.puppiTau2_      = jet.hasUserFloat("NjettinessAK8Puppi:tau2") ? jet.userFloat("NjettinessAK8Puppi:tau2") : -999;
     output.puppiTau3_      = jet.hasUserFloat("NjettinessAK8Puppi:tau3") ? jet.userFloat("NjettinessAK8Puppi:tau3") : -999;
     output.puppiMass_      = jet.hasUserFloat("ak8PFJetsPuppiSoftDropMass") ? jet.userFloat("ak8PFJetsPuppiSoftDropMass") : -999;
  } 
  else if(setup_ == 2018){
    output.tau1_           = jet.hasUserFloat("ak8PFJetsCHSValueMap:NjettinessAK8CHSTau1") ? jet.userFloat("ak8PFJetsCHSValueMap:NjettinessAK8CHSTau1") : -999;
    output.tau2_           = jet.hasUserFloat("ak8PFJetsCHSValueMap:NjettinessAK8CHSTau2") ? jet.userFloat("ak8PFJetsCHSValueMap:NjettinessAK8CHSTau2") : -999;
    output.tau3_           = jet.hasUserFloat("ak8PFJetsCHSValueMap:NjettinessAK8CHSTau3") ? jet.userFloat("ak8PFJetsCHSValueMap:NjettinessAK8CHSTau3") : -999;
    output.corrPrunedMass_ = jet.hasUserFloat("ak8PFJetsCHSCorrPrunedMass") ? jet.userFloat("ak8PFJetsCHSCorrPrunedMass") : -999;
    
    output.prunedMass_     = jet.hasUserFloat("ak8PFJetsCHSValueMap:ak8PFJetsCHSPrunedMass")   ? jet.userFloat("ak8PFJetsCHSValueMap:ak8PFJetsCHSPrunedMass")   : -999;
    output.softDropMass_   = jet.hasUserFloat("ak8PFJetsCHSValueMap:ak8PFJetsCHSSoftDropMass") ? jet.userFloat("ak8PFJetsCHSValueMap:ak8PFJetsCHSSoftDropMass") : -999;
    
    output.puppiTau1_      = jet.hasUserFloat("NjettinessAK8Puppi:tau1") ? jet.userFloat("NjettinessAK8Puppi:tau1") : -999;
    output.puppiTau2_      = jet.hasUserFloat("NjettinessAK8Puppi:tau2") ? jet.userFloat("NjettinessAK8Puppi:tau2") : -999;
    output.puppiTau3_      = jet.hasUserFloat("NjettinessAK8Puppi:tau3") ? jet.userFloat("NjettinessAK8Puppi:tau3") : -999;
    output.puppiMass_      = jet.hasUserFloat("ak8PFJetsPuppiSoftDropMass") ? jet.userFloat("ak8PFJetsPuppiSoftDropMass") : -999;
  }
  else {
    edm::LogError("TreePlanter") << "Do not know what to do with the year " << setup_ << endl;
  }




  return output; 
}

template<typename T, typename PAR>
phys::Boson<PAR> TreePlanter::fillBoson(const pat::CompositeCandidate & v, int type, bool requireQualityCriteria) const {

  // Filter on the quality of daughters, right now only for ll pairs
  if(requireQualityCriteria && v.hasUserFloat("GoodLeptons") && !v.userFloat("GoodLeptons")) return phys::Boson<PAR>();

  PAR d0 = fill(*dynamic_cast<const T*>(v.daughter(0)->masterClone().get()));
  PAR d1 = fill(*dynamic_cast<const T*>(v.daughter(1)->masterClone().get()));
  
    if(d0.id() == 0 || d1.id() == 0) edm::LogError("TreePlanter") << "TreePlanter: VB candidate does not have a matching good particle!";
  
  phys::Boson<PAR> physV(d0, d1, type);
  
  // Add FSR
  for(unsigned int i = 2; i < v.numberOfDaughters(); ++i){
    phys::Particle photon(phys::Particle::convert(v.daughter(i)->p4()), 0, 22);
    const pat::PFParticle* fsr = dynamic_cast<const pat::PFParticle*>(v.daughter(i));
    physV.addFSR(fsr->userFloat("leptIdx"), photon);
  }
  
  return physV;
}


// Right now, only for ZZ
template<typename T1, typename T2>
phys::DiBoson<phys::Lepton,phys::Lepton> TreePlanter::fillDiBoson(const pat::CompositeCandidate& edmVV) const{

  int regionWord = computeRegionFlag(edmVV);
  Channel channel = NONE;
  if(test_bit(regionWord,ZZ) || test_bit(regionWord,ZZOnShell)) channel = ZZ;
  else if(test_bit(regionWord,CRZLLos_2P2F) || test_bit(regionWord,CRZLLos_3P1F) || test_bit(regionWord,CRZLLos_2P2F_ZZOnShell) || test_bit(regionWord,CRZLLos_3P1F_ZZOnShell)) channel = ZLL;

  if(channel == NONE) {cout << "Channel cannot be identified, aborting..." << endl; abort();}
  
  const pat::CompositeCandidate* edmV0   = dynamic_cast<const pat::CompositeCandidate*>(edmVV.daughter("Z1")->masterClone().get());      
  const pat::CompositeCandidate* edmV1   = dynamic_cast<const pat::CompositeCandidate*>(edmVV.daughter("Z2")->masterClone().get());
  
  phys::Boson<phys::Lepton> V0;
  phys::Boson<phys::Lepton> V1;
   

  // The first boson is always a good Z, also in the CR. For the other particle assign 23 if it is a true Z from SR
  // or 26 if the two additional leptons comes from LL.

  int idV1 = channel != ZLL ? 23 : 26;

  if(dynamic_cast<const T1*>(edmV0->daughter(0)->masterClone().get()) && dynamic_cast<const T2*>(edmV1->daughter(0)->masterClone().get())){
    V0 = fillBoson<T1,phys::Lepton>(*edmV0, 23, false);
    V1 = fillBoson<T2,phys::Lepton>(*edmV1, idV1, false);
  }
  else if(dynamic_cast<const T2*>(edmV0->daughter(0)->masterClone().get()) && dynamic_cast<const T1*>(edmV1->daughter(0)->masterClone().get())){
    V0 = fillBoson<T1,phys::Lepton>(*edmV1, idV1, false);
    V1 = fillBoson<T2,phys::Lepton>(*edmV0, 23, false);
  }
  else{
    edm::LogError("TreePlanter") << "Do not know what to cast in fillDiBosons" << endl;
    return phys::DiBoson<phys::Lepton,phys::Lepton>();
  }
  
  phys::DiBoson<phys::Lepton,phys::Lepton> VV(V0, V1);

  VV.isBestCand_                = edmVV.hasUserFloat("isBestCand")             ? edmVV.userFloat("isBestCand")             : false;
  VV.passFullSel_               = edmVV.hasUserFloat("SR")                     ? edmVV.userFloat("SR")                     : false;
  VV.isBestCRZLLos_2P2F_        = edmVV.hasUserFloat("isBestCRZLLos_2P2F")     ? edmVV.userFloat("isBestCRZLLos_2P2F")     : false;
  VV.passSelZLL_2P2F_           = edmVV.hasUserFloat("CRZLLos_2P2F")           ? edmVV.userFloat("CRZLLos_2P2F")           : false;
  VV.isBestCRZLLos_3P1F_        = edmVV.hasUserFloat("isBestCRZLLos_3P1F")     ? edmVV.userFloat("isBestCRZLLos_3P1F")     : false;
  VV.passSelZLL_3P1F_           = edmVV.hasUserFloat("CRZLLos_3P1F")           ? edmVV.userFloat("CRZLLos_3P1F")           : false;
  VV.passSRZZOnShell_           = edmVV.hasUserFloat("SR_ZZOnShell")           ? edmVV.userFloat("SR_ZZOnShell")           : false;
  VV.passSelZLL_2P2F_ZZOnShell_ = edmVV.hasUserFloat("CRZLLos_2P2F_ZZOnShell") ? edmVV.userFloat("CRZLLos_2P2F_ZZOnShell") : false;
  VV.passSelZLL_3P1F_ZZOnShell_ = edmVV.hasUserFloat("CRZLLos_3P1F_ZZOnShell") ? edmVV.userFloat("CRZLLos_3P1F_ZZOnShell") : false;   
  VV.regionWord_  = regionWord;
  VV.triggerWord_ = triggerWord_;
  VV.passTrigger_ = filterController_.passTrigger(channel, triggerWord_); // triggerWord_ needs to be filled beforehand (as it is).

  // MELA info
  MELA_.p_JJVBF_BKG_MCFM_JECNominal_ = edmVV.userFloat("p_JJVBF_BKG_MCFM_JECNominal");
  MELA_.p_JJVBF_BKG_MCFM_JECUp_      = edmVV.userFloat("p_JJVBF_BKG_MCFM_JERUp"     );
  MELA_.p_JJVBF_BKG_MCFM_JECDn_      = edmVV.userFloat("p_JJVBF_BKG_MCFM_JERDn"     );

  MELA_.p_JJQCD_BKG_MCFM_JECNominal_ = edmVV.userFloat("p_JJQCD_BKG_MCFM_JECNominal");
  MELA_.p_JJQCD_BKG_MCFM_JECUp_      = edmVV.userFloat("p_JJQCD_BKG_MCFM_JERUp"     );
  MELA_.p_JJQCD_BKG_MCFM_JECDn_      = edmVV.userFloat("p_JJQCD_BKG_MCFM_JERDn"     );

  MELA_.p_JJEW_BKG_MCFM_JECNominal_  = edmVV.userFloat("p_JJEW_BKG_MCFM_JECNominal" );
  MELA_.p_JJEW_BKG_MCFM_JECUp_       = edmVV.userFloat("p_JJEW_BKG_MCFM_JERUp"      );
  MELA_.p_JJEW_BKG_MCFM_JECDn_       = edmVV.userFloat("p_JJEW_BKG_MCFM_JERDn"      );


  return VV;
} //filldibosons end



std::vector<phys::Boson<phys::Jet> > TreePlanter::fillHadBosons(const edm::Handle<edm::View<pat::CompositeCandidate> > & edmBosons, int type) const {
  
  std::vector<phys::Boson<phys::Jet> > physBosons;
  
  foreach(const pat::CompositeCandidate& v, *edmBosons){
    phys::Boson<phys::Jet> physV = fillBoson<pat::Jet, phys::Jet> (v, type, true);
    if(physV.isValid()) physBosons.push_back(physV);
  }
  return physBosons;
}


std::vector<phys::Boson<phys::Lepton> > TreePlanter::fillLepBosons(const edm::Handle<edm::View<pat::CompositeCandidate> > & edmBosons, int type) const {
  
  std::vector<phys::Boson<phys::Lepton> > physBosons;
  
  foreach(const pat::CompositeCandidate& v, *edmBosons){
    phys::Boson<phys::Lepton> physV;
    int rawchannel = v.daughter(0)->pdgId() * v.daughter(1)->pdgId();
    if      (rawchannel == -pow(11,2)) physV = fillBoson<pat::Electron, phys::Lepton>(v, type, true);
    else if (rawchannel == -pow(13,2)) physV = fillBoson<pat::Muon,     phys::Lepton>(v, type, true);
    else {cout << "TreePlanter: unexpected boson final state: " << rawchannel << " ... going to abort.. " << endl; abort();}
    
    if(physV.isValid()) physBosons.push_back(physV);
  }
  return physBosons;
}



// Right now, only for ZZ
std::vector<phys::DiBoson<phys::Lepton,phys::Lepton> > TreePlanter::fillDiBosons(const edm::Handle<edm::View<pat::CompositeCandidate> > & edmDiBosons) const{

  std::vector<phys::DiBoson<phys::Lepton,phys::Lepton> > physDiBosons;
  
  foreach(const pat::CompositeCandidate& edmVV, *edmDiBosons){
    phys::DiBoson<phys::Lepton,phys::Lepton> physVV;
    
    int finalStateZ1 = abs(edmVV.daughter(0)->daughter(0)->pdgId())+abs(edmVV.daughter(0)->daughter(1)->pdgId());
    int finalStateZ2 = abs(edmVV.daughter(1)->daughter(0)->pdgId())+abs(edmVV.daughter(1)->daughter(1)->pdgId());

    int rawchannel = finalStateZ1+finalStateZ2;
    
    if     (rawchannel == 52)         physVV = fillDiBoson<pat::Muon,     pat::Muon    >(edmVV);
    else if(rawchannel == 44)         physVV = fillDiBoson<pat::Electron, pat::Electron>(edmVV);
    else if(rawchannel == 48)  {
      if(finalStateZ1 < finalStateZ2) physVV = fillDiBoson<pat::Electron, pat::Muon    >(edmVV);
      else                            physVV = fillDiBoson<pat::Muon, pat::Electron    >(edmVV);
    }
    else {cout << "TreePlanter: unexpected diboson final state: " << rawchannel << " ... going to abort.. " << endl; abort();}
    
    if(physVV.isValid()) physDiBosons.push_back(physVV);
    
  }


  return physDiBosons;
}

std::vector<std::pair<phys::Boson<phys::Lepton>, phys::Lepton> > TreePlanter::fillZLCandidates(const edm::Handle<edm::View<pat::CompositeCandidate> > & edmZLs) const{

  std::vector<std::pair<phys::Boson<phys::Lepton>, phys::Lepton> > physZLs;

  if(edmZLs->size() != 1) return physZLs; // FIXME: make physZLs an obj not a container.
  if(!filterController_.passTrigger(ZL, triggerWord_)) return physZLs;

  foreach(const pat::CompositeCandidate& edmZL, *edmZLs){
    
    const pat::CompositeCandidate* edmV0 = dynamic_cast<const pat::CompositeCandidate*>(edmZL.daughter(0)->masterClone().get());      
    phys::Boson<phys::Lepton> V0;
    
    if     (abs(edmV0->daughter(0)->pdgId()) == 11) V0 = fillBoson<pat::Electron,phys::Lepton>(*edmV0, 23, false);
    else if(abs(edmV0->daughter(0)->pdgId()) == 13) V0 = fillBoson<pat::Muon,phys::Lepton>(*edmV0, 23, false);
    else{
      edm::LogError("TreePlanter") << "Do not know what to cast in fillZLCandidates, Z part" << endl;
      abort();
    }
    
    phys::Lepton lep;
    if     (abs(edmZL.daughter(1)->pdgId()) == 11)
      lep = fill(*dynamic_cast<const pat::Electron*>(edmZL.daughter(1)->masterClone().get()));
    else if(abs(edmZL.daughter(1)->pdgId()) == 13)
      lep = fill(*dynamic_cast<const pat::Muon*>(edmZL.daughter(1)->masterClone().get()));
    else{
      edm::LogError("TreePlanter") << "Do not know what to cast in fillZLCandidates, LEP part" << endl;
      abort();
    }

    if(V0.isValid() && lep.isValid()) physZLs.push_back(std::make_pair(V0,lep));
  }

  return physZLs;  
}



phys::DiBoson<phys::Lepton,phys::Lepton> 
TreePlanter::fillZWCandidate(const edm::Handle<edm::View<pat::CompositeCandidate> > & edmZWs) const{

  // Only one ZW is allowed. ZL_ is empty when edmZL.size() != 1. 
  if(edmZWs->size() != 1 || ZL_.empty()) return  phys::DiBoson<phys::Lepton,phys::Lepton>();
  //if(!filterController_.passTrigger(ZL, triggerWord_)) return physZLs;


  // for ZW,  daughter(0) is the ZL while daughter(1) is the MET
  const pat::CompositeCandidate edmZW = edmZWs->front();

  // ----- Build the Z boson -----
  const pat::CompositeCandidate* edmZ  = dynamic_cast<const pat::CompositeCandidate*>(edmZW.daughter(0)->daughter(0)->masterClone().get());
  phys::Boson<phys::Lepton> physZ;
    
  if     (abs(edmZ->daughter(0)->pdgId()) == 11) physZ = fillBoson<pat::Electron, phys::Lepton>(*edmZ, 23, false);
  else if(abs(edmZ->daughter(0)->pdgId()) == 13) physZ = fillBoson<pat::Muon    , phys::Lepton>(*edmZ, 23, false);
  else{
    edm::LogError("TreePlanter") << "Do not know what to cast in fillZWCandidates, Z part" << endl;
    abort();
  }
  // -----------------------------

  // ----- Build the W boson -----
  phys::Lepton lep;
  if     (abs(edmZW.daughter(0)->daughter(1)->pdgId()) == 11)
    lep = fill(*dynamic_cast<const pat::Electron*>(edmZW.daughter(0)->daughter(1)->masterClone().get()));
  else if(abs(edmZW.daughter(0)->daughter(1)->pdgId()) == 13)
    lep = fill(*dynamic_cast<const pat::Muon*>    (edmZW.daughter(0)->daughter(1)->masterClone().get()));
  else{
    edm::LogError("TreePlanter") << "Do not know what to cast in fillZWCandidates, LEP part" << endl;
    abort();
  }

  phys::Lepton met = phys::Lepton(phys::Particle::convert(edmZW.daughter(1)->p4()));

  phys::Boson<phys::Lepton> physW(lep, met, copysign(24,-1*lep.id()));
  // -----------------------------

  
  // ----- Build the ZW -----
  if(physZ.isValid() && physW.isValid() ){
    phys::DiBoson<phys::Lepton,phys::Lepton> ZW(physZ, physW);
    
    int regionWord = 0;
    set_bit(regionWord,30); // use 30 for WZ
    
    ZW.isBestCand_                = true;
    ZW.passFullSel_               = true;
    ZW.isBestCRZLLos_2P2F_        = false;
    ZW.passSelZLL_2P2F_           = false;
    ZW.isBestCRZLLos_3P1F_        = false;
    ZW.passSelZLL_3P1F_           = false;
    ZW.passSRZZOnShell_           = false;
    ZW.passSelZLL_2P2F_ZZOnShell_ = false;
    ZW.passSelZLL_3P1F_ZZOnShell_ = false;
    ZW.regionWord_                = regionWord;
    ZW.triggerWord_               = triggerWord_;
    ZW.passTrigger_               = filterController_.passTrigger(ZZ, triggerWord_); // use same trigger as ZZ
    
    return ZW;
  }
  else return phys::DiBoson<phys::Lepton,phys::Lepton>();
}








int TreePlanter::computeRegionFlag(const pat::CompositeCandidate & vv) const{
  int REGIONFLAG=0;
  
  if(vv.hasUserFloat("isBestCand")         && vv.hasUserFloat("SR")                     && vv.userFloat("isBestCand")         && vv.userFloat("SR"))
    set_bit(REGIONFLAG,ZZ);
  if(vv.hasUserFloat("isBestCRZLLos_2P2F") && vv.hasUserFloat("CRZLLos_2P2F")           && vv.userFloat("isBestCRZLLos_2P2F") && vv.userFloat("CRZLLos_2P2F"))
    set_bit(REGIONFLAG,CRZLLos_2P2F);
  if(vv.hasUserFloat("isBestCRZLLos_3P1F") && vv.hasUserFloat("CRZLLos_3P1F")           && vv.userFloat("isBestCRZLLos_3P1F") && vv.userFloat("CRZLLos_3P1F"))
    set_bit(REGIONFLAG,CRZLLos_3P1F);
  if(vv.hasUserFloat("isBestCand")         && vv.hasUserFloat("SR_ZZOnShell")           && vv.userFloat("isBestCand")         && vv.userFloat("SR_ZZOnShell"))
    set_bit(REGIONFLAG,ZZOnShell);
  if(vv.hasUserFloat("isBestCRZLLos_2P2F") && vv.hasUserFloat("CRZLLos_2P2F_ZZOnShell") && vv.userFloat("isBestCRZLLos_2P2F") && vv.userFloat("CRZLLos_2P2F_ZZOnShell"))
    set_bit(REGIONFLAG,CRZLLos_2P2F_ZZOnShell);
  if(vv.hasUserFloat("isBestCRZLLos_3P1F") && vv.hasUserFloat("CRZLLos_3P1F_ZZOnShell") && vv.userFloat("isBestCRZLLos_3P1F") && vv.userFloat("CRZLLos_3P1F_ZZOnShell"))
    set_bit(REGIONFLAG,CRZLLos_3P1F_ZZOnShell);


  return REGIONFLAG;
}

#include "FWCore/Framework/interface/MakerMacros.h"
// ---- define this as a plug-in ----------------------------------------
DEFINE_FWK_MODULE(TreePlanter);
