#include "VVXAnalysis/Producers/plugins/TreePlanter.h"

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
#include "DataFormats/Common/interface/MergeableCounter.h"
#include "DataFormats/Common/interface/Ptr.h"
#include "DataFormats/JetReco/interface/GenJet.h"

#include "CondFormats/JetMETObjects/interface/JetCorrectorParameters.h"
#include "CondFormats/JetMETObjects/interface/SimpleJetResolution.h"

#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"
#include "SimDataFormats/GeneratorProducts/interface/GenRunInfoProduct.h"

#include "AnalysisDataFormats/CMGTools/interface/PFJet.h"
#include "AnalysisDataFormats/CMGTools/interface/BaseMET.h"
#include "VVXAnalysis/Commons/interface/Utilities.h"

#include "ZZAnalysis/AnalysisStep/interface/MCHistoryTools.h"
#include "ZZAnalysis/AnalysisStep/interface/bitops.h"

#include "TTree.h"
#include "TVectorD.h"

#include <boost/foreach.hpp>
#define foreach BOOST_FOREACH

#include <boost/assign/std/vector.hpp>
using namespace boost::assign;

using std::cout;
using std::endl;


TreePlanter::TreePlanter(const edm::ParameterSet &config)
  : PUWeighter_      (PUReweight::LEGACY)
  , filterController_(config)
  , mcHistoryTools_  (0)
  , leptonScaleFactors_(edm::FileInPath("ZZAnalysis/AnalysisStep/test/Macros/scale_factors_muons2012.root").fullPath(),
			edm::FileInPath("ZZAnalysis/AnalysisStep/test/Macros/scale_factors_ele2012.root").fullPath(),
			edm::FileInPath("VVXAnalysis/Commons/data/fakeRates_mu.root").fullPath(),
			edm::FileInPath("VVXAnalysis/Commons/data/fakeRates_el.root").fullPath())
  , signalDefinition_(config.getParameter<int>("signalDefinition"   ))
  , passTrigger_(false)
  , passSkim_(false)
  , triggerWord_(0)
  , preSkimCounter_  (0)
  , postSkimCounter_ (0)
  , postSkimSignalCounter_(0)
  , signalCounter_(0)
  , postSkimSignalEvents_(0)
  , theMuonLabel     (config.getParameter<edm::InputTag>("muons"    ))
  , theElectronLabel (config.getParameter<edm::InputTag>("electrons"))
  , theJetLabel      (config.getParameter<edm::InputTag>("jets"     ))
  , theVhadLabel     (config.getParameter<edm::InputTag>("Vhad"     ))
  , theZZLabel       (config.getParameter<edm::InputTag>("ZZ"       ))
  , theZLLabel       (config.getParameter<edm::InputTag>("ZL"       ))
  , theMETLabel      (config.getParameter<edm::InputTag>("MET"      ))
  , theVertexLabel   (config.getParameter<edm::InputTag>("Vertices" ))
  , sampleName_      (config.getParameter<std::string>("sampleName"))
  , isMC_            (config.getUntrackedParameter<bool>("isMC",false))
  , sampleType_      (config.getParameter<int>("sampleType"))
  , setup_           (config.getParameter<int>("setup"))
  , applyTrigger_    (config.getUntrackedParameter<bool>("TriggerRequired", false)) 
  , applySkim_       (config.getUntrackedParameter<bool>("SkimRequired"   , true)) 
  , applyMCSel_      (config.getUntrackedParameter<bool>("DoMCSelection"  , false)) 
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
    thePUInfoLabel           = config.getUntrackedParameter<edm::InputTag>("PUInfo"         , edm::InputTag("addPileupInfo"));
    theGenCategoryLabel      = config.getUntrackedParameter<edm::InputTag>("GenCategory"    , edm::InputTag("genCategory"));
    theGenCollectionLabel    = config.getUntrackedParameter<edm::InputTag>("GenCollection"  , edm::InputTag("genParticlesPruned"));
    theGenJetCollectionLabel = config.getUntrackedParameter<edm::InputTag>("GenJets"        , edm::InputTag("genCategory","genJets"));
    theGenVBCollectionLabel  = config.getUntrackedParameter<edm::InputTag>("GenVBCollection", edm::InputTag("genCategory","vectorBosons"));
    externalCrossSection_    = config.getUntrackedParameter<double>("XSection",-1);
  }
   
  skimPaths_ = config.getParameter<std::vector<std::string> >("skimPaths");

  jetCorrectorParameters_ = new JetCorrectorParameters(edm::FileInPath("CondFormats/JetMETObjects/data/JetResolutionInputAK5PFCHS.txt").fullPath());
  JER_ = new  SimpleJetResolution(*jetCorrectorParameters_);

  initTree();
}


void TreePlanter::beginJob(){
  theTree->Branch("event"     , &event_); 
  theTree->Branch("run"       , &run_); 
  theTree->Branch("lumiBlock" , &lumiBlock_); 

  theTree->Branch("passTrigger" , &passTrigger_); 
  theTree->Branch("passSkim"    , &passSkim_); 
  theTree->Branch("triggerWord" , &triggerWord_); 

  theTree->Branch("mcprocweight"     , &mcprocweight_);
  theTree->Branch("puweight"         , &puweight_);
  theTree->Branch("genCategory" , &genCategory_);

  theTree->Branch("met"   , &met_);
  theTree->Branch("nvtxs" , &nvtx_);

  theTree->Branch("muons"     , &muons_);
  theTree->Branch("electrons" , &electrons_);
  theTree->Branch("jets"      , &jets_); 
  theTree->Branch("VhadCand"  , &Vhad_);
  theTree->Branch("ZZCand"    , &ZZ_); 

  theTree->Branch("ZLCand"    , &ZL_); 

  theTree->Branch("genParticles"  , &genParticles_);
  theTree->Branch("genVBParticles", &genVBParticles_);
  theTree->Branch("genJets"       , &genJets_);
}

void TreePlanter::endLuminosityBlock(edm::LuminosityBlock const& lumi, edm::EventSetup const& setup)
{
  // Beware: preSkimCounter for H->ZZ->4l means a skim done at path level

  Float_t Nevt_preskim = -1.;
  edm::Handle<edm::MergeableCounter> preSkimCounter;
  if (lumi.getByLabel("preSkimCounter", preSkimCounter)) { // Counter before skim. Does not exist for non-skimmed samples.
    Nevt_preskim = preSkimCounter->value;
  }  
  
  edm::Handle<edm::MergeableCounter> prePathCounter;
  lumi.getByLabel("prePathCounter", prePathCounter);       // Counter of input events in the input pattuple

  // Nevt_gen: this is the number before any skim
  if (Nevt_preskim>=0.) {
    theNumberOfEvents += Nevt_preskim; 
  } else {
    theNumberOfEvents += prePathCounter->value;    
  }

  // Beware: pre/post Skim here means before/after the preselection path at Tree building level!
  edm::Handle<edm::MergeableCounter> counter;
  bool found = lumi.getByLabel("prePreselectionCounter", counter);
  if(found) preSkimCounter_ += counter->value;
  
  found = lumi.getByLabel("postPreselectionCounter", counter);
  if(found) postSkimCounter_ += counter->value;

  found = lumi.getByLabel("signalCounter", counter);
  if(found) signalCounter_ += counter->value;

  found = lumi.getByLabel("postSkimSignalCounter", counter);
  if(found) postSkimSignalCounter_ += counter->value;

  found = lumi.getByLabel("srCounter", counter);
  if(found) eventsInSR_ += counter->value;

  found = lumi.getByLabel("cr2P2FCounter", counter);
  if(found) eventsIn2P2FCR_ += counter->value;

  found = lumi.getByLabel("cr3P1FCounter", counter);
  if(found) eventsIn3P1FCR_ += counter->value;
}


void TreePlanter::endRun(const edm::Run& run, const edm::EventSetup& setup){
  if(isMC_){
    edm::Handle<GenRunInfoProduct> genRunInfo;
    run.getByLabel("generator", genRunInfo);
    theXSections.push_back(genRunInfo->crossSection());
  }
}

void TreePlanter::endJob(){

  if(jetCorrectorParameters_) delete jetCorrectorParameters_;
  if(JER_) delete JER_;

  if(mcHistoryTools_) delete mcHistoryTools_;

  edm::Service<TFileService> fs;
  TTree *countTree = fs->make<TTree>("HollyTree","HollyTree");
  
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

  mcprocweight_       = 1.;
  puweight_           = 1.; 

  genCategory_    = -1;
  nobservedPUInt_ = -1;
  ntruePUInt_     = -1;

  met_    = phys::Particle();
  nvtx_   = -1;

  muons_     = std::vector<phys::Lepton>();
  electrons_ = std::vector<phys::Electron>();
  jets_      = std::vector<phys::Jet>();
  ZZ_        = phys::DiBoson<phys::Lepton  , phys::Lepton>();
  Vhad_      = std::vector<phys::Boson<phys::Jet>      >();   

  ZL_        = std::vector<std::pair<phys::Boson<phys::Lepton>, phys::Lepton> >();

  genParticles_ = std::vector<phys::Particle>();
  genVBParticles_ = std::vector<phys::Boson<phys::Particle> >();
  genJets_ = std::vector<phys::Particle>();
 }


bool TreePlanter::fillEventInfo(const edm::Event& event){

  // Fill some info abut acceptance before cutting away events. Beware: if the signal is defined a-posteriori, we will have a problem. For that case, we need to
  // explicitly check here that we are counting signal and not irreducible background.
  if (isMC_) {
    // Apply MC filter
    if (!filterController_.passMCFilter(event)) return false;

    if(mcHistoryTools_) delete mcHistoryTools_;
    mcHistoryTools_ = new MCHistoryTools(event, sampleName_);
    bool gen_ZZ4lInEtaAcceptance   = false; // All 4 gen leptons in eta acceptance
    bool gen_ZZ4lInEtaPtAcceptance = false; // All 4 gen leptons in eta,pT acceptance
    bool gen_m4l_180               = false; // gen_m4l > 180
    bool gen_ZZInAcceptance        = false; // Unused; old ZZ phase space
    mcHistoryTools_->genAcceptance(gen_ZZInAcceptance, gen_ZZ4lInEtaAcceptance, gen_ZZ4lInEtaPtAcceptance, gen_m4l_180);
    if (gen_ZZ4lInEtaAcceptance)   ++eventsInEtaAcceptance_;  
    if (gen_ZZ4lInEtaPtAcceptance) ++eventsInEtaPtAcceptance_;


    edm::Handle<std::vector<PileupSummaryInfo> > puInfo; event.getByLabel(thePUInfoLabel, puInfo);
    foreach(const PileupSummaryInfo& pui, *puInfo)
      if(pui.getBunchCrossing() == 0) { 
	nobservedPUInt_  = pui.getPU_NumInteractions();
	ntruePUInt_      = pui.getTrueNumInteractions();
	break;
      }
    

    // Info about the MC weight
    puweight_ = PUWeighter_.weight(sampleType_, setup_, ntruePUInt_);
    
    mcprocweight_ = mcHistoryTools_->gethepMCweight();
    
    // Sum of weight, particularly imprtant for MCs that return also negative weights
    // or, in general, weighted events
    sumpuweights_       += puweight_;
    summcprocweights_   += mcprocweight_;
    sumpumcprocweights_ += puweight_*mcprocweight_;
  }
    
  // Check trigger request. Actually, it is a very very loose request, not the actual one, that instead should be
  // asked to the specific final state
  passTrigger_ = filterController_.passTrigger(NONE, event, triggerWord_);
  if (applyTrigger_ && !passTrigger_) return false;
  
  // Check Skim requests
  passSkim_ = filterController_.passSkim(event, triggerWord_);
  if (applySkim_    && !passSkim_)   return false;

  run_       = event.id().run();
  event_     = event.id().event(); 
  lumiBlock_ = event.luminosityBlock();

  edm::Handle<std::vector<cmg::BaseMET> > met;   event.getByLabel(theMETLabel, met);
  met_ = phys::Particle(phys::Particle::convert(met->front().p4()));

  edm::Handle<std::vector<reco::Vertex> > vertices; event.getByLabel(theVertexLabel, vertices);
  nvtx_ = vertices->size();
    
  if(isMC_){
    
    edm::Handle<edm::View<reco::Candidate> > genParticles;
    event.getByLabel(theGenCollectionLabel,  genParticles);
    
    for (edm::View<reco::Candidate>::const_iterator p = genParticles->begin(); p != genParticles->end(); ++p) 
      if (p->status() == 3 && std::distance(genParticles->begin(),p) > 5)
	genParticles_.push_back(phys::Particle(p->p4(), phys::Particle::computeCharge(p->pdgId()), p->pdgId()));
         
  
    edm::Handle<int> genCategory;
    event.getByLabel(theGenCategoryLabel, genCategory);
    genCategory_ = *genCategory;
    if(((genCategory_ ^ signalDefinition_) & signalDefinition_) == 0) ++postSkimSignalEvents_;

    event.getByLabel(theGenVBCollectionLabel,  genParticles);
    
    for(edm::View<reco::Candidate>::const_iterator p = genParticles->begin(); p != genParticles->end(); ++p)
      if(fabs(p->pdgId()) == 24 || p->pdgId() == 23)
	genVBParticles_.push_back(phys::Boson<phys::Particle>(phys::Particle(p->daughter(0)->p4(), phys::Particle::computeCharge(p->daughter(0)->pdgId()), p->daughter(0)->pdgId()),
							      phys::Particle(p->daughter(1)->p4(), phys::Particle::computeCharge(p->daughter(1)->pdgId()), p->daughter(1)->pdgId()),
							      p->pdgId()));
    
    // Get the gen jet collection
    edm::Handle<edm::View<reco::Candidate> > genJets;
    event.getByLabel(theGenJetCollectionLabel,  genJets);

    for(edm::View<reco::Candidate>::const_iterator jet = genJets->begin(); jet != genJets->end(); ++jet)
      genJets_.push_back(phys::Particle(jet->p4(), phys::Particle::computeCharge(jet->pdgId()), jet->pdgId()));
  }

  return true;
}



void TreePlanter::analyze(const edm::Event& event, const edm::EventSetup& setup){

  initTree();

  bool goodEvent = fillEventInfo(event);
  if(!goodEvent) return;
  ++theNumberOfAnalyzedEvents;

  //// For Z+L CRs, we want only events with exactly 1 Z+l candidate.
  ////if (filterController_.channel() == ZL && ???size() != 1) return;


  // Load a bunch of objects from the event
  edm::Handle<pat::MuonCollection>       muons            ; event.getByLabel(theMuonLabel    ,     muons);
  edm::Handle<pat::ElectronCollection>   electrons        ; event.getByLabel(theElectronLabel, electrons);
  edm::Handle<std::vector<cmg::PFJet> >  jets             ; event.getByLabel(theJetLabel     ,      jets);
  edm::Handle<edm::View<pat::CompositeCandidate> > Vhad   ; event.getByLabel(theVhadLabel    ,      Vhad);
  edm::Handle<edm::View<pat::CompositeCandidate> > ZZ     ; event.getByLabel(theZZLabel      ,        ZZ);
  edm::Handle<edm::View<pat::CompositeCandidate> > ZL     ; event.getByLabel(theZLLabel      ,        ZL);

  foreach(const pat::Muon& muon, *muons){
    //if(!muon.userFloat("isGood") || muon.userFloat("CombRelIsoPF") >= 4) continue;  // commented because the combination of the two flags is more restrictive than Z.userfloat("goodLeptons"), hence the matching can fail.
    if(!muon.userFloat("isGood")) continue; 
    phys::Lepton physmuon = fill(muon);
    muons_.push_back(physmuon);
  }

  foreach(const pat::Electron& electron, *electrons){
    //if(!electron.userFloat("isGood") || electron.userFloat("CombRelIsoPF") >= 4) continue; 
    if(!electron.userFloat("isGood")) continue; 
    phys::Electron physelectron =  fill(electron);
    electrons_.push_back(physelectron);
  }

  foreach(const cmg::PFJet& jet, *jets){
    phys::Jet physjet = fill(jet);
    jets_.push_back(physjet);
  }

  // The bosons are selected requiring that their daughters pass the quality criteria to be good daughters
  Vhad_   = fillBosons<cmg::PFJet,phys::Jet>(Vhad, 24);

  // The bosons have NOT any requirement on the quality of their daughters, only the flag is set (because of the same code is usd for CR too)
  std::vector<phys::DiBoson<phys::Lepton,phys::Lepton> > ZZs = fillDiBosons(ZZ);

  // Fill Z+l pairs for fake rate measurements
  ZL_ = fillZLCandidates(ZL);

  if(ZZ->size() > 1) {
    cout << "----------------------------------------------------" << endl;
    cout << "More than one ZZ candidate!! " << ZZ->size() << endl;  
    cout << "Event: " << event_ << endl;
    typedef phys::DiBoson<phys::Lepton,phys::Lepton> ZZlep;
    foreach(const ZZlep& zz , ZZs){
      cout << "....................." << endl;
      cout << zz << " SR? " << test_bit(zz.regionWord_,Channel::ZZ) << " CR2P2F? " << test_bit(zz.regionWord_,CRZLLos_2P2F) << " CR3P1F? " << test_bit(zz.regionWord_,CRZLLos_3P1F) << endl;
      cout << "daughter 0: "   << zz.first() << endl;
      cout << "daughter 0.1: " << zz.first().daughter(0) << " is good? " <<  zz.first().daughter(0).isGood() << " pass full sel? " << zz.first().daughter(0).passFullSel() <<  endl;
      cout << "daughter 0.1: " << zz.first().daughter(1) << " is good? " <<  zz.first().daughter(1).isGood() << " pass full sel? " << zz.first().daughter(1).passFullSel() <<  endl;
      cout << "daughter 1: "   << zz.second() << endl;
      cout << "daughter 1.1: " << zz.second().daughter(0) << " is good? " <<  zz.second().daughter(0).isGood() << " pass full sel? " << zz.second().daughter(0).passFullSel() <<  endl;
      cout << "daughter 1.1: " << zz.second().daughter(1) << " is good? " <<  zz.second().daughter(1).isGood() << " pass full sel? " << zz.second().daughter(1).passFullSel() <<  endl;
      cout << "....................." << endl;
    }
    cout << "----------------------------------------------------" << endl;
  }
  if(ZZs.size() == 1 && ZZs.front().passTrigger()) ZZ_ = ZZs.front();     
  else if(ZL_.empty() && applySkim_) return;

  theTree->Fill();
}


template<typename LEP>
phys::Lepton TreePlanter::fillLepton(const LEP& lepton) const{

  phys::Lepton output(phys::Particle::convert(lepton.p4()),lepton.charge(),lepton.pdgId());
  
  output.dxy_             = lepton.userFloat("dxy"              );               
  output.dz_              = lepton.userFloat("dz"               );                
  output.sip_             = lepton.userFloat("SIP"              );
  output.combRelIso_      = lepton.userFloat("combRelIso"       );
  output.pfChargedHadIso_ = lepton.userFloat("PFChargedHadIso"  );
  output.pfNeutralHadIso_ = lepton.userFloat("PFNeutralHadIso"  );
  output.pfPhotonIso_     = lepton.userFloat("PFPhotonIso"      );
  output.pfCombRelIso_    = lepton.userFloat("combRelIsoPF"     );  
  output.rho_             = lepton.userFloat("rho"              );
  output.isPF_            = lepton.userFloat("isPFMuon"         );
  output.matchHLT_        = lepton.userFloat("HLTMatch"         );
  output.isGood_          = lepton.userFloat("isGood"           );
  output.efficiencySF_    = leptonScaleFactors_.efficiencyScaleFactor(output);
  output.setFakeRateSF(leptonScaleFactors_.fakeRateScaleFactor(output));

  return output; 
}


phys::Lepton TreePlanter::fill(const pat::Electron &electron) const{

  return fillLepton(electron);
  // phys::Electron output(fillLepton(electron));
  
  // output.energy_     = electron.userFloat("energy"    );
  // output.phiWidth_   = electron.userFloat("phiWidth"  );
  // output.etaWidth_   = electron.userFloat("etaWidth"  );
  // output.BDT_        = electron.userFloat("BDT"       );
  // output.isBDT_      = electron.userFloat("isBDT"     );
  // output.missingHit_ = electron.userInt  ("missingHit");
  // output.nCrystals_  = electron.userInt  ("nCrystals" );

  //return output;
}


phys::Lepton TreePlanter::fill(const pat::Muon& mu) const{
  return fillLepton(mu);
}



phys::Jet TreePlanter::fill(const cmg::PFJet &jet) const{
  
  phys::Jet output(phys::Particle::convert(jet.p4()),jet.charge(),1);
  
  output.nConstituents_ = jet.nConstituents();
  output.nCharged_ = jet.nCharged();
  output.nNeutral_ = jet.nNeutral();
    
  output.neutralHadronEnergyFraction_ = jet.component(reco::PFCandidate::ParticleType::h0).fraction();
  output.chargedHadronEnergyFraction_ = jet.component(reco::PFCandidate::ParticleType::h).fraction();
  output.chargedEmEnergyFraction_     = jet.component(reco::PFCandidate::ParticleType::e).fraction();    
  output.neutralEmEnergyFraction_     = jet.component(reco::PFCandidate::ParticleType::gamma).fraction();    
  output.muonEnergyFraction_          = jet.component(reco::PFCandidate::ParticleType::mu).fraction();
  
  output.csvtagger_               = jet.btag("combinedSecondaryVertexBJetTags");
  output.girth_                   = jet.girth();
  output.girth_charged_           = jet.girth_charged();
  output.ptd_                     = jet.ptd();
  output.rms_                     = jet.rms();
  output.beta_                    = jet.beta();
  output.jetArea_                 = jet.jetArea();
  output.secvtxMass_              = jet.secvtxMass();
  output.Lxy_                     = jet.Lxy();
  output.LxyErr_                  = jet.LxyErr();
  output.rawFactor_               = jet.rawFactor();
  output.uncOnFourVectorScale_    = jet.uncOnFourVectorScale();
  output.puMVAFull_               = jet.puMva("full53x");
  output.puMVASimple_             = jet.puMva("simple");
  output.puCutBased_              = jet.puMva("cut-based");
  output.pass_puMVAFull_loose_    = jet.passPuJetId("full53x"  , PileupJetIdentifier::kLoose);
  output.pass_pUMVAFull_medium_   = jet.passPuJetId("full53x"  , PileupJetIdentifier::kMedium);
  output.pass_pUMVAFull_tight_    = jet.passPuJetId("full53x"  , PileupJetIdentifier::kTight); 
  output.pass_puMVASimple_loose_  = jet.passPuJetId("simple"   , PileupJetIdentifier::kLoose); 
  output.pass_puMVASimple_medium_ = jet.passPuJetId("simple"   , PileupJetIdentifier::kMedium); 
  output.pass_puMVASimple_tight_  = jet.passPuJetId("simple"   , PileupJetIdentifier::kTight); 
  output.pass_puCutBased_loose_   = jet.passPuJetId("cut-based", PileupJetIdentifier::kLoose); 
  output.pass_puCutBased_medium_  = jet.passPuJetId("cut-based", PileupJetIdentifier::kMedium);
  output.pass_puCutBased_tight_   = jet.passPuJetId("cut-based", PileupJetIdentifier::kTight);
  
  output.mcPartonFlavour_ = jet.partonFlavour();

  // Fill JER quantities
  if(isMC_){
    std::vector<float> fx, fY;
    fx.push_back(output.eta()); // Jet Eta
    fY.push_back(output.pt());  // Jet PT
    fY.push_back(ntruePUInt_);  // Number of truth pileup
    output.sigma_MC_ = JER_->resolution(fx,fY);
    
    // From: https://twiki.cern.ch/twiki/bin/view/CMS/JetResolution
    if      (abs(output.eta()) < 0.5) {output.jer_c_ = 1.079; output.jer_cdown_ = 1.053; output.jer_cup_ = 1.105;}
    else if (abs(output.eta()) < 1.1) {output.jer_c_ = 1.099; output.jer_cdown_ = 1.071; output.jer_cup_ = 1.127;}
    else if (abs(output.eta()) < 1.7) {output.jer_c_ = 1.121; output.jer_cdown_ = 1.092; output.jer_cup_ = 1.150;}
    else if (abs(output.eta()) < 2.3) {output.jer_c_ = 1.208; output.jer_cdown_ = 1.162; output.jer_cup_ = 1.254;}
    else if (abs(output.eta()) < 2.8) {output.jer_c_ = 1.254; output.jer_cdown_ = 1.192; output.jer_cup_ = 1.316;}
    else if (abs(output.eta()) < 3.2) {output.jer_c_ = 1.395; output.jer_cdown_ = 1.332; output.jer_cup_ = 1.458;}
    else if (abs(output.eta()) < 5.0) {output.jer_c_ = 1.056; output.jer_cdown_ = 0.865; output.jer_cup_ = 1.247;}
  }

  return output; 
}


void TreePlanter::addExtras(phys::Jet &jet, const pat::CompositeCandidate & v, const std::string& userFloatName) const{}

void TreePlanter::addExtras(phys::Lepton& l, const pat::CompositeCandidate & v, const std::string& userFloatName) const {
  // in the future it could become a map inside the phys::Particle class. Right now there is not a realy need to
  // make such a complication.
  if(v.hasUserFloat(userFloatName)){
    l.pfCombRelIsoFSRCorr_ = v.userFloat(userFloatName);
    l.realignFakeRate();
  }
}


template<typename T, typename PAR>
phys::Boson<PAR> TreePlanter::fillBoson(const pat::CompositeCandidate & v, int type, bool requireQualityCriteria) const {

  // Filter on the quality of daughters, right now only for ll pairs
  if(requireQualityCriteria && v.hasUserFloat("GoodLeptons") && !v.userFloat("GoodLeptons")) return phys::Boson<PAR>();

  PAR d0 = fill(*dynamic_cast<const T*>(v.daughter(0)->masterClone().get()));
  addExtras(d0, v ,"d0.combRelIsoPFFSRCorr");
  PAR d1 = fill(*dynamic_cast<const T*>(v.daughter(1)->masterClone().get()));
  addExtras(d1, v ,"d1.combRelIsoPFFSRCorr");

  if(d0.id() == 0 || d1.id() == 0) edm::LogError("TreePlanter") << "TreePlanter: VB candidate does not have a matching good particle!";
  
  phys::Boson<PAR> physV(d0, d1, type);
  
  // Add FSR
  if(v.hasUserFloat("dauWithFSR") && v.userFloat("dauWithFSR") >= 0){
    phys::Particle photon(phys::Particle::convert(v.daughter(2)->p4()), 0, 22);
    physV.addFSR(v.userFloat("dauWithFSR"), photon);
  }
  
  // Add quality of daughters, right now only for ll pairs
  if(v.hasUserFloat("GoodLeptons")) physV.hasGoodDaughters_ = v.userFloat("GoodLeptons");
  
  return physV;
}


// Right now, only for ZZ
template<typename T1, typename T2>
phys::DiBoson<phys::Lepton,phys::Lepton> TreePlanter::fillDiBoson(const pat::CompositeCandidate& edmVV) const{

  int regionWord = computeRegionFlag(edmVV);

  Channel channel = NONE;
  if(test_bit(regionWord,ZZ)) channel = ZZ;
  else if(test_bit(regionWord,CRZLLos_2P2F) || test_bit(regionWord,CRZLLos_3P1F)) channel = ZLL;
  
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

  // Set the trigger bit. Special care should be put for ZLL CR, for the other channels the code below do not add anything to FilterController:passTrigger
  Channel effectiveChannel = channel;
  if      (VV.id() == 44) effectiveChannel = EEEE;  // ZZ->4e
  else if (VV.id() == 48) effectiveChannel = EEMM;  // ZZ->2e2mu
  else if (VV.id() == 52) effectiveChannel = MMMM;  // ZZ->4mu
  else {cout << "Do not know what to do when setting trigger bit in TreePlanter. Unknown ZZ id: " << VV.id() << endl; abort();}
  
  VV.isBestCand_         = edmVV.userFloat("isBestCand");
  VV.passFullSel_        = edmVV.userFloat("FullSelTight");
  VV.isBestCRZLLos_2P2F_ = edmVV.userFloat("isBestCRZLLos_2P2F");
  VV.passSelZLL_2P2F_    = edmVV.userFloat("SelZLL_2P2F");
  VV.isBestCRZLLos_3P1F_ = edmVV.userFloat("isBestCRZLLos_3P1F");
  VV.passSelZLL_3P1F_    = edmVV.userFloat("SelZLL_3P1F");   

  VV.regionWord_  = regionWord;
  VV.triggerWord_ = triggerWord_;
  VV.passTrigger_ = filterController_.passTrigger(effectiveChannel, triggerWord_); // triggerWord_ needs to be filled beforehand (as it is).
  
  return VV;
}

template<typename T, typename PAR>
std::vector<phys::Boson<PAR> > TreePlanter::fillBosons(const edm::Handle<edm::View<pat::CompositeCandidate> > & edmBosons, int type) const {
  
  std::vector<phys::Boson<PAR> > physBosons;
  
  foreach(const pat::CompositeCandidate& v, *edmBosons){
    phys::Boson<PAR> physV = fillBoson<T, PAR>(v, type, true);
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



int TreePlanter::computeRegionFlag(const pat::CompositeCandidate & vv) const{
  int REGIONFLAG=0;

  if(vv.userFloat("isBestCand") && vv.userFloat("FullSelTight"))
    set_bit(REGIONFLAG,ZZ);
  if(vv.userFloat("isBestCRZLLos_2P2F") && vv.userFloat("SelZLL_2P2F"))
    set_bit(REGIONFLAG,CRZLLos_2P2F);
  if(vv.userFloat("isBestCRZLLos_3P1F") && vv.userFloat("SelZLL_3P1F"))
    set_bit(REGIONFLAG,CRZLLos_3P1F);

  //For the SR, also fold information about acceptance in CRflag 
  if (isMC_ && test_bit(REGIONFLAG,ZZ)) {
    bool gen_ZZ4lInEtaAcceptance   = false; // All 4 gen leptons in eta acceptance
    bool gen_ZZ4lInEtaPtAcceptance = false; // All 4 gen leptons in eta,pT acceptance
    bool gen_m4l_180               = false; // gen_m4l > 180
    bool gen_ZZInAcceptance        = false; // Unused; old ZZ phase space
    
    mcHistoryTools_->genAcceptance(gen_ZZInAcceptance, gen_ZZ4lInEtaAcceptance, gen_ZZ4lInEtaPtAcceptance, gen_m4l_180);
    
    if (gen_ZZ4lInEtaAcceptance)   set_bit(REGIONFLAG,28);
    if (gen_ZZ4lInEtaPtAcceptance) set_bit(REGIONFLAG,29);
    if (gen_m4l_180)               set_bit(REGIONFLAG,30);
  }
  return REGIONFLAG;
}



#include "FWCore/Framework/interface/MakerMacros.h"
// ---- define this as a plug-in ----------------------------------------
DEFINE_FWK_MODULE(TreePlanter);
