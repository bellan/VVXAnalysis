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
  , leptonEfficiency_(edm::FileInPath("ZZAnalysis/AnalysisStep/test/Macros/scale_factors_muons2012.root").fullPath(),
		      edm::FileInPath("ZZAnalysis/AnalysisStep/test/Macros/scale_factors_ele2012.root").fullPath())
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
  , theZmmLabel      (config.getParameter<edm::InputTag>("Zmm"      ))
  , theZeeLabel      (config.getParameter<edm::InputTag>("Zee"      ))
  , theWLabel        (config.getParameter<edm::InputTag>("Wjj"      ))
  , theZZ4mLabel     (config.getParameter<edm::InputTag>("ZZ4m"     ))
  , theZZ4eLabel     (config.getParameter<edm::InputTag>("ZZ4e"     ))
  , theZZ2e2mLabel   (config.getParameter<edm::InputTag>("ZZ2e2m"   ))
  , theZllLabel      (config.getParameter<edm::InputTag>("Zll"      ))
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
  , eventsInEtaPtAcceptance_(0){
 
  edm::Service<TFileService> fs;
  theTree = fs->make<TTree>("ElderTree","ElderTree");

  if(isMC_){
    thePUInfoLabel          = config.getUntrackedParameter<edm::InputTag>("PUInfo"         , edm::InputTag("addPileupInfo"));
    theGenCategoryLabel     = config.getUntrackedParameter<edm::InputTag>("GenCategory"    , edm::InputTag("genCategory"));
    theGenCollectionLabel   = config.getUntrackedParameter<edm::InputTag>("GenCollection"  , edm::InputTag("genParticlesPruned"));
    theGenVBCollectionLabel = config.getUntrackedParameter<edm::InputTag>("GenVBCollection", edm::InputTag("genCategory"));
    externalCrossSection_   = config.getUntrackedParameter<double>("XSection",-1);
  }
   
  skimPaths_ = config.getParameter<std::vector<std::string> >("skimPaths");

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
  theTree->Branch("xsec"        , &xsec_);
  theTree->Branch("genCategory" , &genCategory_);

  theTree->Branch("met"   , &met_);
  theTree->Branch("rho"   , &rho_); 
  theTree->Branch("nvtxs" , &nvtx_);

  theTree->Branch("muons"     , &muons_);
  theTree->Branch("electrons" , &electrons_);
  theTree->Branch("jets"      , &jets_); 
  theTree->Branch("ZmmCand"   , &Zmm_); 
  theTree->Branch("ZeeCand"   , &Zee_); 
  theTree->Branch("WjjCand"   , &Wjj_);
  theTree->Branch("ZZ4mCand"  , &ZZ4m_); 
  theTree->Branch("ZZ4eCand"  , &ZZ4e_); 
  theTree->Branch("ZZ2e2mCand", &ZZ2e2m_); 
  theTree->Branch("ZllCand"   , &Zll_); 

  theTree->Branch("genParticles"  , &genParticles_);
  theTree->Branch("genVBParticles", &genVBParticles_);

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
}


void TreePlanter::endRun(const edm::Run& run, const edm::EventSetup& setup){
  if(isMC_){
    edm::Handle<GenRunInfoProduct> genRunInfo;
    run.getByLabel("generator", genRunInfo);
    theXSections.push_back(genRunInfo->crossSection());
  }
}

void TreePlanter::endJob(){

  if(mcHistoryTools_) delete mcHistoryTools_;

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
    
    edm::Service<TFileService> fs;
    TTree *countTree = fs->make<TTree>("HollyTree","HollyTree");
    countTree->Branch("signalDefinition"      , &signalDefinition_);
    countTree->Branch("genEvents"             , &theNumberOfEvents);
    countTree->Branch("analyzedEvents"        , &theNumberOfAnalyzedEvents);
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


    countTree->Fill();
  }
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

  xsec_           = -1.;
  genCategory_    = -1;
  nobservedPUInt_ = -1;
  ntruePUInt_     = -1;

  met_    = phys::Particle();
  nvtx_   = -1;
  rho_    = -1;

  muons_     = std::vector<phys::Lepton>();
  electrons_ = std::vector<phys::Electron>();
  jets_      = std::vector<phys::Jet>();
  Zmm_       = std::vector<phys::Boson<phys::Lepton>   >();   
  Zee_       = std::vector<phys::Boson<phys::Electron> >();   
  Wjj_       = std::vector<phys::Boson<phys::Jet>      >();   
  ZZ4m_      = std::vector<phys::DiBoson<phys::Lepton  , phys::Lepton>   >();
  ZZ4e_      = std::vector<phys::DiBoson<phys::Electron, phys::Electron> >();
  ZZ2e2m_    = std::vector<phys::DiBoson<phys::Electron, phys::Lepton>   >();
  Zll_       = std::vector<phys::DiBoson<phys::Lepton  , phys::Lepton>   >();

  genParticles_ = std::vector<phys::Particle>();
  genVBParticles_ = std::vector<phys::Boson<phys::Particle> >();
 }


bool TreePlanter::fillEventInfo(const edm::Event& event){

  // Fill some info abut acceptance before cutting away events. Beware: if the signal is defined a-posteriori, we will have a problem. For that case, we need to
  // explicitly check here that we are counting signal and not irreducible background.
  if (isMC_) {
    if(mcHistoryTools_) delete mcHistoryTools_;
    mcHistoryTools_ = new MCHistoryTools(event, sampleName_);
    bool gen_ZZ4lInEtaAcceptance   = false; // All 4 gen leptons in eta acceptance
    bool gen_ZZ4lInEtaPtAcceptance = false; // All 4 gen leptons in eta,pT acceptance
    bool gen_m4l_180               = false; // gen_m4l > 180
    bool gen_ZZInAcceptance        = false; // Unused; old ZZ phase space
    mcHistoryTools_->genAcceptance(gen_ZZInAcceptance, gen_ZZ4lInEtaAcceptance, gen_ZZ4lInEtaPtAcceptance, gen_m4l_180);
    if (gen_ZZ4lInEtaAcceptance)   ++eventsInEtaAcceptance_;  
    if (gen_ZZ4lInEtaPtAcceptance) ++eventsInEtaPtAcceptance_;
  }
    
  // Check trigger request. Actually, it is a very very loose request, not the actual one, that instead should be
  // asked to the specific final state
  passTrigger_ = filterController_.passTrigger(NONE, event, triggerWord_);
  if (applyTrigger_ && !passTrigger_) return false;
  
  // Check Skim requests
  passSkim_ = filterController_.passSkim(event, triggerWord_);
  if (applySkim_    && !passSkim_)   return false;
  
  // Apply MC filter
  //if (isMC_ && !(filterController_.passMCFilter(event))) return false;

  run_       = event.id().run();
  event_     = event.id().event(); 
  lumiBlock_ = event.luminosityBlock();

  edm::Handle<std::vector<cmg::BaseMET> > met;   event.getByLabel(theMETLabel, met);
  met_ = phys::Particle(phys::Particle::convert(met->front().p4()));

  edm::Handle<std::vector<reco::Vertex> > vertices; event.getByLabel(theVertexLabel, vertices);
  nvtx_ = vertices->size();
    
  if(isMC_){
    edm::Handle<std::vector<PileupSummaryInfo> > puInfo; event.getByLabel(thePUInfoLabel, puInfo);
    foreach(const PileupSummaryInfo& pui, *puInfo)
      if(pui.getBunchCrossing() == 0) { 
	nobservedPUInt_  = pui.getPU_NumInteractions();
	ntruePUInt_      = pui.getTrueNumInteractions();
	break;
      }
    
    edm::Handle<edm::View<reco::Candidate> > genParticles;
    event.getByLabel(theGenCollectionLabel,  genParticles);
    
    for (edm::View<reco::Candidate>::const_iterator p = genParticles->begin(); p != genParticles->end(); ++p) 
      if (p->status() ==3 && std::distance(genParticles->begin(),p) > 5)
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
    
    // Info about the MC weight
    puweight_ = PUWeighter_.weight(sampleType_, setup_, ntruePUInt_);
    
    mcprocweight_ = mcHistoryTools_->gethepMCweight();

    // Sum of weight, particularly imprtant for MCs that return also negative weights
    // or, in general, weighted events
    sumpuweights_       += puweight_;
    summcprocweights_   += mcprocweight_;
    sumpumcprocweights_ += puweight_*mcprocweight_;
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
  edm::Handle<edm::View<pat::CompositeCandidate> > Zmm    ; event.getByLabel(theZmmLabel     ,       Zmm);
  edm::Handle<edm::View<pat::CompositeCandidate> > Zee    ; event.getByLabel(theZeeLabel     ,       Zee);
  edm::Handle<edm::View<pat::CompositeCandidate> > Wjj    ; event.getByLabel(theWLabel       ,       Wjj);
  edm::Handle<edm::View<pat::CompositeCandidate> > ZZ4m   ; event.getByLabel(theZZ4mLabel    ,      ZZ4m);
  edm::Handle<edm::View<pat::CompositeCandidate> > ZZ4e   ; event.getByLabel(theZZ4eLabel    ,      ZZ4e);
  edm::Handle<edm::View<pat::CompositeCandidate> > ZZ2e2m ; event.getByLabel(theZZ2e2mLabel  ,    ZZ2e2m);
  // Collections for CR
  edm::Handle<edm::View<pat::CompositeCandidate> > Zll    ; event.getByLabel(theZllLabel     ,       Zll);
  

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
  Zmm_    = fillBosons<pat::Muon,phys::Lepton>(Zmm, 23);

  Zee_    = fillBosons<pat::Electron,phys::Electron>(Zee, 23);

  Wjj_    = fillBosons<cmg::PFJet,phys::Jet>(Wjj, 24);


  // The bosons have NOT any requirement on the quality of their daughters, only the flag is set (because of the same code is usd for CR too)

  ZZ4m_   = fillDiBosons<pat::Muon,phys::Lepton,pat::Muon,phys::Lepton>(MMMM,ZZ4m);

  ZZ4e_   = fillDiBosons<pat::Electron,phys::Electron,pat::Electron,phys::Electron>(EEEE,ZZ4e);

  ZZ2e2m_ = fillDiBosons<pat::Electron,phys::Electron,pat::Muon,phys::Lepton>(EEMM,ZZ2e2m);

  Zll_    = fillDiBosons<phys::Lepton,phys::Lepton>(ZLL,Zll);


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
  output.pfCombRelIso_    = lepton.userFloat("CombRelIsoPF"     );
  output.rho_             = lepton.userFloat("rho"              );
  output.isPF_            = lepton.userFloat("isPFMuon"         );
  output.matchHLT_        = lepton.userFloat("HLTMatch"         );
  output.isGood_          = lepton.userFloat("isGood"           );
  output.efficiencySF_    = leptonEfficiency_.scaleFactor(output);
  
     
  return output; 
}


phys::Electron TreePlanter::fill(const pat::Electron &electron) const{

  phys::Electron output(fillLepton(electron));
  
  output.energy_     = electron.userFloat("energy"    );
  output.phiWidth_   = electron.userFloat("phiWidth"  );
  output.etaWidth_   = electron.userFloat("etaWidth"  );
  output.BDT_        = electron.userFloat("BDT"       );
  output.isBDT_      = electron.userFloat("isBDT"     );
  output.missingHit_ = electron.userInt  ("missingHit");
  output.nCrystals_  = electron.userInt  ("nCrystals" );

  return output;
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
  if(v.hasUserFloat("dauWithFSR") && v.userFloat("dauWithFSR") >= 0){
    phys::Particle photon(phys::Particle::convert(v.daughter(2)->p4()), 0, 22);
    physV.addFSR(v.userFloat("dauWithFSR"), photon);
  }
  
  // Add quality of daughters, right now only for ll pairs
  if(v.hasUserFloat("GoodLeptons")) physV.hasGoodDaughters_ = v.userFloat("GoodLeptons");
  

  return physV;
}


// Right now, only for ZZ
template<typename T1, typename PAR1, typename T2, typename PAR2>
phys::DiBoson<PAR1,PAR2> TreePlanter::fillDiBoson(Channel channel, const pat::CompositeCandidate& edmVV) const{
  
  const pat::CompositeCandidate* edmV0;
  const pat::CompositeCandidate* edmV1;
  
  if (channel != ZL) { // Regular 4l candidates
    edmV0   = dynamic_cast<const pat::CompositeCandidate*>(edmVV.daughter("Z1")->masterClone().get());      
    edmV1   = dynamic_cast<const pat::CompositeCandidate*>(edmVV.daughter("Z2")->masterClone().get());
    
  } else {              // Special handling of Z+l candidates 
    edmV0   = dynamic_cast<const pat::CompositeCandidate*>(edmVV.daughter(0)->masterClone().get());      
    edmV1   = dynamic_cast<const pat::CompositeCandidate*>(edmVV.daughter(1)->masterClone().get());
  }
  
  phys::Boson<PAR1> V0;
  phys::Boson<PAR2> V1;
  
  // The first boson is always a good Z, also in the CR. For the other particle assign 23 if it is a true Z from SR
  // or 26 if the two additional leptons comes from LL, 27 if the additional lepton is only one.
  int idV1 = (channel != ZL && channel != ZLL) ? 23 : (channel == ZLL ? 26 : (channel == ZL ? 27 : 0));

  if(dynamic_cast<const T1*>(edmV0->daughter(0)->masterClone().get()) && dynamic_cast<const T2*>(edmV1->daughter(0)->masterClone().get())){
    V0 = fillBoson<T1,PAR1>(*edmV0, 23, false);
    V1 = fillBoson<T2,PAR2>(*edmV1, idV1, false);
  }
  else if(dynamic_cast<const T2*>(edmV0->daughter(0)->masterClone().get()) && dynamic_cast<const T1*>(edmV1->daughter(0)->masterClone().get())){
    V0 = fillBoson<T1,PAR1>(*edmV1, idV1, false);
    V1 = fillBoson<T2,PAR2>(*edmV0, 23, false);
  }
  else{
    edm::LogError("TreePlanter") << "Do not know what to cast in fillDiBosons" << endl;
    return phys::DiBoson<PAR1,PAR2>();
  }
  
  phys::DiBoson<PAR1,PAR2> VV(V0, V1);
  VV.isBestCand_  = edmVV.userFloat("isBestCand");
  VV.passFullSel_ = edmVV.userFloat("FullSel");
  VV.regionWord_  = computeCRFlag(channel,edmVV);
  VV.triggerWord_ = triggerWord_;
  VV.passTrigger_ = filterController_.passTrigger(channel, triggerWord_); // triggerWord_ needs to be filled beforehand (as it is).
  
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
template<typename T1, typename PAR1, typename T2, typename PAR2>
std::vector<phys::DiBoson<PAR1,PAR2> > TreePlanter::fillDiBosons(Channel channel, const edm::Handle<edm::View<pat::CompositeCandidate> > & edmDiBosons) const{

  std::vector<phys::DiBoson<PAR1,PAR2> > physDiBosons;
  
  foreach(const pat::CompositeCandidate& edmVV, *edmDiBosons){
    phys::DiBoson<PAR1,PAR2> physVV = fillDiBoson<T1,PAR1,T2,PAR2>(channel, edmVV);
    if(physVV.isValid()) physDiBosons.push_back(physVV);
  }

  return physDiBosons;
}



// ZLL
template<typename PAR1, typename PAR2>
std::vector<phys::DiBoson<PAR1,PAR2> > TreePlanter::fillDiBosons(Channel channel, const edm::Handle<edm::View<pat::CompositeCandidate> > & edmDiBosons) const{

  std::vector<phys::DiBoson<PAR1,PAR2> > physDiBosons;
  
  foreach(const pat::CompositeCandidate& edmVV, *edmDiBosons){
    
    phys::DiBoson<PAR1,PAR2> physVV;

    int count = 0;
    
    int regionWord = computeCRFlag(channel,edmVV);
    // 4 mu cases
    if(test_bit(regionWord,CRMMMMss ) || test_bit(regionWord,CRMMMMos )){
      physVV = fillDiBoson<pat::Muon,phys::Lepton,pat::Muon,phys::Lepton>(channel, edmVV);
      ++count;
    }  
    
    // 4 e cases
    if(test_bit(regionWord, CREEEEss) || test_bit(regionWord, CREEEEos)){
      physVV = fillDiBoson<pat::Electron,phys::Lepton,pat::Electron,phys::Lepton>(channel, edmVV);      
      ++count;
    }  
    
    // 2 e 2 mu cases
    if(test_bit(regionWord, CREEMMss) || test_bit(regionWord, CREEMMos)){
      physVV = fillDiBoson<pat::Electron,phys::Lepton,pat::Muon,phys::Lepton>(channel, edmVV);
      ++count;
    }  
    
    // 2 mu 2 e cases
    if(test_bit(regionWord, CRMMEEss) || test_bit(regionWord, CRMMEEos)){
      physVV = fillDiBoson<pat::Muon,phys::Lepton,pat::Electron,phys::Lepton>(channel, edmVV);
      ++count;
    }  
    
    if(physVV.isValid()) physDiBosons.push_back(physVV);    
  }
  
  return physDiBosons;
}







int TreePlanter::computeCRFlag(Channel channel, const pat::CompositeCandidate & vv) const{
  int CRFLAG=0;
    
  if (channel==ZLL) {
    if(vv.userFloat("isBestCRZLL")&&vv.userFloat("CRZLL"))
      set_bit(CRFLAG,CRZLL);
    if(vv.userFloat("isBestCRMMMMss")&&vv.userFloat("CRLLLL"))
      set_bit(CRFLAG,CRMMMMss);
    if(vv.userFloat("isBestCRMMMMos")&&vv.userFloat("CRLLLL"))
      set_bit(CRFLAG,CRMMMMos);
    if(vv.userFloat("isBestCREEEEss")&&vv.userFloat("CRLLLL"))
      set_bit(CRFLAG,CREEEEss);
    if(vv.userFloat("isBestCREEEEos")&&vv.userFloat("CRLLLL"))
      set_bit(CRFLAG,CREEEEos);
    if(vv.userFloat("isBestCREEMMss")&&vv.userFloat("CRLLLL"))
      set_bit(CRFLAG,CREEMMss);
    if(vv.userFloat("isBestCREEMMos")&&vv.userFloat("CRLLLL"))
      set_bit(CRFLAG,CREEMMos);
    if(vv.userFloat("isBestCRMMEEss")&&vv.userFloat("CRLLLL"))
      set_bit(CRFLAG,CRMMEEss);
    if(vv.userFloat("isBestCRMMEEos")&&vv.userFloat("CRLLLL"))
      set_bit(CRFLAG,CRMMEEos);
    if(vv.userFloat("isBestCRZLL")&&vv.userFloat("CRZLLHiSIP")) 
      set_bit(CRFLAG,CRZLLHiSIP);
    if(vv.userFloat("isBestCRZMM")&&vv.userFloat("CRZLLHiSIP"))
      set_bit(CRFLAG,CRZLLHiSIPMM);
    if(vv.userFloat("isBestCRZLLHiSIPKin")&&vv.userFloat("CRZLLHiSIPKin"))
      set_bit(CRFLAG,CRZLLHiSIPKin);  
  }

  //For the SR, also fold information about acceptance in CRflag 
  if (isMC_ && (channel==EEEE||channel==MMMM||channel==EEMM)) {
    bool gen_ZZ4lInEtaAcceptance   = false; // All 4 gen leptons in eta acceptance
    bool gen_ZZ4lInEtaPtAcceptance = false; // All 4 gen leptons in eta,pT acceptance
    bool gen_m4l_180               = false; // gen_m4l > 180
    bool gen_ZZInAcceptance        = false; // Unused; old ZZ phase space

    mcHistoryTools_->genAcceptance(gen_ZZInAcceptance, gen_ZZ4lInEtaAcceptance, gen_ZZ4lInEtaPtAcceptance, gen_m4l_180);

    if (gen_ZZ4lInEtaAcceptance)   set_bit(CRFLAG,28);
    if (gen_ZZ4lInEtaPtAcceptance) set_bit(CRFLAG,29);
    if (gen_m4l_180)               set_bit(CRFLAG,30);
  }
  return CRFLAG;
}



#include "FWCore/Framework/interface/MakerMacros.h"
// ---- define this as a plug-in ----------------------------------------
DEFINE_FWK_MODULE(TreePlanter);
