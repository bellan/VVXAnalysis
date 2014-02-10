#include "VVXAnalysis/Producers/plugins/TreePlanter.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/Framework/interface/LuminosityBlock.h"
#include "FWCore/Framework/interface/Run.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/Common/interface/View.h"
#include "DataFormats/Common/interface/MergeableCounter.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"

#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"
#include "SimDataFormats/GeneratorProducts/interface/GenRunInfoProduct.h"

#include "AnalysisDataFormats/CMGTools/interface/PFJet.h"
#include "AnalysisDataFormats/CMGTools/interface/BaseMET.h"
#include "VVXAnalysis/DataFormats/interface/Utilities.h"

#include "ZZAnalysis/AnalysisStep/interface/MCHistoryTools.h"

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
  , theMuonLabel     (config.getParameter<edm::InputTag>("muons"    ))
  , theElectronLabel (config.getParameter<edm::InputTag>("electrons"))
  , theJetLabel      (config.getParameter<edm::InputTag>("jets"     ))
  , theZmmLabel      (config.getParameter<edm::InputTag>("Zmm"      ))
  , theZeeLabel      (config.getParameter<edm::InputTag>("Zee"      ))
  , theWLabel        (config.getParameter<edm::InputTag>("Wjj"      ))
  , theMETLabel      (config.getParameter<edm::InputTag>("MET"      ))
  , theVertexLabel   (config.getParameter<edm::InputTag>("Vertices" ))
  , isMC_            (config.getUntrackedParameter<bool>("isMC",false))
  , sampleType_      (config.getParameter<int>("sampleType"))
  , setup_           (config.getParameter<int>("setup"))
  , theNumberOfEvents(0){
 
  edm::Service<TFileService> fs;
  theTree = fs->make<TTree>("ElderTree","ElderTree");

  if(isMC_){
    thePUInfoLabel        = config.getUntrackedParameter<edm::InputTag>("PUInfo"       , edm::InputTag("addPileupInfo"));
    theGenCategoryLabel   = config.getUntrackedParameter<edm::InputTag>("GenCategory"  , edm::InputTag("genCategory"));
    theGenCollectionLabel = config.getUntrackedParameter<edm::InputTag>("GenCollection", edm::InputTag("genParticlesPruned"));
  }

  initTree();
}


void TreePlanter::beginJob(){
  theTree->Branch("event"     , &event_); 
  theTree->Branch("run"       , &run_); 
  theTree->Branch("lumiBlock" , &lumiBlock_); 

  theTree->Branch("mcprocweight", &mcprocweight_);
  theTree->Branch("puweight"    , &puweight_);
  theTree->Branch("xsec"        , &xsec_);
  theTree->Branch("genCategory" , &genCategory_);

  theTree->Branch("met"   , &met_);
  theTree->Branch("rho"   , &rho_); 
  theTree->Branch("nvtxs" , &nvtx_);

  theTree->Branch("muons"    , &muons_);
  theTree->Branch("electrons", &electrons_);
  theTree->Branch("jets"     , &jets_); 
  theTree->Branch("Zmm"      , &Zmm_); 
  theTree->Branch("Zee"      , &Zee_); 
  theTree->Branch("Wjj"      , &Wjj_);

  theTree->Branch("genParticles", &genParticles_);
}

void TreePlanter::endLuminosityBlock(edm::LuminosityBlock const& lumi, edm::EventSetup const& setup)
{
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
}


void TreePlanter::endRun(const edm::Run& run, const edm::EventSetup& setup){
  edm::Handle<GenRunInfoProduct> genRunInfo;
  run.getByLabel("generator", genRunInfo);
  theXSections.push_back(genRunInfo->crossSection());
}

void TreePlanter::endJob(){

  double sum = 0;
  foreach(const double &xsec, theXSections) sum += xsec;
  
  TVectorD * crossSection =  new TVectorD(1);
  crossSection[0] = sum/theXSections.size();
  TVectorD * genEvents =  new TVectorD(1);
  genEvents[0] = theNumberOfEvents;
  cout << "================================================" << endl;
  cout << "Total number of generated events: " << (*genEvents)[0] 
       <<", cross-section: " << (*crossSection)[0]           << endl;
  cout << "================================================" << endl;
  theTree->GetUserInfo()->Add(genEvents);
  theTree->GetUserInfo()->Add(crossSection);
}



void TreePlanter::initTree(){
  event_     = -1;
  run_       = -1;
  lumiBlock_ = -1;

  mcprocweight_   =  1;
  puweight_       =  1; 
  xsec_           = -1;
  genCategory_    = -1;
  nobservedPUInt_ = -1;
  ntruePUInt_     = -1;

  met_    = phys::Particle();
  nvtx_   = -1;
  rho_    = -1;

  muons_     = std::vector<phys::Lepton>();
  electrons_ = std::vector<phys::Electron>();
  jets_      = std::vector<phys::Jet>();
  Zmm_       = std::vector<phys::Boson<phys::Lepton> >();   
  Zee_       = std::vector<phys::Boson<phys::Electron> >();   
  Wjj_       = std::vector<phys::Boson<phys::Jet> >();   

  genParticles_ = std::vector<phys::Particle>();
 }


void TreePlanter::fillEventInfo(const edm::Event& event){

  // What is missing:
  // rho, weight (MC and PU), xsection, process type

  run_       = event.id().run();
  event_     = event.id().event(); 
  lumiBlock_ = event.luminosityBlock();

  edm::Handle<std::vector<cmg::BaseMET> > met;   event.getByLabel(theMETLabel, met);
  met_ = phys::Particle(phys::Particle::convert(met->front().p4()));

  edm::Handle<std::vector<reco::Vertex> > vertices; event.getByLabel(theVertexLabel, vertices);
  nvtx_ = vertices->size();
    
  if (isMC_) {
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

    // Info about the MC weight
    puweight_ = PUWeighter_.weight(sampleType_, setup_, ntruePUInt_);
    
    MCHistoryTools mch(event);
    mcprocweight_ = mch.gethepMCweight();

  }
}



void TreePlanter::analyze(const edm::Event& event, const edm::EventSetup& setup){
  initTree();
  
  fillEventInfo(event);

  // Load a bunch of objects from the event
  edm::Handle<pat::MuonCollection>       muons        ; event.getByLabel(theMuonLabel    ,     muons);
  edm::Handle<pat::ElectronCollection>   electrons    ; event.getByLabel(theElectronLabel, electrons);
  edm::Handle<std::vector<cmg::PFJet> >  jets         ; event.getByLabel(theJetLabel     ,      jets);
  edm::Handle<edm::View<pat::CompositeCandidate> > Zmm; event.getByLabel(theZmmLabel     ,       Zmm);
  edm::Handle<edm::View<pat::CompositeCandidate> > Zee; event.getByLabel(theZeeLabel     ,       Zee);
  edm::Handle<edm::View<pat::CompositeCandidate> > Wjj; event.getByLabel(theWLabel       ,       Wjj);


  foreach(const pat::Muon& muon, *muons){
    phys::Lepton physmuon = fillLepton(muon);
    muons_.push_back(physmuon);
  }

  foreach(const pat::Electron& electron, *electrons){
    phys::Electron physelectron =  fillElectron(electron);
    electrons_.push_back(physelectron);
  }

  foreach(const cmg::PFJet& jet, *jets){
    phys::Jet physjet = fillJet(jet);
    jets_.push_back(physjet);
  }

  Zmm_ = fillBosons(Zmm, muons_);
 
  Zee_ = fillBosons(Zee, electrons_);

  Wjj_ = fillBosons(Wjj, jets_, 24);



  theTree->Fill();
}


template<typename LEP>
phys::Lepton TreePlanter::fillLepton(const LEP& lepton) const{

  phys::Lepton output(phys::Particle::convert(lepton.p4()),lepton.charge(),13);
  
  output.dxy_             = lepton.userFloat("dxy"            );               
  output.dz_              = lepton.userFloat("dz"             );                
  output.sip_             = lepton.userFloat("SIP"            );
  output.combRelIso_      = lepton.userFloat("combRelIso"     );
  output.pfChargedHadIso_ = lepton.userFloat("PFChargedHadIso");
  output.pfNeutralHadIso_ = lepton.userFloat("PFNeutralHadIso");
  output.pfPhotonIso_     = lepton.userFloat("PFPhotonIso"    );
  output.pfCombRelIso_    = lepton.userFloat("CombRelIsoPF"   );
  output.rho_             = lepton.userFloat("rho"            );
  output.isPF_            = lepton.userInt  ("isPFMuon"       );
  output.matchHLT_        = lepton.userInt  ("HLTMatch"       );
     
  return output; 
}


phys::Electron TreePlanter::fillElectron(const pat::Electron &electron) const{

  phys::Electron output(fillLepton(electron));
  
  output.energy_     = electron.userFloat("energy"    );
  output.phiWidth_   = electron.userFloat("phiWidth"  );
  output.etaWidth_   = electron.userFloat("etaWidth"  );
  output.BDT_        = electron.userFloat("BDT"       );
  output.isBDT_      = electron.userInt  ("isBDT"     );
  output.missingHit_ = electron.userInt  ("missingHit");
  output.nCrystals_  = electron.userInt  ("nCrystals" );

  return output;
}


phys::Jet TreePlanter::fillJet(const cmg::PFJet &jet) const{
  
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
  output.puMVAFull_               = jet.puMva("full");
  output.puMVASimple_             = jet.puMva("simple");
  output.puCutBased_              = jet.puMva("cut-based");
  output.pass_puMVAFull_loose_    = jet.passPuJetId("full"     , PileupJetIdentifier::kLoose);
  output.pass_pUMVAFull_medium_   = jet.passPuJetId("full"     , PileupJetIdentifier::kMedium);
  output.pass_pUMVAFull_tight_    = jet.passPuJetId("full"     , PileupJetIdentifier::kTight); 
  output.pass_puMVASimple_loose_  = jet.passPuJetId("simple"   , PileupJetIdentifier::kLoose); 
  output.pass_puMVASimple_medium_ = jet.passPuJetId("simple"   , PileupJetIdentifier::kMedium); 
  output.pass_puMVASimple_tight_  = jet.passPuJetId("simple"   , PileupJetIdentifier::kTight); 
  output.pass_puCutBased_loose_   = jet.passPuJetId("cut-based", PileupJetIdentifier::kLoose); 
  output.pass_puCutBased_medium_  = jet.passPuJetId("cut-based", PileupJetIdentifier::kMedium);
  output.pass_puCutBased_tight_   = jet.passPuJetId("cut-based", PileupJetIdentifier::kTight);
  
  output.mcPartonFlavour_ = jet.partonFlavour();
 
  return output; 
}

template<typename PAR>
std::vector<phys::Boson<PAR> > TreePlanter::fillBosons(const edm::Handle<edm::View<pat::CompositeCandidate> > & edmBosons, const std::vector<PAR> & physDaughtersCand, int type) const {
  
  std::vector<phys::Boson<PAR> > physBosons;
  
  foreach(const pat::CompositeCandidate& v, *edmBosons){
    
    PAR d0,d1;
    
    foreach(const PAR& particle, physDaughtersCand){
      if( isAlmostEqual(particle.p4().Pt(), v.daughter(0)->pt()) && isAlmostEqual(particle.p4().Eta(), v.daughter(0)->eta())) d0 = particle;
      if( isAlmostEqual(particle.p4().Pt(), v.daughter(1)->pt()) && isAlmostEqual(particle.p4().Eta(), v.daughter(1)->eta())) d1 = particle;
    }
    
    if(d0.id() == 0 || d1.id() == 0) edm::LogError("TreePlanter") << "TreePlanter: VB candidate does not have a matching good particle!";
  
    phys::Boson<PAR> physV(d0, d1, type);
    physBosons.push_back(physV);
  }
  return physBosons;
}



#include "FWCore/Framework/interface/MakerMacros.h"
// ---- define this as a plug-in ----------------------------------------
DEFINE_FWK_MODULE(TreePlanter);
