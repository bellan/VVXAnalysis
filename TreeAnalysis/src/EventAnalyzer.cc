/** \class EventAnalyzer
 *  Base class for event analyzers. Analyzers have to inherit from this class 
 *  and implement the pure virtual function analyze(), called each event.
 *
 *  $Date: 2013/03/15 13:37:31 $
 *  $Revision: 1.4 $
 *  \author R. Bellan - UNITO <riccardo.bellan@cern.ch>
 */

#include "VVXAnalysis/TreeAnalysis/interface/EventAnalyzer.h"
#include "VVXAnalysis/Commons/interface/Colours.h"
#include "VVXAnalysis/Commons/interface/Comparators.h"
#include "VVXAnalysis/Commons/interface/Constants.h"

#include <TChain.h>
#include <TFile.h>
#include <TH1F.h>
#include <TH1I.h>
#include <TString.h>

#include <boost/lexical_cast.hpp>
#include <boost/foreach.hpp>
#define foreach BOOST_FOREACH

using std::cout;
using std::endl;
using std::flush;
using namespace colour;


// ------------------------------------------------------------------------------------------ //
// ------------------------------- C'tor/Des'tor/Init the tree ------------------------------ //
// ------------------------------------------------------------------------------------------ //



EventAnalyzer::EventAnalyzer(SelectorBase& aSelector,
			     std::string filename, 
			     double lumi, 
			     double externalXSection, bool doBasicPlots)
  : select(aSelector)
  , doBasicPlots_(doBasicPlots)
  , theMCInfo(filename, lumi, externalXSection)
  , theWeight(1.)
  , theCutCounter(0.)
  , theInputWeightedEvents(0.){

  TChain *tree = new TChain("treePlanter/ElderTree");
  tree->Add(filename.c_str());

  if (tree == 0) std::cout<<Important("Error in EventAnalyzer ctor:")<<" The tree has a null pointer."<<std::endl;

  Init(tree);
}



EventAnalyzer::~EventAnalyzer(){
  if (!theTree) return;
  delete theTree->GetCurrentFile();
}



void EventAnalyzer::Init(TTree *tree)
{
  TH1F::SetDefaultSumw2(kTRUE);
  
  // Set branch addresses and branch pointers
  if (!tree) return;
  theTree = tree;
  fCurrent = -1;
  
  // Muons   
  muons     = 0; b_muons     = 0; theTree->SetBranchAddress("muons"    , &muons    , &b_muons    );
  
  // Electrons   
  electrons = 0; b_electrons = 0; theTree->SetBranchAddress("electrons", &electrons, &b_electrons);
  
  // Jets   
  pjets = 0;      b_pjets = 0;    theTree->SetBranchAddress("jets", &pjets, &b_pjets);
  jets  = new std::vector<phys::Jet>(); centralJets  = new std::vector<phys::Jet>();

  // Bosons   
  Zmm = new std::vector<phys::Boson<phys::Lepton> >()  ; ZmmCand = 0; b_ZmmCand = 0; theTree->SetBranchAddress("ZmmCand", &ZmmCand, &b_ZmmCand);
  Zee = new std::vector<phys::Boson<phys::Electron> >(); ZeeCand = 0; b_ZeeCand = 0; theTree->SetBranchAddress("ZeeCand", &ZeeCand, &b_ZeeCand);
  Wjj = new std::vector<phys::Boson<phys::Jet> > ()    ; WjjCand = 0; b_WjjCand = 0; theTree->SetBranchAddress("WjjCand", &WjjCand, &b_WjjCand);

  // DiBosons
  ZZ4m   = 0; b_ZZ4m   = 0; theTree->SetBranchAddress("ZZ4mCand"  , &ZZ4m  , &b_ZZ4m  );
  ZZ4e   = 0; b_ZZ4e   = 0; theTree->SetBranchAddress("ZZ4eCand"  , &ZZ4e  , &b_ZZ4e  );
  ZZ2e2m = 0; b_ZZ2e2m = 0; theTree->SetBranchAddress("ZZ2e2mCand", &ZZ2e2m, &b_ZZ2e2m);

  ZZ = 0;


  // Gen Particles   
  genParticles   = 0;                                                b_genParticles   = 0; theTree->SetBranchAddress("genParticles"  , &genParticles  , &b_genParticles);
  genVBParticles = new std::vector<phys::Boson<phys::Particle> > (); b_genVBParticles = 0; theTree->SetBranchAddress("genVBParticles", &genVBParticles, &b_genVBParticles);

  // MET
  met = new phys::Particle();
  b_met         = 0; theTree->SetBranchAddress("met"    , &met    ,  &b_met    );
  
  // Other events variables
  b_event       = 0; theTree->SetBranchAddress("event"    , &event    , &b_event    );
  b_run         = 0; theTree->SetBranchAddress("run"      , &run      , &b_run      );
  b_lumiBlock   = 0; theTree->SetBranchAddress("lumiBlock", &lumiBlock, &b_lumiBlock);

  b_nvtx        = 0; theTree->SetBranchAddress("nvtxs"  , &nvtx   , &b_nvtx   );
  b_rho         = 0; theTree->SetBranchAddress("rho"    , &rho    , &b_rho    );
  
  // MC related variables
  b_puweight     = 0; theTree->SetBranchAddress("puweight"    , &theMCInfo.puweight_     , &b_puweight );
  b_mcprocweight = 0; theTree->SetBranchAddress("mcprocweight", &theMCInfo.mcprocweight_ , &b_mcprocweight);
  b_genCategory  = 0; theTree->SetBranchAddress("genCategory" , &genCategory             , &b_genCategory  );

  // Info about selections
  b_passTrigger = 0; theTree->SetBranchAddress("passTrigger", &passTrigger, &b_passTrigger); 
  b_passSkim    = 0; theTree->SetBranchAddress("passSkim"   , &passSkim   , &b_passSkim   ); 
  b_triggerWord = 0; theTree->SetBranchAddress("triggerWord", &triggerWord, &b_triggerWord); 
}



// ------------------------------------------------------------------------------------------ //
// ------------------------------------- Master the loop ------------------------------------ //
// ------------------------------------------------------------------------------------------ //



Int_t EventAnalyzer::GetEntry(Long64_t entry){
  // Read contents of entry.
      
  if (!theTree) return 0;
  
  int e =  theTree->GetEntry(entry);
  
  stable_sort(muons->begin(),     muons->end(),     phys::PtComparator());
  stable_sort(electrons->begin(), electrons->end(), phys::PtComparator());
  stable_sort(pjets->begin(),     pjets->end(),     phys::PtComparator());

  // Some selection on jets
  jets->clear(); centralJets->clear();
  foreach(const phys::Jet &jet, *pjets)
    if(jet.pt() > 30){
      if(fabs(jet.eta()) < 4.7) jets->push_back(jet);
      if(fabs(jet.eta()) < 2.5) centralJets->push_back(jet);
    }

  

  Zmm->clear(); Zee->clear(); Wjj->clear();

  foreach(const phys::Boson<phys::Lepton> z, *ZmmCand)
    if(select(z)) Zmm->push_back(z);
  foreach(const phys::Boson<phys::Electron> z, *ZeeCand)
    if(select(z)) Zee->push_back(z);
  foreach(const phys::Boson<phys::Jet> w, *WjjCand)
    if(select(w)) Wjj->push_back(w);

  stable_sort(Zmm->begin(), Zmm->end(), phys::PtComparator());
  stable_sort(Zee->begin(), Zee->end(), phys::PtComparator());
  stable_sort(Wjj->begin(), Wjj->end(), phys::PtComparator());
  
  int totCand = ZZ4m->size() + ZZ4e->size() + ZZ2e2m->size();
  theHistograms.fill("nZZCandidates",     "Number of good candidates in the event", 10, 0, 10, totCand       , 1);
  theHistograms.fill("nZZ4eCandidates",   "Number of good candidates in the event", 10, 0, 10, ZZ4e->size()  , 1);
  theHistograms.fill("nZZ4mCandidates",   "Number of good candidates in the event", 10, 0, 10, ZZ4m->size()  , 1);
  theHistograms.fill("nZZ2e2mCandidates", "Number of good candidates in the event", 10, 0, 10, ZZ2e2m->size(), 1);
  
  int triggers = 0;
  if(totCand > 0){
    
    if(ZZ4m->size() == 1 && ZZ4m->front().passTrigger()){ ++triggers;
      ZZ = new phys::DiBoson<phys::Lepton  , phys::Lepton>(ZZ4m->front().clone<phys::Lepton,phys::Lepton>());
    }
    if(ZZ4e->size() == 1 && ZZ4e->front().passTrigger()){ ++triggers;
      ZZ = new phys::DiBoson<phys::Lepton  , phys::Lepton>(ZZ4e->front().clone<phys::Lepton,phys::Lepton>());
    }
    if(ZZ2e2m->size() == 1 && ZZ2e2m->front().passTrigger()){ ++triggers;
      ZZ = new phys::DiBoson<phys::Lepton  , phys::Lepton>(ZZ2e2m->front().clone<phys::Lepton,phys::Lepton>());
    }
  }

  theHistograms.fill("GoodTriggerableCands", "Number of good triggerable candidates in the event", 10, 0, 10, triggers, 1);
  if(triggers != 1) return 0;
  
  theWeight = theMCInfo.weight(*ZZ);
  
  theInputWeightedEvents += theWeight;

  return e;
}



Long64_t EventAnalyzer::LoadTree(Long64_t entry){

  if (!theTree) return -5;
  Long64_t centry = theTree->LoadTree(entry);

  if (centry < 0) return centry;
  if (!theTree->InheritsFrom(TChain::Class()))  return centry;

  TChain *chain = (TChain*)theTree;
  if (chain->GetTreeNumber() != fCurrent) fCurrent = chain->GetTreeNumber();
  
  return centry;
}



void EventAnalyzer::loop(const std::string outputfile){

  if (theTree == 0) return;

  Long64_t nentries = theTree->GetEntries();  

  begin();

  for (Long64_t jentry=0; jentry<nentries; ++jentry) {

    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    if (!GetEntry(jentry)) continue;

    if (cut() < 0) continue;
    theCutCounter += theWeight;

    theHistograms.fill("weight",100, 0, 200, theWeight);

    if(doBasicPlots_) fillBasicPlots();
    analyze();
  }

  TFile fout(TString(outputfile),"RECREATE");
  fout.cd(); 
  theHistograms.write(fout);

  end(fout);
  fout.Close();
  //cout<<"Events in input: " << Green(theInputWeightedEvents)<< endl;
  cout<<"Events passing all cuts: "<< Green(theCutCounter) << endl;
}



// ------------------------------------------------------------------------------------------ //
// -------------------------------------- Preselection -------------------------------------- //
// ------------------------------------------------------------------------------------------ //


Int_t EventAnalyzer::cut() {
  
  bool pass = true;
  
  return pass ? 1 : -1;
}


// ------------------------------------------------------------------------------------------ //
// --------------------------------------- Histograms --------------------------------------- //
// ------------------------------------------------------------------------------------------ //



void EventAnalyzer::fillBasicPlots(){
  theHistograms.fill<TH1I>("nvtx"     , "Number of vertices" , 100, 0, 100, nvtx             , theWeight);
  theHistograms.fill      ("rho"      , "Mean energy density", 100, 0, 50 , rho              , theWeight);
  theHistograms.fill      ("met"      , "Missing energy"     , 200, 0, 800, met->pt()        , theWeight);
  theHistograms.fill<TH1I>("nmuons"    ,"Number of muons"    ,  10, 0, 10 , muons->size()    , theWeight);
  theHistograms.fill<TH1I>("nelectrons","Number of electrons",  10, 0, 10 , electrons->size(), theWeight);

  foreach(const phys::Lepton& mu , *muons)     fillLeptonPlots  ("mu",  mu  );
  foreach(const phys::Electron& e, *electrons) fillElectronPlots("e" ,  e   );
  foreach(const phys::Jet& jet   , *jets)      fillJetPlots     ("jet", jet );
}



void EventAnalyzer::fillParticlePlots(const std::string &type, const phys::Particle & lepton){
  theHistograms.fill(type+"_pt" ,    "p_{T} spectrum", 100,   0   , 500   ,lepton.pt()    , theWeight);
  theHistograms.fill(type+"_eta",    "#eta spectrum" , 100,  -2.5 ,   2.5 ,lepton.eta()   , theWeight);
  theHistograms.fill(type+"_phi",    "#phi spectrum" ,  50,  -3.15,   3.15,lepton.phi()   , theWeight);
  theHistograms.fill(type+"_charge", "charge"        ,  50,  -25  ,  25   ,lepton.charge(), theWeight);
}



void EventAnalyzer::fillLeptonPlots(const std::string &type, const phys::Lepton & lepton){
  
  fillParticlePlots(type, lepton);

  theHistograms.fill(type+"_dxy"             , "d_{xy}"         , 200,   0,   0.5, lepton.dxy()            , theWeight);   
  theHistograms.fill(type+"_dz"              , "d_{z}"          , 200,   0,   0.5, lepton.dz()             , theWeight);
  theHistograms.fill(type+"_sip"             , "sip"            , 150,   0,  15  , lepton.sip()            , theWeight); 
  theHistograms.fill(type+"_combRelIso"      , "combRelIso"     , 200,   0,   1  , lepton.combRelIso()     , theWeight);       
  theHistograms.fill(type+"_pfCombRelIso"    , "pfCombRelIso"   , 200,   0,   1  , lepton.pfChargedHadIso(), theWeight);           
  theHistograms.fill(type+"_pfNeutralHadIso" , "pfNeutralHadIso", 200,   0,   1  , lepton.pfNeutralHadIso(), theWeight);        
  theHistograms.fill(type+"_pfChargedHadIso" , "pfChargedHadIso", 200,   0,   1  , lepton.pfPhotonIso()    , theWeight);        
  theHistograms.fill(type+"_pfPhotonHadIso"  , "pfPhotonHadIso" , 200,   0,   1  , lepton.pfCombRelIso()   , theWeight);         
  theHistograms.fill(type+"_rho"             , "rho"            ,  50,   0,  50  , lepton.rho()            , theWeight);                    
}



void EventAnalyzer::fillElectronPlots (const std::string &type, const phys::Electron &electron){
  fillLeptonPlots           (type, electron);
  fillExtraPlotsForElectrons(type, electron);
}



void EventAnalyzer::fillExtraPlotsForElectrons(const std::string &type, const phys::Electron &electron){

  theHistograms.fill      (type+"_energy"    , "energy spectrum"   , 100,  0   , 500   , electron.energy()    , theWeight);  
  theHistograms.fill      (type+"_phiWidth"  , "#phi width"        ,  50,  0   ,   2   , electron.phiWidth()  , theWeight);  
  theHistograms.fill      (type+"_etaWidth"  , "#eta width"        ,  50,  0   ,   2   , electron.etaWidth()  , theWeight);  
  theHistograms.fill      (type+"_BDT"       , "BDT"               , 100, -1   ,   1   , electron.BDT()       , theWeight);  
  theHistograms.fill<TH1I>(type+"_isBDT"     , "is BDT?"           ,   2,  0   ,   2   , electron.isBDT()     , theWeight);  // void??
  theHistograms.fill<TH1I>(type+"_missingHit", "Missing hits"      ,  50,  0   ,  50   , electron.missingHit(), theWeight);  // void?? 
  theHistograms.fill<TH1I>(type+"_nCrystals" , "Number of Crystals", 200,  0   , 200   , electron.nCrystals() , theWeight);  
}



void EventAnalyzer::fillJetPlots(const std::string &type, const phys::Jet      &jet){

  fillParticlePlots(type, jet);


  theHistograms.fill<TH1I>(type+"_nConstituents", 120, 0, 120, jet.nConstituents(),theWeight);
  theHistograms.fill<TH1I>(type+"_nCharged"     , 100, 0, 100, jet.nCharged(),theWeight);
  theHistograms.fill<TH1I>(type+"_nNeutral"     , 100, 0, 100, jet.nNeutral(),theWeight);

  theHistograms.fill(type+"_neutralHadronEnergyFraction", 100, 0, 1, jet.neutralHadronEnergyFraction(), theWeight);
  theHistograms.fill(type+"_chargedHadronEnergyFraction", 100, 0, 1, jet.chargedHadronEnergyFraction(), theWeight);
  theHistograms.fill(type+"_chargedEmEnergyFraction"    , 100, 0, 1, jet.chargedEmEnergyFraction    (), theWeight);
  theHistograms.fill(type+"_neutralEmEnergyFraction"    , 100, 0, 1, jet.neutralEmEnergyFraction    (), theWeight);
  theHistograms.fill(type+"_muonEnergyFraction"         , 100, 0, 1, jet.muonEnergyFraction         (), theWeight);

  theHistograms.fill(type+"_csvtagger"    , 200,  0, 1   , jet.csvtagger    (), theWeight);
  theHistograms.fill(type+"_girth"        , 200,  0, 1   , jet.girth        (), theWeight);
  theHistograms.fill(type+"_girth_charged", 200,  0, 1   , jet.girth_charged(), theWeight);
  theHistograms.fill(type+"_ptd"          , 100,  0, 1   , jet.ptd          (), theWeight);
  theHistograms.fill(type+"_rms"          ,  50,  0, 0.15, jet.rms          (), theWeight);
  theHistograms.fill(type+"_beta"         ,  50,  0, 1   , jet.beta         (), theWeight);
  theHistograms.fill(type+"_jetArea"      ,  60,  0, 1.2 , jet.jetArea      (), theWeight);
  theHistograms.fill(type+"_secvtxMass"   ,  50,  0, 5   , jet.secvtxMass   (), theWeight);
  theHistograms.fill(type+"_Lxy"          ,  50,  0, 10  , jet.Lxy          (), theWeight);
  theHistograms.fill(type+"_LxyErr"       ,  50,  0, 1.5 , jet.LxyErr       (), theWeight);
  theHistograms.fill(type+"_rawFactor"    ,  50,  0, 2.5 , jet.rawFactor    (), theWeight);

  theHistograms.fill(type+"_uncOnFourVectorScale", 50, 0, 0.1, jet.uncOnFourVectorScale(), theWeight);
  theHistograms.fill(type+"_puMVAFull"  , 100, -1, 1, jet.puMVAFull  (), theWeight);
  theHistograms.fill(type+"_puMVASimple", 100, -1, 1, jet.puMVASimple(), theWeight);
  theHistograms.fill(type+"_puCutBased" , 100, -1, 1, jet.puCutBased (), theWeight); // void??

  theHistograms.fill<TH1I>(type+"_pass_puMVAFull_loose"   , 2, 0, 2, jet.pass_puMVAFull_loose   (), theWeight);
  theHistograms.fill<TH1I>(type+"_pass_pUMVAFull_medium"  , 2, 0, 2, jet.pass_pUMVAFull_medium  (), theWeight);
  theHistograms.fill<TH1I>(type+"_pass_pUMVAFull_tight"   , 2, 0, 2, jet.pass_pUMVAFull_tight   (), theWeight);
				                          	                             
  theHistograms.fill<TH1I>(type+"_pass_puMVASimple_loose" , 2, 0, 2, jet.pass_puMVASimple_loose (), theWeight);
  theHistograms.fill<TH1I>(type+"_pass_puMVASimple_medium", 2, 0, 2, jet.pass_puMVASimple_medium(), theWeight);
  theHistograms.fill<TH1I>(type+"_pass_puMVASimple_tight" , 2, 0, 2, jet.pass_puMVASimple_tight (), theWeight);
				                          	                             
  theHistograms.fill<TH1I>(type+"_pass_puCutBased_loose"  , 2, 0, 2, jet.pass_puCutBased_loose  (), theWeight); // void??
  theHistograms.fill<TH1I>(type+"_pass_puCutBased_medium" , 2, 0, 2, jet.pass_puCutBased_medium (), theWeight); // void??
  theHistograms.fill<TH1I>(type+"_pass_puCutBased_tight"  , 2, 0, 2, jet.pass_puCutBased_tight  (), theWeight); // void??
}


