/** \class EventAnalyzer
 *  Base class for event analyzers. Analyzers have to inherit from this class 
 *  and implement the pure virtual function analyze(), called each event.
 *
 *  $Date: 2013/03/15 13:37:31 $
 *  $Revision: 1.4 $
 *  \author R. Bellan - UNITO <riccardo.bellan@cern.ch>
 */

#include "EventAnalyzer.h"

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



// ------------------------------------------------------------------------------------------ //
// ------------------------------- C'tor/Des'tor/Init the tree ------------------------------ //
// ------------------------------------------------------------------------------------------ //



EventAnalyzer::EventAnalyzer(std::string filename, 
			     double lumi, 
			     double externalXSection, bool doBasicPlots)
  : doBasicPlots_(doBasicPlots)
  , theWeight(1.)
  , theSampleWeight(1.)
  , theCutCounter(0)
  , puweight(1.){

  TChain *tree = new TChain("treePlanter/ElderTree");
  tree->Add(filename.c_str());

  if (tree == 0) std::cout<<"Error in EventAnalyzer ctor"<<std::endl;

  //theSampleWeight = extractMCWeight(filename, lumi, externalXSection);
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
  jets = 0;      b_jets = 0;      theTree->SetBranchAddress("jets", &jets, &b_jets);

  // Jets   
  genParticles = 0; b_genParticles = 0; theTree->SetBranchAddress("genParticles", &genParticles, &b_genParticles);
  
  // MET
  met = new phys::Particle();
  b_met         = 0; theTree->SetBranchAddress("met"    , &met    ,  &b_met    );
  
  // Other events variables
  b_event       = 0; theTree->SetBranchAddress("event"    , &event    , &b_event    );
  b_run         = 0; theTree->SetBranchAddress("run"      , &run      , &b_run      );
  b_lumiBlock   = 0; theTree->SetBranchAddress("lumiBlock", &lumiBlock, &b_lumiBlock);

  b_nvtx        = 0; theTree->SetBranchAddress("nvtxs"  , &nvtx   , &b_nvtx   );
  b_rho         = 0; theTree->SetBranchAddress("rho"    , &rho    , &b_rho    );
  b_weight      = 0; theTree->SetBranchAddress("weight" , &weight , &b_weight );
  
  
  // MC related variables
  cout << "Weight from the event sample type: " << theSampleWeight << ", total weight (including PU reweight, if applicable): " << theWeight << endl;
  if(theSampleWeight != 1) b_puweight    = 0; theTree->SetBranchAddress("puweight" , &puweight , &b_puweight );
  b_xsec        = 0; theTree->SetBranchAddress("xsec"   , &xsec   , &b_xsec  );
  b_genCategory = 0; theTree->SetBranchAddress("genCategory"   , &genCategory   , &b_genCategory  );
   
   // Global variables of the sample
   //b_totalEvents = 0; theTree->SetBranchAddress("totalEvents", &totalEvents  , &b_totalEvents  );
}



// ------------------------------------------------------------------------------------------ //
// ------------------------------------- Master the loop ------------------------------------ //
// ------------------------------------------------------------------------------------------ //



Int_t EventAnalyzer::GetEntry(Long64_t entry){
  // Read contents of entry.
      
  if (!theTree) return 0;
  
  int e =  theTree->GetEntry(entry);
  
  stable_sort(muons->begin(),     muons->end(),     PtComparator());
  stable_sort(electrons->begin(), electrons->end(), PtComparator());
  stable_sort(jets->begin(),      jets->end(),      PtComparator());

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
  Long64_t nbytes = 0, nb = 0;

  begin();

  for (Long64_t jentry=0; jentry<nentries;jentry++) {

    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = GetEntry(jentry);  nbytes += nb; 

    theWeight = theSampleWeight;
    if(theWeight != 1) theWeight *= puweight;

    if (cut() < 0) continue;

    if(doBasicPlots_) fillBasicPlots();
    analyze();
  }

  TFile fout(TString(outputfile),"RECREATE");
  fout.cd(); 
  histograms.write(fout);

  end(fout);
  fout.Close();
  cout<<"Event passing all cuts: " << theCutCounter << endl;
}



// ------------------------------------------------------------------------------------------ //
// -------------------------------------- Preselection -------------------------------------- //
// ------------------------------------------------------------------------------------------ //



Int_t EventAnalyzer::cut(){
  
  bool pass = true;

  if(pass) ++theCutCounter;

  return pass ? 1 : -1;
}



// ------------------------------------------------------------------------------------------ //
// --------------------------------------- Histograms --------------------------------------- //
// ------------------------------------------------------------------------------------------ //



void EventAnalyzer::fillBasicPlots(){
  histograms.fill<TH1I>("nvtx"     , "Number of vertices" , 100, 0, 100, nvtx             , theWeight);
  histograms.fill      ("rho"      , "Mean energy density", 100, 0, 50 , rho              , theWeight);
  histograms.fill      ("met"      , "Missing energy"     , 200, 0, 800, met->pt()        , theWeight);
  histograms.fill<TH1I>("nmuons"    ,"Number of muons"    ,  10, 0, 10 , muons->size()    , theWeight);
  histograms.fill<TH1I>("nelectrons","Number of electrons",  10, 0, 10 , electrons->size(), theWeight);

  foreach(const phys::Lepton& mu , *muons)     fillLeptonPlots  ("mu",  mu  );
  foreach(const phys::Electron& e, *electrons) fillElectronPlots("e" ,  e   );
  foreach(const phys::Jet& jet   , *jets)      fillJetPlots     ("jet", jet );
}



void EventAnalyzer::fillParticlePlots(const std::string &type, const phys::Particle & lepton){
  histograms.fill(type+"_pt" , "p_{T} spectrum", 100,   0   , 500   ,lepton.pt(),   theWeight);
  histograms.fill(type+"_eta", "#eta spectrum" , 100,  -2.5 ,   2.5 ,lepton.eta(), theWeight);
  histograms.fill(type+"_phi", "#phi spectrum" ,  50,  -3.15,   3.15,lepton.phi(), theWeight);
}



void EventAnalyzer::fillLeptonPlots(const std::string &type, const phys::Lepton & lepton){
  
  fillParticlePlots(type, lepton);

  histograms.fill(type+"_dxy"             , "d_{xy}"         , 200,   0,   0.5, lepton.dxy()            , theWeight);   
  histograms.fill(type+"_dz"              , "d_{z}"          , 200,   0,   0.5, lepton.dz()             , theWeight);
  histograms.fill(type+"_sip"             , "sip"            , 100, -10,  10  , lepton.sip()            , theWeight); 
  histograms.fill(type+"_combRelIso"      , "combRelIso"     , 200,   0,   1  , lepton.combRelIso()     , theWeight);       
  histograms.fill(type+"_pfCombRelIso"    , "pfCombRelIso"   , 200,   0,   1  , lepton.pfChargedHadIso(), theWeight);           
  histograms.fill(type+"_pfNeutralHadIso" , "pfNeutralHadIso", 200,   0,   1  , lepton.pfNeutralHadIso(), theWeight);        
  histograms.fill(type+"_pfChargedHadIso" , "pfChargedHadIso", 200,   0,   1  , lepton.pfPhotonIso()    , theWeight);        
  histograms.fill(type+"_pfPhotonHadIso"  , "pfPhotonHadIso" , 200,   0,   1  , lepton.pfCombRelIso()   , theWeight);         
  histograms.fill(type+"_rho"             , "rho"            , 200,   0, 200  , lepton.rho()            , theWeight);                    
}



void EventAnalyzer::fillElectronPlots (const std::string &type, const phys::Electron &electron){
  fillLeptonPlots           (type, electron);
  fillExtraPlotsForElectrons(type, electron);
}



void EventAnalyzer::fillExtraPlotsForElectrons(const std::string &type, const phys::Electron &electron){

  histograms.fill      (type+"_energy"    , "energy spectrum"   , 100,  0   , 500   , electron.energy()    , theWeight);  
  histograms.fill      (type+"_phiWidth"  , "#phi width"        ,  50, -3.15,   3.15, electron.phiWidth()  , theWeight);  
  histograms.fill      (type+"_etaWidth"  , "#eta width"        , 100, -2.5 ,   2.5 , electron.etaWidth()  , theWeight);  
  histograms.fill      (type+"_BDT"       , "BDT"               , 100, -1   ,   1   , electron.BDT()       , theWeight);  
  histograms.fill<TH1I>(type+"_isBDT"     , "is BDT?"           ,   2,  0   ,   2   , electron.isBDT()     , theWeight);  
  histograms.fill<TH1I>(type+"_missingHit", "Missing hits"      ,  50,  0   ,  50   , electron.missingHit(), theWeight);  
  histograms.fill<TH1I>(type+"_nCrystals" , "Number of Crystals",  50,  0   ,  50   , electron.nCrystals() , theWeight);  
}



void EventAnalyzer::fillJetPlots(const std::string &type, const phys::Jet      &jet){

  histograms.fill<TH1I>(type+"_nConstituents", 100, 0, 100,jet.nConstituents(),theWeight);
  histograms.fill<TH1I>(type+"_nCharged"     , 100, 0, 100,jet.nCharged(),theWeight);
  histograms.fill<TH1I>(type+"_nNeutral"     , 100, 0, 100,jet.nNeutral(),theWeight);

  histograms.fill(type+"_neutralHadronEnergyFraction", 100, 0, 1, jet.neutralHadronEnergyFraction(), theWeight);
  histograms.fill(type+"_chargedHadronEnergyFraction", 100, 0, 1, jet.chargedHadronEnergyFraction(), theWeight);
  histograms.fill(type+"_chargedEmEnergyFraction"    , 100, 0, 1, jet.chargedEmEnergyFraction    (), theWeight);
  histograms.fill(type+"_neutralEmEnergyFraction"    , 100, 0, 1, jet.neutralEmEnergyFraction    (), theWeight);
  histograms.fill(type+"_muonEnergyFraction"         , 100, 0, 1, jet.muonEnergyFraction         (), theWeight);

  histograms.fill(type+"_csvtagger"    , 200, -1, 1 , jet.csvtagger    (), theWeight);
  histograms.fill(type+"_girth"        , 200,  0, 1 , jet.girth        (), theWeight);
  histograms.fill(type+"_girth_charged", 200,  0, 1 , jet.girth_charged(), theWeight);
  histograms.fill(type+"_ptd"          , 1, 1, 1    , jet.ptd          (), theWeight);
  histograms.fill(type+"_rms"          , 1, 1, 1    , jet.rms          (), theWeight);
  histograms.fill(type+"_beta"         , 1, 1, 1    , jet.beta         (), theWeight);
  histograms.fill(type+"_jetArea"      , 1, 1, 1    , jet.jetArea      (), theWeight);
  histograms.fill(type+"_secvtxMass"   , 1, 1, 1    , jet.secvtxMass   (), theWeight);
  histograms.fill(type+"_Lxy"          , 1, 1, 1    , jet.Lxy          (), theWeight);
  histograms.fill(type+"_LxyErr"       , 1, 1, 1    , jet.LxyErr       (), theWeight);
  histograms.fill(type+"_rawFactor"    , 1, 1, 1    , jet.rawFactor    (), theWeight);

  histograms.fill(type+"_uncOnFourVectorScale", 1, 1, 1, jet.uncOnFourVectorScale(), theWeight);
  histograms.fill(type+"_puMVAFull"  , 1, 1, 1, jet.puMVAFull  (), theWeight);
  histograms.fill(type+"_puMVASimple", 1, 1, 1, jet.puMVASimple(), theWeight);
  histograms.fill(type+"_puCutBased" , 1, 1, 1, jet.puCutBased (), theWeight);

  histograms.fill<TH1I>(type+"_pass_puMVAFull_loose"   , 2, 0, 2, jet.pass_puMVAFull_loose   (), theWeight);
  histograms.fill<TH1I>(type+"_pass_pUMVAFull_medium"  , 2, 0, 2, jet.pass_pUMVAFull_medium  (), theWeight);
  histograms.fill<TH1I>(type+"_pass_pUMVAFull_tight"   , 2, 0, 2, jet.pass_pUMVAFull_tight   (), theWeight);
				                          	                             
  histograms.fill<TH1I>(type+"_pass_puMVASimple_loose" , 2, 0, 2, jet.pass_puMVASimple_loose (), theWeight);
  histograms.fill<TH1I>(type+"_pass_puMVASimple_medium", 2, 0, 2, jet.pass_puMVASimple_medium(), theWeight);
  histograms.fill<TH1I>(type+"_pass_puMVASimple_tight" , 2, 0, 2, jet.pass_puMVASimple_tight (), theWeight);
				                          	                             
  histograms.fill<TH1I>(type+"_pass_puCutBased_loose"  , 2, 0, 2, jet.pass_puCutBased_loose  (), theWeight);
  histograms.fill<TH1I>(type+"_pass_puCutBased_medium" , 2, 0, 2, jet.pass_puCutBased_medium (), theWeight);
  histograms.fill<TH1I>(type+"_pass_puCutBased_tight"  , 2, 0, 2, jet.pass_puCutBased_tight  (), theWeight);
}
