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

#include <algorithm>

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

SelectorBase::~SelectorBase(){};


EventAnalyzer::EventAnalyzer(SelectorBase& aSelector,
			     const AnalysisConfiguration& configuration)
  : select(aSelector)
  , maxNumEvents_(configuration.getParameter<int>("maxNumEvents"))
  , doBasicPlots_(configuration.getParameter<bool>("doBasicPlots"))
  , doSF         (configuration.getParameter<bool>("doSF"))
  , regions_     (configuration.getParameter<std::vector<phys::RegionTypes>>("regions"))
  , theSampleInfo(configuration.getParameter<std::string>("filename"), 
		  configuration.getParameter<double>("lumi"), 
		  configuration.getParameter<double>("externalXSection"),
		  configuration.getParameter<bool>("blinded"),
		  configuration.getParameter<bool>("applyFRSF"),
		  configuration.getParameter<bool>("forcePosWeight")
		  )
  , theWeight(1.)
  , theCutCounter(0.)
  , theInputWeightedEvents(0.)
  , genCategory(-128){

  if(configuration.getParameter<int>("year") != theSampleInfo.setup() && theSampleInfo.isMC())
    cout << colour::Warning("Possible mismatch") << ": simulation scenario is " << Green(configuration.getParameter<int>("year")) << ", chosen sample is " << Green(theSampleInfo.setup()) << endl;

  TChain *tree = new TChain("treePlanter/ElderTree");
  tree->Add(configuration.getParameter<std::string>("filename").c_str());

  if (tree == 0) std::cout<<Important("Error in EventAnalyzer ctor:")<<" The tree has a null pointer."<<std::endl;

  Init(tree);
  fileName =  configuration.getParameter<std::string>("filename");

  cout<<Yellow("Analyzing "+fileName+" ... please wait... \n")<<endl ;
}


EventAnalyzer::~EventAnalyzer(){
  delete theTree->GetCurrentFile();
}


void EventAnalyzer::Init(TTree *tree)
{   
  // Set branch addresses and branch pointers
  if (!tree) return;
  theTree = tree;
  fCurrent = -1;

  Long64_t selectedEvents = theTree->GetEntries();
  std::cout<<"\nThis sample contains " << Green(selectedEvents) << " selected events.\n"   << std::endl;

  
  // Muons   
  muons     = 0; b_muons     = 0; theTree->SetBranchAddress("muons"    , &muons    , &b_muons    );
  
  // Electrons   
  electrons = 0; b_electrons = 0; theTree->SetBranchAddress("electrons", &electrons, &b_electrons);
  
  // Jets   
  pjets = 0;      b_pjets = 0;    theTree->SetBranchAddress("jets", &pjets, &b_pjets);
  jets  = new std::vector<phys::Jet>(); centralJets  = new std::vector<phys::Jet>();

  pjetsAK8 = 0;      b_pjetsAK8 = 0;    theTree->SetBranchAddress("jetsAK8", &pjetsAK8, &b_pjetsAK8);
  jetsAK8  = new std::vector<phys::Jet>();
  
  // Photons
  photons = 0;      b_photons = 0;    theTree->SetBranchAddress("photons", &photons, &b_photons);
  jetsAK8  = new std::vector<phys::Jet>();

  // Bosons   
  Vhad = new std::vector<phys::Boson<phys::Jet> > ()    ; VhadCand = 0; b_VhadCand = 0; theTree->SetBranchAddress("VhadCand", &VhadCand, &b_VhadCand);

  // DiBoson, if in SR, or Z+ll if in CR
  ZZ   = new phys::DiBoson<phys::Lepton, phys::Lepton>(); b_ZZ   = 0; theTree->SetBranchAddress("ZZCand"  , &ZZ  , &b_ZZ  );

  // DiBoson, if in SR, or Z+ll if in CR
  ZW   = new phys::DiBoson<phys::Lepton, phys::Lepton>(); b_ZW   = 0; theTree->SetBranchAddress("ZWCand"  , &ZW  , &b_ZW  );


  // Z->ll
  Z   = new phys::Boson<phys::Lepton>(); b_Z   = 0; theTree->SetBranchAddress("ZCand"  , &Z  , &b_Z  );


  // Z+L 
  ZL = new ZLCompositeCandidate()    ; b_ZLCand = 0; theTree->SetBranchAddress("ZLCand", &ZL, &b_ZLCand);


  // Gen Particles   
  genParticles   = 0;                                                b_genParticles   = 0; theTree->SetBranchAddress("genParticles"  , &genParticles  , &b_genParticles);
  genTaus        = 0;                                                b_genTaus        = 0; theTree->SetBranchAddress("genTaus"       , &genTaus       , &b_genTaus);
  genVBParticles = new std::vector<phys::Boson<phys::Particle> > (); b_genVBParticles = 0; theTree->SetBranchAddress("genVBParticles", &genVBParticles, &b_genVBParticles);
  
  // Gen Jets
  pgenJets   = 0;                                                    b_pgenJets   = 0; theTree->SetBranchAddress("genJets"  , &pgenJets  , &b_pgenJets);
  genJets  = new std::vector<phys::Particle>(); centralGenJets  = new std::vector<phys::Particle>();

  pgenJetsAK8 = 0;                                                   b_pgenJetsAK8 = 0; theTree->SetBranchAddress("genJetsAK8"  , &pgenJetsAK8  , &b_pgenJetsAK8);
  genJetsAK8  = new std::vector<phys::Particle>();


  // MET
  met = new phys::Particle();
  b_met         = 0; theTree->SetBranchAddress("met"    , &met    ,  &b_met    );
  
  // Other events variables
  b_event       = 0; theTree->SetBranchAddress("event"    , &event    , &b_event    );
  b_run         = 0; theTree->SetBranchAddress("run"      , &run      , &b_run      );
  b_lumiBlock   = 0; theTree->SetBranchAddress("lumiBlock", &lumiBlock, &b_lumiBlock);

  b_nvtx        = 0; theTree->SetBranchAddress("nvtxs"  , &nvtx   , &b_nvtx   );
  
  b_genCategory  = 0; theTree->SetBranchAddress("genCategory" , &genCategory             , &b_genCategory  );

  // MC related variables
  b_genEventWeights = 0; theTree->SetBranchAddress("genEventWeights", &theSampleInfo.genEventWeights_ , &b_genEventWeights);

  // Info about selections
  b_passTrigger = 0; theTree->SetBranchAddress("passTrigger", &passTrigger, &b_passTrigger); 
  b_passSkim    = 0; theTree->SetBranchAddress("passSkim"   , &passSkim   , &b_passSkim   ); 
  b_triggerWord = 0; theTree->SetBranchAddress("triggerWord", &triggerWord, &b_triggerWord); 
  b_regionWord = 0; theTree->SetBranchAddress("regionWord", &pregionWord, &b_regionWord); 
  
  // MELA discriminators
  mela = new phys::MELA();
  b_mela = 0;  theTree->SetBranchAddress("MELA"    , &mela    ,  &b_mela    );

}

// ------------------------------------------------------------------------------------------ //
// ------------------------------------- Master the loop ------------------------------------ //
// ------------------------------------------------------------------------------------------ //



Int_t EventAnalyzer::GetEntry(Long64_t entry){
  // Read contents of entry.
  if (!theTree) return 0;
  
  int e =  theTree->GetEntry(entry);
  
  // Check if the request on region tye matches with the categorization of the event
  regionWord = std::bitset<32>(pregionWord);

  // FIXME: rise _1P or _1F bits
  if(std::any_of(regions_.begin(), regions_.end(), [](phys::RegionTypes r){ return r==phys::MC; })){
    region_ = phys::MC;
    goto _continueEvent;
  }
  for(phys::RegionTypes region : regions_){
    if(regionWord.test(region)){
      region_ = region;
      goto _continueEvent;
    }
  }
  return 0;
  
  _continueEvent:
  theHistograms = &(mapRegionHisto_[region_]);
  

  if(muons)       stable_sort(muons->begin(),       muons->end(),       phys::PtComparator());
  if(electrons)   stable_sort(electrons->begin(),   electrons->end(),   phys::PtComparator());
  if(pjets)       stable_sort(pjets->begin(),       pjets->end(),       phys::PtComparator());
  if(pgenJets)    stable_sort(pgenJets->begin(),    pgenJets->end(),    phys::PtComparator());
  if(pjetsAK8)    stable_sort(pjetsAK8->begin(),    pjetsAK8->end(),    phys::PtComparator());
  if(pgenJetsAK8) stable_sort(pgenJetsAK8->begin(), pgenJetsAK8->end(), phys::PtComparator());
  if(photons)     stable_sort(photons->begin(),   photons->end(),   phys::PtComparator());
	
  // Some selection on jets
  jets->clear(); centralJets->clear(); 
  if(pjets){
    foreach(const phys::Jet &jet, *pjets)
      if(jet.pt() > 30){
	if(fabs(jet.eta()) < 4.7) jets->push_back(jet);
	if(fabs(jet.eta()) < 2.4) centralJets->push_back(jet);
      }
  }
  genJets->clear(); centralGenJets->clear();
  if(pgenJets){
    foreach(const phys::Particle &jet, *pgenJets)
      if(jet.pt() > 30){
	if(fabs(jet.eta()) < 4.7) genJets->push_back(jet);
	if(fabs(jet.eta()) < 2.4) centralGenJets->push_back(jet);
      }
  }
  genVBHelper_.analyze(*genParticles, *genVBParticles);

  // Some selection on jets
  jetsAK8->clear();
  if(pjetsAK8){
    foreach(const phys::Jet &jet, *pjetsAK8)
      if(jet.pt() > 30 && fabs(jet.eta()) < 4.7) jetsAK8->push_back(jet);
  }
  genJetsAK8->clear();
  if(pgenJetsAK8){
    foreach(const phys::Particle &jet, *pgenJetsAK8)
      if(jet.pt() > 30 && fabs(jet.eta()) < 4.7) genJetsAK8->push_back(jet);
  }  
  
  Vhad->clear();
  if(VhadCand){
    foreach(const phys::Boson<phys::Jet> v, *VhadCand)
      if(select(v)) Vhad->push_back(v);
  }
  stable_sort(Vhad->begin(), Vhad->end(), phys::PtComparator());
  
  
  if(region_ == phys::MC){
    if(!ZZ->isValid()){
      if(ZZ) delete ZZ;
      ZZ = new phys::DiBoson<phys::Lepton, phys::Lepton>();
    }
    if(!ZW->isValid()){ 
      if(ZW) delete ZW;
      ZW = new phys::DiBoson<phys::Lepton, phys::Lepton>();
    }
    if(!Z->isValid()){ 
      if(Z) delete Z;
      Z = new phys::Boson<phys::Lepton>();
    }
    if(!ZL->first.isValid()){ 
      if(ZL) delete ZL;
      ZL = new std::pair<phys::Boson<phys::Lepton>, phys::Lepton>();
    }

  }
  
  
  addOptions();
  
  if(region_ < phys::SR3P)
    theWeight = theSampleInfo.weight(*ZZ);

  else if(region_ >= phys::SR3P && region_ < phys::CRLFR)
    theWeight = theSampleInfo.weight(*ZW);

  else if(region_ == phys::CRLFR)
    theWeight = theSampleInfo.weight(ZL->first);
  
  else if(region_ > phys::CRLFR && region_ <= phys::CR2P_1F)
    theWeight = theSampleInfo.weight(*Z);
  
  else{
    if(region_ != phys::MC){
      std::cout<<"Do not know what weight to set. Aborting... "  << endl;
      std::abort();
    }
  }
    

  theHistograms->fill("weight_full"  , "All weights applied"                                    , 1200, -2, 10, theWeight);
  theHistograms->fill("weight_bare"  , "All weights, but efficiency and fake rate scale factors", 1200, -2, 10, theSampleInfo.weight());
  theHistograms->fill("weight_pu"    , "Weight from PU reweighting procedure"                   , 1200, -2, 10, theSampleInfo.puWeight());
  theHistograms->fill("weight_sample", "Weight from cross-section and luminosity"               , 1200, -2, 10, theSampleInfo.sampleWeight());
  theHistograms->fill("weight_mcProc", "Weight from MC intrinsic event weight"                  , 1200, -2, 10, theSampleInfo.mcWeight());
  theHistograms->fill("weight_efficiencySF", "Weight from data/MC lepton efficiency"            , 1200, -2, 10, ZZ->efficiencySF());
  theHistograms->fill("weight_fakeRateSF"  , "Weight from fake rate scale factor"               , 1200, -2, 10, ZZ->fakeRateSF());
  
  
  theInputWeightedEvents += theWeight;
  
  topology = std::bitset<16>(genCategory);
  
  
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


  Long64_t nentries = maxNumEvents_ > 0 ? maxNumEvents_ : theTree->GetEntries();  

  begin();

  for (Long64_t jentry=0; jentry<nentries; ++jentry) {
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    if (!GetEntry(jentry)) continue;
    if (cut() < 0) continue;
    theCutCounter += theWeight;
    if(doBasicPlots_) fillBasicPlots();
    analyze();
  }
  
  for(std::pair<phys::RegionTypes, Histogrammer> regHist : mapRegionHisto_){
    std::string regionName = phys::regionType(regHist.first);
    region_ = regHist.first;  // In case some Analyzer wants to do something special in certain regions
    theHistograms = & regHist.second;  // So that any operation on the histograms is done on the correct ones
    
    std::string fout_formatted( Form(outputfile.c_str(), regionName.c_str()) );
    TFile fout(TString(fout_formatted),"RECREATE");
    fout.cd();
    
    end(fout);
    theHistograms->write(fout);
    
    fout.Close();
  }
  
  finish();
  
  std::string regionsString;
  for(phys::RegionTypes r : regions_) regionsString += regionType(r)+';';
  regionsString.pop_back();
  cout<<"Analyzed events in the chosen region(s) (" << Blue(regionsString) << "): " << Green(theInputWeightedEvents)<< endl;

  cout<<"Events passing all cuts: "<< Green(theCutCounter) << endl; // FIXME
}



// ------------------------------------------------------------------------------------------ //
// -------------------------------------- Preselection -------------------------------------- //
// ------------------------------------------------------------------------------------------ //


Int_t EventAnalyzer::cut() {
  
  bool pass = true;
  
  return pass ? 1 : -1;
}

//
//
//

// void EventAnalyzer::applyLeptonScaleFactors(){

//   // Protection
//   if(!doSF) return;

//   theWeight = theSampleInfo.weight();
  
//   if(region_ == phys::CR2P2F || region_ == phys::CR3P1F || region_ == phys::CR2P2F_HZZ || region_ == phys::CR3P1F_HZZ){
    
    
//     if(!ZZ->first().daughterPtr(0)->passFullSel())   theWeight*= (leptonScaleFactors_.fakeRateScaleFactor(*ZZ->first().daughterPtr(0))).first;	
//     if(!ZZ->first().daughterPtr(1)->passFullSel())   theWeight*= (leptonScaleFactors_.fakeRateScaleFactor(*ZZ->first().daughterPtr(1))).first;	
//     if(!ZZ->second().daughterPtr(0)->passFullSel())  theWeight*= (leptonScaleFactors_.fakeRateScaleFactor(*ZZ->second().daughterPtr(0))).first;
//     if(!ZZ->second().daughterPtr(1)->passFullSel())  theWeight*= (leptonScaleFactors_.fakeRateScaleFactor(*ZZ->second().daughterPtr(1))).first; 
    
//     if(region_ == phys::CR2P2F || region_ == phys::CR2P2F_HZZ ) theWeight*=-1;
    
//   }
  
//   if(theSampleInfo.isMC()){
    
    
//     std::pair<double,double> lepSF;
    
//     lepSF=leptonScaleFactors_.efficiencyScaleFactor(*ZZ->first().daughterPtr(0));
//     (ZZ->first().daughterPtr(0))->setEfficenySFUnc(lepSF.second);
//     theWeight*=lepSF.first;
    
//     lepSF=leptonScaleFactors_.efficiencyScaleFactor(*ZZ->first().daughterPtr(1));
//     (ZZ->first().daughterPtr(1))->setEfficenySFUnc(lepSF.second);
//     theWeight*=lepSF.first;
    
//     lepSF=leptonScaleFactors_.efficiencyScaleFactor(*ZZ->second().daughterPtr(0));
//     ( ZZ->second().daughterPtr(0))->setEfficenySFUnc(lepSF.second);
//     theWeight*=lepSF.first;
    
//     lepSF=leptonScaleFactors_.efficiencyScaleFactor(*ZZ->second().daughterPtr(1));
//     (ZZ->second().daughterPtr(1))->setEfficenySFUnc(lepSF.second);
//     theWeight*=lepSF.first;
//   }

// }





// ------------------------------------------------------------------------------------------ //
// --------------------------------------- Histograms --------------------------------------- //
// ------------------------------------------------------------------------------------------ //



void EventAnalyzer::fillBasicPlots(){
  theHistograms->fill<TH1I>("nvtx"     , "Number of vertices" , 100, 0, 100, nvtx             , theWeight);
  theHistograms->fill      ("met"      , "Missing energy"     , 200, 0, 800, met->pt()        , theWeight);
  theHistograms->fill<TH1I>("nmuons"    ,"Number of muons"    ,  10, 0, 10 , muons->size()    , theWeight);
  theHistograms->fill<TH1I>("nelectrons","Number of electrons",  10, 0, 10 , electrons->size(), theWeight);

  foreach(const phys::Lepton& mu , *muons)     fillLeptonPlots("mu",  mu  );
  foreach(const phys::Electron& e, *electrons) fillLeptonPlots("e" ,  e   );
  foreach(const phys::Jet& jet   , *jets)      fillJetPlots   ("jet", jet );
}



void EventAnalyzer::fillParticlePlots(const std::string &type, const phys::Particle & lepton){
  theHistograms->fill(type+"_pt" ,    "p_{T} spectrum", 100,   0   , 500   ,lepton.pt()    , theWeight);
  theHistograms->fill(type+"_eta",    "#eta spectrum" , 100,  -2.5 ,   2.5 ,lepton.eta()   , theWeight);
  theHistograms->fill(type+"_phi",    "#phi spectrum" ,  50,  -3.15,   3.15,lepton.phi()   , theWeight);
  theHistograms->fill(type+"_charge", "charge"        ,  50,  -25  ,  25   ,lepton.charge(), theWeight);
}



void EventAnalyzer::fillLeptonPlots(const std::string &type, const phys::Lepton & lepton){
  
  fillParticlePlots(type, lepton);

  theHistograms->fill(type+"_dxy"             , "d_{xy}"         , 200,   0,   0.5, lepton.dxy()            , theWeight);
  theHistograms->fill(type+"_dz"              , "d_{z}"          , 200,   0,   0.5, lepton.dz()             , theWeight);
  theHistograms->fill(type+"_sip"             , "sip"            , 150,   0,  15  , lepton.sip()            , theWeight);
}



void EventAnalyzer::fillJetPlots(const std::string &type, const phys::Jet      &jet){

  fillParticlePlots(type, jet);

  theHistograms->fill(type+"_csvtagger"    , 200,  0, 1   , jet.csvtagger    (), theWeight);
  theHistograms->fill(type+"_girth"        , 200,  0, 1   , jet.girth        (), theWeight);
  theHistograms->fill(type+"_girth_charged", 200,  0, 1   , jet.girth_charged(), theWeight);
  theHistograms->fill(type+"_ptd"          , 100,  0, 1   , jet.ptd          (), theWeight);
  theHistograms->fill(type+"_jetArea"      ,  60,  0, 1.2 , jet.jetArea      (), theWeight);
  theHistograms->fill(type+"_secvtxMass"   ,  50,  0, 5   , jet.secvtxMass   (), theWeight);
  theHistograms->fill(type+"_rawFactor"    ,  50,  0, 2.5 , jet.rawFactor    (), theWeight);

  theHistograms->fill(type+"_jecUncertainty", 50, 0, 0.1, jet.jecUncertainty(), theWeight);
}


