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


EventAnalyzer::EventAnalyzer(SelectorBase& aSelector,
			     const AnalysisConfiguration& configuration)
  : select(aSelector)
  ,  leptonScaleFactors_("../../ZZAnalysis/AnalysisStep/data/LeptonEffScaleFactors/ScaleFactors_mu_Moriond2017.root",
			 "../../ZZAnalysis/AnalysisStep/data/LeptonEffScaleFactors/ScaleFactors_non_gap_ele_Moriond2017_v2.root",
			 "../../ZZAnalysis/AnalysisStep/data/LeptonEffScaleFactors/ScaleFactors_gap_ele_Moriond2017_v2.root",
			 "../../ZZAnalysis/AnalysisStep/data/LeptonEffScaleFactors/ScaleFactors_RECO_ele_Moriond2017_v1.root",
			 "../../VVXAnalysis/Commons/data/fakeRates.root",
			 "../../VVXAnalysis/Commons/data/fakeRates.root")
  , doBasicPlots_(configuration.getParameter<bool>("doBasicPlots"))
  , doSF         (configuration.getParameter<bool>("doSF"))
  , region_      (configuration.getParameter<phys::RegionTypes>("region"))
  , theMCInfo    (configuration.getParameter<std::string>("filename"), 
		  configuration.getParameter<double>("lumi"), 
		  configuration.getParameter<double>("externalXSection"))
  , theWeight(1.)
  , theCutCounter(0.)
  , theInputWeightedEvents(0.)
  , unweightedEventsInSR(0)
  , unweightedEventsIn2P2FCR(0)
  , unweightedEventsIn3P1FCR(0)
  , genCategory(-128){

  TChain *tree = new TChain("treePlanter/ElderTree");
  tree->Add(configuration.getParameter<std::string>("filename").c_str());


  if (tree == 0) std::cout<<Important("Error in EventAnalyzer ctor:")<<" The tree has a null pointer."<<std::endl;

  Init(tree);

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

  // Muons   
  muons     = 0; b_muons     = 0; theTree->SetBranchAddress("muons"    , &muons    , &b_muons    );
  
  // Electrons   
  electrons = 0; b_electrons = 0; theTree->SetBranchAddress("electrons", &electrons, &b_electrons);
  
  // Jets   
  pjets = 0;      b_pjets = 0;    theTree->SetBranchAddress("jets", &pjets, &b_pjets);
  jets  = new std::vector<phys::Jet>(); centralJets  = new std::vector<phys::Jet>();
  pileUpIds  = new std::vector<Int_t>();

  // Bosons   
  Vhad = new std::vector<phys::Boson<phys::Jet> > ()    ; VhadCand = 0; b_VhadCand = 0; theTree->SetBranchAddress("VhadCand", &VhadCand, &b_VhadCand);

  // DiBoson, if in SR, or Z+ll if in CR
  ZZ   = new phys::DiBoson<phys::Lepton, phys::Lepton>(); b_ZZ   = 0; theTree->SetBranchAddress("ZZCand"  , &ZZ  , &b_ZZ  );

  // Z+L 
  ZL = new ZLCompositeCandidates()    ; ZLCand = 0; b_ZLCand = 0; theTree->SetBranchAddress("ZLCand", &ZLCand, &b_ZLCand);


  // Gen Particles   
  genParticles   = 0;                                                b_genParticles   = 0; theTree->SetBranchAddress("genParticles"  , &genParticles  , &b_genParticles);
  genVBParticles = new std::vector<phys::Boson<phys::Particle> > (); b_genVBParticles = 0; theTree->SetBranchAddress("genVBParticles", &genVBParticles, &b_genVBParticles);
  
  // Gen Jets
  pgenJets   = 0;                                                    b_pgenJets   = 0; theTree->SetBranchAddress("genJets"  , &pgenJets  , &b_pgenJets);
  genJets  = new std::vector<phys::Particle>(); centralGenJets  = new std::vector<phys::Particle>();

  // MET
  met = new phys::Particle();
  b_met         = 0; theTree->SetBranchAddress("met"    , &met    ,  &b_met    );
  
  // Other events variables
  b_event       = 0; theTree->SetBranchAddress("event"    , &event    , &b_event    );
  b_run         = 0; theTree->SetBranchAddress("run"      , &run      , &b_run      );
  b_lumiBlock   = 0; theTree->SetBranchAddress("lumiBlock", &lumiBlock, &b_lumiBlock);

  b_nvtx        = 0; theTree->SetBranchAddress("nvtxs"  , &nvtx   , &b_nvtx   );
  

  // MC related variables
  b_puweight     = 0; theTree->SetBranchAddress("puweight"    , &theMCInfo.puweight_     , &b_puweight );
  b_puweightUp   = 0; theTree->SetBranchAddress("puweightUp"  , &theMCInfo.puweightUp_   , &b_puweightUp );
  b_puweightDn   = 0; theTree->SetBranchAddress("puweightDn"  , &theMCInfo.puweightDn_   , &b_puweightDn );
  b_mcprocweight = 0; theTree->SetBranchAddress("mcprocweight", &theMCInfo.mcprocweight_ , &b_mcprocweight);
  b_genCategory  = 0; theTree->SetBranchAddress("genCategory" , &genCategory             , &b_genCategory  );

  // Info about selections
  b_passTrigger = 0; theTree->SetBranchAddress("passTrigger", &passTrigger, &b_passTrigger); 
  b_passSkim    = 0; theTree->SetBranchAddress("passSkim"   , &passSkim   , &b_passSkim   ); 
  b_triggerWord = 0; theTree->SetBranchAddress("triggerWord", &triggerWord, &b_triggerWord); 
  
  //K factors
  
  b_kFactor_ggZZ     =0;  theTree->SetBranchAddress("kFactor_ggZZ"    ,  &kFactor_ggZZ          , &b_kFactor_ggZZ);
  b_kFactor_qqZZM    =0;  theTree->SetBranchAddress("kFactor_qqZZM"   ,  &kFactor_qqZZM         , &b_kFactor_qqZZM);
  b_kFactor_qqZZPt   =0;  theTree->SetBranchAddress("kFactor_qqZZPt"  ,  &kFactor_qqZZPt        , &b_kFactor_qqZZPt);
  b_kFactor_qqZZdPhi =0;  theTree->SetBranchAddress("kFactor_qqZZdPhi",  &kFactor_qqZZdPhi      , &b_kFactor_qqZZdPhi);
  b_kFactor_EWKqqZZ  =0;  theTree->SetBranchAddress("kFactor_EWKqqZZ" ,  &kFactor_EWKqqZZ       , &b_kFactor_EWKqqZZ);

  //theretican uncertainties

  b_LHEPDFScale                       =0; theTree->SetBranchAddress("LHEPDFScale"                     , &LHEPDFScale                      , &b_LHEPDFScale                     );
  b_LHEweight_QCDscale_muR1_muF1      =0; theTree->SetBranchAddress("LHEweight_QCDscale_muR1_muF1"    , &LHEweight_QCDscale_muR1_muF1     , &b_LHEweight_QCDscale_muR1_muF1    );
  b_LHEweight_QCDscale_muR1_muF2      =0; theTree->SetBranchAddress("LHEweight_QCDscale_muR1_muF2"    , &LHEweight_QCDscale_muR1_muF2     , &b_LHEweight_QCDscale_muR1_muF2    );
  b_LHEweight_QCDscale_muR1_muF0p5    =0; theTree->SetBranchAddress("LHEweight_QCDscale_muR1_muF0p5"  , &LHEweight_QCDscale_muR1_muF0p5   , &b_LHEweight_QCDscale_muR1_muF0p5  );
  b_LHEweight_QCDscale_muR2_muF1      =0; theTree->SetBranchAddress("LHEweight_QCDscale_muR2_muF1"    , &LHEweight_QCDscale_muR2_muF1     , &b_LHEweight_QCDscale_muR2_muF1    );
  b_LHEweight_QCDscale_muR2_muF2      =0; theTree->SetBranchAddress("LHEweight_QCDscale_muR2_muF2"    , &LHEweight_QCDscale_muR2_muF2     , &b_LHEweight_QCDscale_muR2_muF2    );
  b_LHEweight_QCDscale_muR2_muF0p5    =0; theTree->SetBranchAddress("LHEweight_QCDscale_muR2_muF0p5"  , &LHEweight_QCDscale_muR2_muF0p5   , &b_LHEweight_QCDscale_muR2_muF0p5  );
  b_LHEweight_QCDscale_muR0p5_muF1    =0; theTree->SetBranchAddress("LHEweight_QCDscale_muR0p5_muF1"  , &LHEweight_QCDscale_muR0p5_muF1   , &b_LHEweight_QCDscale_muR0p5_muF1  );
  b_LHEweight_QCDscale_muR0p5_muF2    =0; theTree->SetBranchAddress("LHEweight_QCDscale_muR0p5_muF2"  , &LHEweight_QCDscale_muR0p5_muF2   , &b_LHEweight_QCDscale_muR0p5_muF2  );
  b_LHEweight_QCDscale_muR0p5_muF0p5  =0; theTree->SetBranchAddress("LHEweight_QCDscale_muR0p5_muF0p5", &LHEweight_QCDscale_muR0p5_muF0p5 , &b_LHEweight_QCDscale_muR0p5_muF0p5);
  b_LHEweight_PDFVariation_Up         =0; theTree->SetBranchAddress("LHEweight_PDFVariation_Up"       , &LHEweight_PDFVariation_Up        , &b_LHEweight_PDFVariation_Up       );
  b_LHEweight_PDFVariation_Dn         =0; theTree->SetBranchAddress("LHEweight_PDFVariation_Dn"       , &LHEweight_PDFVariation_Dn        , &b_LHEweight_PDFVariation_Dn       );
  b_LHEweight_AsMZ_Up                 =0; theTree->SetBranchAddress("LHEweight_AsMZ_Up"               , &LHEweight_AsMZ_Up                , &b_LHEweight_AsMZ_Up               );
  b_LHEweight_AsMZ_Dn                 =0; theTree->SetBranchAddress("LHEweight_AsMZ_Dn"               , &LHEweight_AsMZ_Dn                , &b_LHEweight_AsMZ_Dn               );
  
  
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
  stable_sort(pgenJets->begin(),  pgenJets->end(),  phys::PtComparator());

  // Some selection on jets
  Int_t idx = 0;
  jets->clear(); centralJets->clear(); pileUpIds->clear();

  foreach(const phys::Jet &jet, *pjets){
    if(jet.fullPuId(0)){//2 tight, 1 medium, 0 loose
      if(jet.pt() > 30){
	if(fabs(jet.eta()) < 4.7) jets->push_back(jet);
	if(fabs(jet.eta()) < 2.4) centralJets->push_back(jet);
        }
    }
    else pileUpIds->push_back(idx);   
    idx++; 
  }
  genJets->clear(); centralGenJets->clear();
  foreach(const phys::Particle &jet, *pgenJets)
    if(jet.pt() > 30){
      //  bool leptonMatch = false; //Already done in signaldefinition. Here leptons are also loose. 
      //foreach(const phys::Particle &gen, *genParticles){
      //	if(physmath::deltaR(gen,jet) < 0.4 && (abs(gen.id()) == 11 || abs(gen.id()) == 13)) leptonMatch = true;
      //}
      // if(!leptonMatch){
	if(fabs(jet.eta()) < 4.7) genJets->push_back(jet);
	if(fabs(jet.eta()) < 2.4) centralGenJets->push_back(jet);
    }
  //}
Vhad->clear();

foreach(const phys::Boson<phys::Jet> v, *VhadCand)
if(select(v)) Vhad->push_back(v);

stable_sort(Vhad->begin(), Vhad->end(), phys::PtComparator());

ZL->clear();
foreach(const ZLCompositeCandidate& zl, *ZLCand)
if( zl.second.sip() < 4) ZL->push_back(zl); //To be moved in cfg
if(region_ == phys::MC){
  if(!ZZ->isValid()){
    if(ZZ) delete ZZ;
    ZZ = new phys::DiBoson<phys::Lepton, phys::Lepton>();
  }
 }

 std::vector<phys::Lepton> Leps;

 Leps.push_back(*ZZ->first().daughterPtr(0));
 Leps.push_back(*ZZ->first().daughterPtr(1));
 Leps.push_back(*ZZ->second().daughterPtr(0));
 Leps.push_back(*ZZ->second().daughterPtr(1));

 stable_sort(Leps.begin(),     Leps.end(),     phys::PtComparator());
 

 if((abs(Leps.at(1).id())==11) && (Leps.at(1).pt()<12))  return 0;
 
 addOptions();
 
  // Check if the request on region tye matches with the categorization of the event
  regionWord = std::bitset<128>(ZZ->region());
  // check bits accordingly to ZZAnalysis/AnalysisStep/interface/FinalStates.h

  if(region_ == phys::SR     && !regionWord.test(26)) return 0;
  if(region_ == phys::CR2P2F && !regionWord.test(24)) return 0;
  if(region_ == phys::CR3P1F && !regionWord.test(25)) return 0;

  if(region_ == phys::SR_HZZ     && !regionWord.test(3))  return 0;
  if(region_ == phys::CR2P2F_HZZ && !regionWord.test(22)) return 0;
  if(region_ == phys::CR3P1F_HZZ && !regionWord.test(23)) return 0;

  if(!doSF){
  theWeight = theMCInfo.weight(*ZZ);

  }
  else  theWeight = theMCInfo.weight();

  if(doSF && (region_ == phys::CR2P2F || region_ == phys::CR3P1F || region_ == phys::CR2P2F_HZZ || region_ == phys::CR3P1F_HZZ)){
      if(!ZZ->first().daughterPtr(0)->passFullSel())   theWeight*= (leptonScaleFactors_.fakeRateScaleFactor(*ZZ->first().daughterPtr(0))).first;	
      if(!ZZ->first().daughterPtr(1)->passFullSel())   theWeight*= (leptonScaleFactors_.fakeRateScaleFactor(*ZZ->first().daughterPtr(1))).first;	
      if(!ZZ->second().daughterPtr(0)->passFullSel())  theWeight*= (leptonScaleFactors_.fakeRateScaleFactor(*ZZ->second().daughterPtr(0))).first;
      if(!ZZ->second().daughterPtr(1)->passFullSel())  theWeight*= (leptonScaleFactors_.fakeRateScaleFactor(*ZZ->second().daughterPtr(1))).first; 
     
      if(region_ == phys::CR2P2F || region_ == phys::CR2P2F_HZZ ) theWeight*=-1;
    
    }

  if(doSF && theMCInfo.isMC() ){
    
    
    std::pair<double,double> lepSF;
    
    lepSF=leptonScaleFactors_.efficiencyScaleFactor(*ZZ->first().daughterPtr(0));
    (ZZ->first().daughterPtr(0))->setEfficenySFUnc(lepSF.second);
    theWeight*=lepSF.first;
    
    lepSF=leptonScaleFactors_.efficiencyScaleFactor(*ZZ->first().daughterPtr(1));
    (ZZ->first().daughterPtr(1))->setEfficenySFUnc(lepSF.second);
    theWeight*=lepSF.first;
    
    lepSF=leptonScaleFactors_.efficiencyScaleFactor(*ZZ->second().daughterPtr(0));
    ( ZZ->second().daughterPtr(0))->setEfficenySFUnc(lepSF.second);
    theWeight*=lepSF.first;
    
    lepSF=leptonScaleFactors_.efficiencyScaleFactor(*ZZ->second().daughterPtr(1));
    (ZZ->second().daughterPtr(1))->setEfficenySFUnc(lepSF.second);
    theWeight*=lepSF.first;
  }

  
  theHistograms.fill("weight_full"  , "All weights applied"                                    , 1200, -2, 10, theWeight);
  theHistograms.fill("weight_bare"  , "All weights, but efficiency and fake rate scale factors", 1200, -2, 10, theMCInfo.weight());
  theHistograms.fill("weight_pu"    , "Weight from PU reweighting procedure"                   , 1200, -2, 10, theMCInfo.puWeight());
  theHistograms.fill("weight_sample", "Weight from cross-section and luminosity"               , 1200, -2, 10, theMCInfo.sampleWeight());
  theHistograms.fill("weight_mcProc", "Weight from MC intrinsic event weight"                  , 1200, -2, 10, theMCInfo.mcProcWeight());
  theHistograms.fill("weight_efficiencySF", "Weight from data/MC lepton efficiency"            , 1200, -2, 10, ZZ->efficiencySF());
  theHistograms.fill("weight_fakeRateSF"  , "Weight from fake rate scale factor"               , 1200, -2, 10, ZZ->fakeRateSF());
  
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


  Long64_t nentries = /*10000; */theTree->GetEntries();  
  unweightedEventsInSR     = tree()->GetEntries("ZZCand.passSRZZOnShell_");
  unweightedEventsIn2P2FCR = tree()->GetEntries("ZZCand.passSelZLL_2P2F_ZZOnShell_");
  unweightedEventsIn3P1FCR = tree()->GetEntries("ZZCand.passSelZLL_3P1F_ZZOnShell_");
  begin();
  int looptimes=0;
  for (Long64_t jentry=0; jentry<nentries; ++jentry) {
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    if (!GetEntry(jentry)) continue;
    if (cut() < 0) continue;
    theCutCounter += theWeight;
    if(doBasicPlots_) fillBasicPlots();
    analyze();
    cout<<"\n loop #: "<< looptimes << "\n" << endl;
    looptimes++;
  }


  TFile fout(TString(outputfile),"RECREATE");
  fout.cd(); 

  end(fout);
  theHistograms.write(fout);

  fout.Close();
  cout<<"Events originally in input for the chosen region (" << Blue(regionType(region_)) << "): " << Green(theInputWeightedEvents)<< endl;
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
  theHistograms.fill(type+"_jetArea"      ,  60,  0, 1.2 , jet.jetArea      (), theWeight);
  theHistograms.fill(type+"_secvtxMass"   ,  50,  0, 5   , jet.secvtxMass   (), theWeight);
  theHistograms.fill(type+"_rawFactor"    ,  50,  0, 2.5 , jet.rawFactor    (), theWeight);

  theHistograms.fill(type+"_jecUncertainty", 50, 0, 0.1, jet.jecUncertainty(), theWeight);
}


