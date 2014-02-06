#include "EventAnalyzer.h"

#include <TChain.h>
#include <TFile.h>
#include <TH1F.h>
#include <TH1I.h>
#include <TString.h>

using std::cout;
using std::endl;
using std::flush;

#include <boost/foreach.hpp>
#define foreach BOOST_FOREACH

#include <boost/lexical_cast.hpp>

EventAnalyzer::EventAnalyzer(std::string filename, 
			     double lumi, 
			     double externalXSection)
  : theWeight(1.)
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
  // Set the environment to read one entry
  if (!theTree) return -5;
  Long64_t centry = theTree->LoadTree(entry);
  if (centry < 0) return centry;
  if (!theTree->InheritsFrom(TChain::Class()))  return centry;
  TChain *chain = (TChain*)theTree;
  if (chain->GetTreeNumber() != fCurrent) {
    fCurrent = chain->GetTreeNumber();
    Notify();
   }
   return centry;
}

void EventAnalyzer::Init(TTree *tree)
{
  TH1F::SetDefaultSumw2(kTRUE);

   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set object pointer
   
   // Set branch addresses and branch pointers
   if (!tree) return;
   theTree = tree;
   fCurrent = -1;
   //theTree->SetMakeClass(1); // Turn it on only if there are not omplex object in the tree
   
   // Muons   
   muons     = 0; b_muons     = 0; theTree->SetBranchAddress("muons"    , &muons    , &b_muons    );

   // Electrons   
   electrons = 0; b_electrons = 0; theTree->SetBranchAddress("electrons", &electrons, &b_electrons);

   // Jets   
   jets = 0;      b_jets = 0;      theTree->SetBranchAddress("jets", &jets, &b_jets);

   // MET
   met = new phys::Particle();
   b_met         = 0; theTree->SetBranchAddress("met"    , &met    ,  &b_met    );
   
   b_nvtx        = 0; theTree->SetBranchAddress("nvtxs"  , &nvtx   , &b_nvtx   );
   b_rho         = 0; theTree->SetBranchAddress("rho"    , &rho    , &b_rho    );
   b_weight      = 0; theTree->SetBranchAddress("weight" , &weight , &b_weight );
   cout << "Weight from the event sample type: " << theSampleWeight << ", total weight (including PU reweight, if applicable): " << theWeight << endl;
   if(theSampleWeight != 1) b_puweight    = 0; theTree->SetBranchAddress("puweight" , &puweight , &b_puweight );
   b_xsec        = 0; theTree->SetBranchAddress("xsec"   , &xsec   , &b_xsec  );
   b_totalEvents = 0; theTree->SetBranchAddress("totalEvents", &totalEvents  , &b_totalEvents  );


   thePlots["nelectrons"] = new TH1I("nelectrons","Number of electrons",10,0,10);
   thePlots["nmuons"]     = new TH1I("nmuons"    ,"Number of muons"    ,10,0,10);  
   
   thePlots["nvtx"]       = new TH1I("nvtx"      ,"Number of vertices" ,100,0,100);  
   thePlots["rho"]        = new TH1F("rho"      ,"Mean energy density" ,100,0,50);  
   thePlots["met"]        = new TH1F("met"      ,"Missing energy"      ,200,0,800);  

   book("mu");
   book("e");
   bookExtraForElectrons("e");

   begin();
   book();

   Notify();
}

Bool_t EventAnalyzer::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void EventAnalyzer::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!theTree) return;
   theTree->Show(entry);
}

// FIXME! Split
void EventAnalyzer::book(const std::string &type){
   thePlots[type+"_pt"]    = new TH1F(TString(type+"_pt") , "p_{T} spectrum", 100,   0   , 500   );  
   thePlots[type+"_eta"]   = new TH1F(TString(type+"_eta"), "#eta spectrum" , 100,  -2.5 ,   2.5 );  
   thePlots[type+"_phi"]   = new TH1F(TString(type+"_phi"), "#phi spectrum" ,  50,  -3.15,   3.15);   
   thePlots[type+"_dxy"]   = new TH1F(TString(type+"_dxy"), "d_{xy}"        , 200,   0   ,   0.5 );               
   thePlots[type+"_dz"]    = new TH1F(TString(type+"_dz") , "d_{z}"         , 200,   0   ,   0.5 );                
   thePlots[type+"_sip"]   = new TH1F(TString(type+"_sip"), "sip"           , 100, -10   ,  10   ); 
   thePlots[type+"_combRelIso"]      = new TH1F(TString(type+"_combRelIso")      , "combRelIso"     , 200,   0,   1); 
   thePlots[type+"_pfCombRelIso"]    = new TH1F(TString(type+"_pfCombRelIso")    , "pfCombRelIso"   , 200,   0,   1);   
   thePlots[type+"_pfNeutralHadIso"] = new TH1F(TString(type+"_pfNeutralHadIso") , "pfNeutralHadIso", 200,   0,   1);   
   thePlots[type+"_pfChargedHadIso"] = new TH1F(TString(type+"_pfChargedHadIso") , "pfChargedHadIso", 200,   0,   1);  
   thePlots[type+"_pfPhotonHadIso"]  = new TH1F(TString(type+"_pfPhotonHadIso")  , "pfPhotonHadIso" , 200,   0,   1);    
   thePlots[type+"_rho"]             = new TH1F(TString(type+"_rho")             , "rho"            , 200,   0, 200);   
}


void EventAnalyzer::bookExtraForElectrons(const std::string &type){
  thePlots[type+"_energy"]     = new TH1F(TString(type+"_energy")    , "energy spectrum"   , 100,  0   , 500   );  
  thePlots[type+"_phiWidth"]   = new TH1F(TString(type+"_phiWidth")  , "#phi width"        ,  50, -3.15,   3.15);   
  thePlots[type+"_etaWidth"]   = new TH1F(TString(type+"_etaWidth")  , "#eta width"        , 100, -2.5 ,   2.5 );  
  thePlots[type+"_BDT"]        = new TH1F(TString(type+"_BDT")       , "BDT"               , 100, -1   ,   1   );   
  thePlots[type+"_isBDT"]      = new TH1I(TString(type+"_isBDT")     , "is BDT?"           ,   2,  0   ,   2   );  
  thePlots[type+"_missingHit"] = new TH1I(TString(type+"_missingHit"), "Missing hits"      ,  50,  0   ,  50   );   
  thePlots[type+"_nCrystals"]  = new TH1I(TString(type+"_nCrystals") , "Number of Crystals",  50,  0   ,  50   );   
}


void EventAnalyzer::bookExtraForJets(const std::string &type){

    thePlots[type+"_nConstituents"] = new TH1I(TString(type+"_nConstituents"), "nConstituents", 100, 0, 100); 
    thePlots[type+"_nCharged"]      = new TH1I(TString(type+"_nCharged")     , "nCharged"     , 100, 0, 100);      
    thePlots[type+"_nNeutral"]      = new TH1I(TString(type+"_nNeutral")     , "nNeutral"     , 100, 0, 100);      

    thePlots[type+"_neutralHadronEnergyFraction"] = new TH1F(TString(type+"_neutralHadronEnergyFraction"), "neutralHadronEnergyFraction", 100, 0, 1);   
    thePlots[type+"_chargedHadronEnergyFraction"] = new TH1F(TString(type+"_chargedHadronEnergyFraction"), "chargedHadronEnergyFraction", 100, 0, 1);
    thePlots[type+"_chargedEmEnergyFraction"]     = new TH1F(TString(type+"_chargedEmEnergyFraction"    ), "chargedEmEnergyFraction"    , 100, 0, 1);    
    thePlots[type+"_neutralEmEnergyFraction"]     = new TH1F(TString(type+"_neutralEmEnergyFraction"    ), "neutralEmEnergyFraction"    , 100, 0, 1);    
    thePlots[type+"_muonEnergyFraction"]          = new TH1F(TString(type+"_muonEnergyFraction"         ), "muonEnergyFraction"         , 100, 0, 1); 

    thePlots[type+"_csvtagger"]     = new TH1F(TString(type+"_csvtagger"    ), "csvtagger"    , 200, -1, 1);             
    thePlots[type+"_girth"]         = new TH1F(TString(type+"_girth"        ), "girth"        , 200,  0, 1);    
    thePlots[type+"_girth_charged"] = new TH1F(TString(type+"_girth_charged"), "girth_charged", 200,  0, 1);
    thePlots[type+"_ptd"]           = new TH1F(TString(type+"_ptd"          ), "ptd"          , 1, 1, 1);          
    thePlots[type+"_rms"]           = new TH1F(TString(type+"_rms"          ), "rms"          , 1, 1, 1);          
    thePlots[type+"_beta"]          = new TH1F(TString(type+"_beta"         ), "beta"         , 1, 1, 1);         
    thePlots[type+"_jetArea"]       = new TH1F(TString(type+"_jetArea"      ), "jetArea"      , 1, 1, 1);      
    thePlots[type+"_secvtxMass"]    = new TH1F(TString(type+"_secvtxMass"   ), "secvtxMass"   , 1, 1, 1);   
    thePlots[type+"_Lxy"]           = new TH1F(TString(type+"_Lxy"          ), "Lxy"          , 1, 1, 1);          	 
    thePlots[type+"_LxyErr"]        = new TH1F(TString(type+"_LxyErr"       ), "LxyErr"       , 1, 1, 1);        
    thePlots[type+"_rawFactor"]     = new TH1F(TString(type+"_rawFactor"    ), "rawFactor"    , 1, 1, 1);    

    thePlots[type+"_uncOnFourVectorScale"] = new TH1F(TString(type+"_uncOnFourVectorScale"), "uncOnFourVectorScale", 1, 1, 1);
    thePlots[type+"_puMVAFull"]   = new TH1F(TString(type+"_puMVAFull"  ), "puMVAFull"  , 1, 1, 1);   
    thePlots[type+"_puMVASimple"] = new TH1F(TString(type+"_puMVASimple"), "puMVASimple", 1, 1, 1); 
    thePlots[type+"_puCutBased"]  = new TH1F(TString(type+"_puCutBased" ), "puCutBased" , 1, 1, 1);  

    thePlots[type+"_pass_puMVAFull_loose"]    = new TH1I(TString(type+"_pass_puMVAFull_loose"   ), "pass_puMVAFull_loose"   , 2, 0, 2);   
    thePlots[type+"_pass_pUMVAFull_medium"]   = new TH1I(TString(type+"_pass_pUMVAFull_medium"  ), "pass_pUMVAFull_medium"  , 2, 0, 2);  
    thePlots[type+"_pass_pUMVAFull_tight"]    = new TH1I(TString(type+"_pass_pUMVAFull_tight"   ), "pass_pUMVAFull_tight"   , 2, 0, 2);   
									                                                    
    thePlots[type+"_pass_puMVASimple_loose"]  = new TH1I(TString(type+"_pass_puMVASimple_loose" ), "pass_puMVASimple_loose" , 2, 0, 2); 
    thePlots[type+"_pass_puMVASimple_medium"] = new TH1I(TString(type+"_pass_puMVASimple_medium"), "pass_puMVASimple_medium", 2, 0, 2);
    thePlots[type+"_pass_puMVASimple_tight"]  = new TH1I(TString(type+"_pass_puMVASimple_tight" ), "pass_puMVASimple_tight" , 2, 0, 2); 
									                                                    
    thePlots[type+"_pass_puCutBased_loose"]   = new TH1I(TString(type+"_pass_puCutBased_loose"  ), "pass_puCutBased_loose"  , 2, 0, 2);  
    thePlots[type+"_pass_puCutBased_medium"]  = new TH1I(TString(type+"_pass_puCutBased_medium" ), "pass_puCutBased_medium" , 2, 0, 2); 
    thePlots[type+"_pass_puCutBased_tight"]   = new TH1I(TString(type+"_pass_puCutBased_tight"  ), "pass_puCutBased_tight"  , 2, 0, 2);  
}

Int_t EventAnalyzer::cut(){
  
  bool pass = true;

  if(pass) ++theCutCounter;

  return pass ? 1 : -1;
}


void EventAnalyzer::loop(const std::string outputfile){
  if (theTree == 0) return;

  Long64_t nentries = theTree->GetEntries();  
  Long64_t nbytes = 0, nb = 0;

  for (Long64_t jentry=0; jentry<nentries;jentry++) {

    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = GetEntry(jentry);  nbytes += nb; 
    theWeight = theSampleWeight;
    if(theWeight != 1) theWeight *= puweight;

    if (cut() < 0) continue;

    basicPlots();
    analyze();

    histos.fill<TH1F>("pippo",100,0,100,met->pt());

  }
  TFile fout(TString(outputfile),"RECREATE");
  fout.cd(); 
  for(std::map<std::string,TH1*>::const_iterator h = thePlots.begin(); h != thePlots.end(); ++h)
    h->second->Write();
  histos.write(fout);

  end(fout);
  fout.Close();
  cout<<"Event passing all cuts: " << theCutCounter << endl;
  
}

void EventAnalyzer::basicPlots(){
  thePlots["nvtx"]->Fill(nvtx,theWeight);
  thePlots["rho"]->Fill(rho,theWeight);
  thePlots["met"]->Fill(met->pt(),theWeight);
  
  thePlots["nmuons"]->Fill(muons->size(),theWeight);
  foreach(const phys::Lepton& mu, *muons) fillLeptonPlots("mu", mu);
  
  thePlots["nelectrons"]->Fill(electrons->size(),theWeight);
  foreach(const phys::Electron& e, *electrons) fillLeptonPlots("e",e);
}

void EventAnalyzer::fillParticlePlots(const std::string &type, const phys::Particle & lepton){
  thePlots[type+"_pt"]->Fill(lepton.pt(),   theWeight);
  thePlots[type+"_eta"]->Fill(lepton.eta(), theWeight);
  thePlots[type+"_phi"]->Fill(lepton.phi(), theWeight);
}

void EventAnalyzer::fillLeptonPlots(const std::string &type, const phys::Lepton & lepton){
  
  fillParticlePlots(type, lepton);

  thePlots[type+"_dxy"]->Fill(lepton.dxy(), theWeight);   
  thePlots[type+"_dz"] ->Fill(lepton.dz(), theWeight);
  thePlots[type+"_sip"]->Fill(lepton.sip(), theWeight); 

  thePlots[type+"_combRelIso"]     ->Fill(lepton.combRelIso()     , theWeight);       
  thePlots[type+"_pfCombRelIso"]   ->Fill(lepton.pfChargedHadIso(), theWeight);           
  thePlots[type+"_pfNeutralHadIso"]->Fill(lepton.pfNeutralHadIso(), theWeight);        
  thePlots[type+"_pfChargedHadIso"]->Fill(lepton.pfPhotonIso()    , theWeight);        
  thePlots[type+"_pfPhotonHadIso"] ->Fill(lepton.pfCombRelIso()   , theWeight);         
  thePlots[type+"_rho"]            ->Fill(lepton.rho()            , theWeight);                    
}

void EventAnalyzer::fillElectronPlots (const std::string &type, const phys::Electron &electron){
  fillLeptonPlots  (type, electron);
  fillElectronPlots(type, electron);
}

void EventAnalyzer::fillExtraPlotsForElectrons(const std::string &type, const phys::Electron &electron){

  thePlots[type+"_energy"]    ->Fill(electron.energy()    , theWeight);  
  thePlots[type+"_phiWidth"]  ->Fill(electron.phiWidth()  , theWeight);  
  thePlots[type+"_etaWidth"]  ->Fill(electron.etaWidth()  , theWeight);  
  thePlots[type+"_BDT"]       ->Fill(electron.BDT()       , theWeight);  
  thePlots[type+"_isBDT"]     ->Fill(electron.isBDT()     , theWeight);  
  thePlots[type+"_missingHit"]->Fill(electron.missingHit(), theWeight);  
  thePlots[type+"_nCrystals"] ->Fill(electron.nCrystals() , theWeight);  
}

void EventAnalyzer::fillJetPlots(const std::string &type, const phys::Jet      &jet){}
//   fillParticlePlots(type, jet);
  

//     nConstituents() 
//     nCharged()      
//     nNeutral()      

//     neutralHadronEnergyFraction()   
//     chargedHadronEnergyFraction()
//     chargedEmEnergyFraction()    
//     neutralEmEnergyFraction()    
//     muonEnergyFraction()         

//     csvtagger()             
//     girth()        
//     girth_charged()
//     ptd()          
//     rms()          
//     beta()         
//     jetArea()      
//     secvtxMass()   
//     Lxy()          	 
//     LxyErr()        
//     rawFactor()    
//     uncOnFourVectorScale()
//     puMVAFull()   
//     puMVASimple() 
//     puCutBased()  
					 				      
//     pass_puMVAFull_loose()   
//     pass_pUMVAFull_medium()  
//     pass_pUMVAFull_tight()   
    			     
//     pass_puMVASimple_loose() 
//     pass_puMVASimple_medium()
//     pass_puMVASimple_tight() 
    			     
//     pass_puCutBased_loose()  
//     pass_puCutBased_medium() 
//     pass_puCutBased_tight()  
//}
