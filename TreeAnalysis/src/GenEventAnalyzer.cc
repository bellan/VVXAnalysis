/** \class GenEventAnalyzer
 *  Base class for MadGraph event analyzers. Analyzers have to inherit from this class 
 *  and implement the pure virtual function analyze(), called each event.
 *
 */

#include "VVXAnalysis/TreeAnalysis/interface/GenEventAnalyzer.h"
#include "VVXAnalysis/Commons/interface/Colours.h"
#include "VVXAnalysis/Commons/interface/Comparators.h"
#include "VVXAnalysis/Commons/interface/Constants.h"
#include "VVXAnalysis/Commons/interface/Utilities.h"
#include "VVXAnalysis/DataFormats/interface/TypeDefs.h"
#include "VVXAnalysis/DataFormats/interface/DiBoson.h"

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
using namespace phys;
using namespace std;




// ------------------------------------------------------------------------------------------ //
// ------------------------------- C'tor/Des'tor/Init the tree ------------------------------ //
// ------------------------------------------------------------------------------------------ //



GenEventAnalyzer::GenEventAnalyzer(std::string filename, 
				   double lumi,
				   float xsec
				   )
  : lumi_(lumi)
  , xsec_(xsec)
  , mcprocweight_(0)
  , passSignal(0)
  , passSignalTightFiducialRegion(0)
  , passHiggsFiducialRegion(0)

{

  TChain *tree = new TChain("Cypress"); 
  tree->Add(filename.c_str());

  if (tree == 0) std::cout<<Important("Error in GenEventAnalyzer ctor:")<<" The tree has a null pointer."<<std::endl;
  Init(tree);

  TFile *in = new TFile(filename.c_str());
  TH1F *hNInputEvents = (TH1F*)in->Get("NumberOfEvents");
  numberOfAnalyzedEvents_ = hNInputEvents->GetBinContent(1);
  numberOfInputEvents_    = hNInputEvents->GetBinContent(2);
  theWeight = lumi_*xsec_/numberOfAnalyzedEvents_;
  cout << "Standard weight: " << theWeight << endl;
  

  TH1F *hSumProcWeight = (TH1F*)in->Get("summcprocweights");
  double summcprocweights = hSumProcWeight->GetBinContent(1);
  cout<< "summcprocweights: " << summcprocweights <<endl;
  theWeight = lumi_*xsec_/summcprocweights;


  cout << "Number of analyzed events: "          << Green(numberOfAnalyzedEvents_) << endl
       << "Number of selected events in input: " << Green(numberOfInputEvents_)    << endl
       << "Cross section: " << Green(xsec_) << Green(" pb") << " luminosity: " << Green(lumi_) << Green("/pb") << " sample weight: " << Green(theWeight) << endl;

  delete hNInputEvents;
  delete in;
}



GenEventAnalyzer::~GenEventAnalyzer(){
  if (!theTree) return;
  delete theTree->GetCurrentFile();
}



void GenEventAnalyzer::Init(TTree *tree){

  
  // Set branch addresses and branch pointers
  if (!tree) return;
  theTree = tree;
  fCurrent = -1;
 
  // Gen Particles   

  //genParticlesIn   = 0; b_genParticlesIn = 0; theTree->SetBranchAddress("genParticlesIn", &genParticlesIn, &b_genParticlesIn);
  genParticles   = 0;                                                b_genParticles   = 0; theTree->SetBranchAddress("genParticles"  , &genParticles  , &b_genParticles);
  genVBParticles = new std::vector<phys::Boson<phys::Particle> > (); b_genVBParticles = 0; theTree->SetBranchAddress("genVBParticles", &genVBParticles, &b_genVBParticles);
  pgenJets       = 0;                                                b_pgenJets       = 0; theTree->SetBranchAddress("genJets"       , &pgenJets      , &b_pgenJets);

  genJets  = new std::vector<phys::Particle>(); centralGenJets  = new std::vector<phys::Particle>();

  b_mcprocweight = 0; theTree->SetBranchAddress("mcprocweight"        , &mcprocweight_        , &b_mcprocweight        );  
}



// ------------------------------------------------------------------------------------------ //
// ------------------------------------- Master the loop ------------------------------------ //
// ------------------------------------------------------------------------------------------ //



Int_t GenEventAnalyzer::GetEntry(Long64_t entry){
  // Read contents of entry.
      
  if (!theTree) return 0;
  
  int e =  theTree->GetEntry(entry);

  stable_sort(pgenJets->begin(),  pgenJets->end(),  phys::PtComparator());
  genJets->clear(); centralGenJets->clear();
  foreach(const phys::Particle &jet, *pgenJets)
    if(jet.pt() > 30){
      bool leptonMatch = false;
      foreach(const phys::Particle &gen, *genParticles)
	if(physmath::deltaR(gen,jet) < 0.4 && (abs(gen.id()) == 11 || abs(gen.id()) == 13)) leptonMatch = true;
      
      if(!leptonMatch){
	if(fabs(jet.eta()) < 4.7) genJets->push_back(jet);
	if(fabs(jet.eta()) < 2.4) centralGenJets->push_back(jet);
      }
    }


  return e;
}



Long64_t GenEventAnalyzer::LoadTree(Long64_t entry){

  if (!theTree) return -5;
  Long64_t centry = theTree->LoadTree(entry);

  if (centry < 0) return centry;
  if (!theTree->InheritsFrom(TChain::Class()))  return centry;

  TChain *chain = (TChain*)theTree;
  if (chain->GetTreeNumber() != fCurrent) fCurrent = chain->GetTreeNumber();
  
  return centry;
}


void GenEventAnalyzer::loop(const std::string outputfile){

  if (theTree == 0) return;

  Long64_t nentries = theTree->GetEntries();  
  Long64_t nbytes = 0, nb = 0;

  begin();

  for (Long64_t jentry=0; jentry<nentries; ++jentry) {
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;

    nb = GetEntry(jentry);  nbytes += nb; 

    cout <<"New weight: " << weight() << ", mcprocweight: " << mcprocweight_ << endl;

    if (cut() < 0) continue;
    analyze();
  }

  TFile fout(TString(outputfile),"RECREATE");
  fout.cd(); 
  theHistograms.write(fout);

  end(fout);
  fout.Close();
 
}


// ------------------------------------------------------------------------------------------ //
// --------------------------------------- Analyze ------------------------------------------ //
// ------------------------------------------------------------------------------------------ //



void GenEventAnalyzer::makeBasicPlots(const std::string &selection, const zz::SignalTopology& zzSignalTopology){

  phys::BosonParticle Z0 = std::get<1>(zzSignalTopology);
  phys::BosonParticle Z1 = std::get<2>(zzSignalTopology);

  if(!Z0.isValid() or !Z1.isValid()) return;

  phys::DiBoson<phys::Particle,phys::Particle> ZZ(Z0,Z1);
  
  theHistograms.fill(selection+"_Z0mass"    , selection+" Z0 mass", 100,  0,150,  Z0.mass(), weight());
  theHistograms.fill(selection+"_Z1mass"    , selection+" Z1 mass", 100,  0,150,  Z1.mass(), weight());
  theHistograms.fill(selection+"_Z0pt"      , selection+" Z0 pT"  , 30,   0,300,  Z0.pt(), weight());
  theHistograms.fill(selection+"_Z1pt"      , selection+" Z1 pT"  , 30,   0,300,  Z1.pt(), weight());
  theHistograms.fill(selection+"_ZZmass"    , selection+" ZZ mass", 75, 100,1000, ZZ.mass(), weight());
  theHistograms.fill(selection+"_ZZmassZoom", selection+" ZZ mass", 10, 100,140,  ZZ.mass(), weight());
  theHistograms.fill(selection+"_ZZpt"      , selection+" ZZ pT"  , 50,   0,500,  ZZ.pt(), weight());
  
  theHistograms.fill(selection+"_DeltaPhiZ0Z1", selection+" #Delta #Phi", 30,  0, M_PI, fabs(physmath::deltaPhi(Z0.phi(),Z1.phi())),weight());
  theHistograms.fill(selection+"_DeltaRZ0Z1"  , selection+" #Delta R"   , 100, 0,10   , fabs(physmath::deltaR(Z0,Z1)), weight());

  std::vector<phys::Particle>  leptons;
  for(int i = 0; i < 2; ++i) {leptons.push_back(Z0.daughter(i));leptons.push_back(Z1.daughter(i));}
  stable_sort(leptons.begin(),  leptons.end(),  phys::PtComparator());
  theHistograms.fill(selection+"_leadingLeptonPt", selection+"_Leading lepton pt", 30,0,300, leptons.at(0).pt(), weight());
  

  theHistograms.fill(selection+"_nJets"       , selection+" number of jets", 10, 0, 10, genJets->size(), weight());
  theHistograms.fill(selection+"_nCentralJets", selection+" number of jets", 10, 0, 10, centralGenJets->size(), weight());
  if(genJets->size() > 1)        theHistograms.fill(selection+"_DeltaEtaJJ", selection+" #Delta #eta(j,j)", 25,0,7.5, abs(genJets->at(0).eta()-genJets->at(1).eta()), weight());
  if(centralGenJets->size() > 1) theHistograms.fill(selection+"_mJJ"       , selection+" m_{jj}"          , 20,0,400, (centralGenJets->at(0).p4()+centralGenJets->at(1).p4()).M(), weight());


}


void GenEventAnalyzer::analyze() {
  
  
  zz::SignalTopology zzSignalTopology = zz::getSignalTopology(*genParticles, *genJets, *genJets); // FIXME AK8 as third argument!
  
                                                    makeBasicPlots("All",zzSignalTopology);
  if(std::get<0>(zzSignalTopology) > 0)            {makeBasicPlots("Signal",zzSignalTopology);                    ++passSignal;}
  if(zz::inTightFiducialRegion(zzSignalTopology))  {makeBasicPlots("SignalTightFiducialRegion",zzSignalTopology); ++passSignalTightFiducialRegion;}
  if(std::get<0>(zzSignalTopology) == 0)            makeBasicPlots("NonSignal",zzSignalTopology);
  if(zz::inHiggsFiducialRegion(zzSignalTopology))  {makeBasicPlots("HiggsFiducialRegion",zzSignalTopology);       ++passHiggsFiducialRegion;}


  if(zz::inHiggsFiducialRegion(zzSignalTopology) > 0){
    phys::BosonParticle Z1 = std::get<2>(zzSignalTopology);
    if (Z1.mass() < 1)
      cout << "----------------- "<< Z1.mass()<<" ---------------------------" << endl
	   << Z1.daughter(0) << " status: " << Z1.daughter(0).genStatusFlags().test(GenStatusBit::fromHardProcess) << endl
	   << Z1.daughter(1) << " status: " << Z1.daughter(1).genStatusFlags().test(GenStatusBit::fromHardProcess) << endl
	   << endl;
  }
}




// ------------------------------------------------------------------------------------------ //
// -------------------------------------- Preselection -------------------------------------- //
// ------------------------------------------------------------------------------------------ //


Int_t GenEventAnalyzer::cut() {
  
  bool pass = true;
  
  return pass ? 1 : -1;
}



void GenEventAnalyzer::end(TFile &){
  cout << "Signal cross section: "                       << passSignal/float(numberOfAnalyzedEvents_)*xsec_                    << " pb" << endl
       << "Signal Tight Fiducial Region cross section: " << passSignalTightFiducialRegion/float(numberOfAnalyzedEvents_)*xsec_ << " pb" << endl
       << "Higgs Fiducial Region cross section: "        << passHiggsFiducialRegion/float(numberOfAnalyzedEvents_)*xsec_       << " pb" << endl;

}
