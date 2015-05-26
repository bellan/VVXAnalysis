/** \class GenEventAnalyzer
 *  Base class for MadGraph event analyzers. Analyzers have to inherit from this class 
 *  and implement the pure virtual function analyze(), called each event.
 *
 */

#include "VVXAnalysis/TreeAnalysis/interface/GenEventAnalyzer.h"
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

{

  TChain *tree = new TChain("Cypress"); 
  tree->Add(filename.c_str());

  if (tree == 0) std::cout<<Important("Error in GenEventAnalyzer ctor:")<<" The tree has a null pointer."<<std::endl;
  Init(tree);
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
  

  


}



// ------------------------------------------------------------------------------------------ //
// ------------------------------------- Master the loop ------------------------------------ //
// ------------------------------------------------------------------------------------------ //



Int_t GenEventAnalyzer::GetEntry(Long64_t entry){
  // Read contents of entry.
      
  if (!theTree) return 0;
  
  int e =  theTree->GetEntry(entry);

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


void GenEventAnalyzer::analyze() {
  
  
  foreach(const Particle& gen, *genParticles) {
    cout << gen << endl;
  }    

}




// ------------------------------------------------------------------------------------------ //
// -------------------------------------- Preselection -------------------------------------- //
// ------------------------------------------------------------------------------------------ //


Int_t GenEventAnalyzer::cut() {
  
  bool pass = true;
  
  return pass ? 1 : -1;
}




