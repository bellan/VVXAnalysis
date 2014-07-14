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
				   double lumi 
				   )
 //  : theGenMCInfo(filename, lumi)
//   , theWeight(1.)
//   , theCutCounter(0.)
// theInputWeightedEvents(0.)
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



void GenEventAnalyzer::Init(TTree *tree)
{
  //  TH1F::SetDefaultSumw2(kTRUE);
  
  // Set branch addresses and branch pointers
  if (!tree) return;
  theTree = tree;
  fCurrent = -1;
 
  // Gen Particles   
  genParticles     = 0; 
  b_genParticles   = 0; 
  theTree->SetBranchAddress("genParticles"  , &genParticles  , &b_genParticles);
  genParticlesIn   = 0; 
  b_genParticlesIn = 0; 
  theTree->SetBranchAddress("genParticlesIn", &genParticlesIn, &b_genParticlesIn);
 
}



// ------------------------------------------------------------------------------------------ //
// ------------------------------------- Master the loop ------------------------------------ //
// ------------------------------------------------------------------------------------------ //



Int_t GenEventAnalyzer::GetEntry(Long64_t entry){
  // Read contents of entry.
      
  if (!theTree) return 0;
  
  int e =  theTree->GetEntry(entry);


 //  theWeight = lumi*XS/;
//   theInputWeightedEvents += theWeight;

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
 //    theCutCounter += theWeight;

    //  theHistograms.fill("weight",100, 0, 200, theWeight);

    analyze();
  }

  TFile fout(TString(outputfile),"RECREATE");
  fout.cd(); 
  theHistograms.write(fout);

  end(fout);
  fout.Close();
  //cout<<"Events in input: " << Green(theInputWeightedEvents)<< endl;
  // cout<<"Events passing all cuts: "<< Green(theCutCounter) << endl;
}


// ------------------------------------------------------------------------------------------ //
// --------------------------------------- Analyze ------------------------------------------ //
// ------------------------------------------------------------------------------------------ //


void GenEventAnalyzer::analyze() {
  
 if (genParticles->size() >= 3) {

  theHistograms.fill("Number of events", "Number of events", 10, 0, 10, 0);//, theWeight);

   std::vector<const Particle* > Genj;
   std::vector<const Particle* > Genq;
   std::vector<const Particle* > Genl;
   std::vector<const Particle* > GenZ;
   std::vector<const Particle* > GenW;
   
   
   foreach(const Particle& b, *genParticles) {
     
     int s_id = b.id();
     int id   = abs(s_id);
     
     if ( (id < 7 || id == 21) && id > 5 ) Genj.push_back(&b);  //quark and gluons
     
     if ( id < 7  && id > 5 ) Genq.push_back(&b);               //quark
     
     else if ( id >= 11 && id <= 16 ) Genl.push_back(&b);        //leptons
     
     else if (id == 23)  GenZ.push_back(&b);                     // Z
     
     else if ( id == 24 ) GenW.push_back(&b);                    // W

     else cout << "I've found id = " << id << endl;
      
   } 



    cout << "---------- genParticles: history information ----------" << endl; 
    cout << "Number of generated Z= "        << GenZ.size() << endl;
    cout << "Number of generated W= "        << GenW.size() << endl;
    cout << "Number of generated q= "        << Genq.size() << endl;
    cout << "Number of generated q and g= "  << Genj.size() << endl;
    cout << "Number of generated l= "        << Genl.size() << endl;
    
 
      if ( GenZ.size() >= 2 && GenW.size() >= 1 ) {
	
	const Particle* Z0gen = GenZ.at(0);
	const Particle* Z1gen = GenZ.at(1);
	const Particle* Wgen  = GenW.at(0);
	
	
	//------------Mass--------------
	
	theHistograms.fill("Z0Gen_Mass", "Z0Gen_Mass", 200, 0, 200, Z0gen->p4().M());
	theHistograms.fill("Z1Gen_Mass", "Z1Gen_Mass", 200, 0, 200, Z1gen->p4().M());
	theHistograms.fill("WGen_Mass" , "WGen_Mass" , 200, 0, 200, Wgen->p4().M());

	theHistograms.fill("Z0Gen_Mass", "Z0Gen_Mass", 200, 0, 200, Z0gen->p4().M());
	theHistograms.fill("Z1Gen_Mass", "Z1Gen_Mass", 200, 0, 200, Z1gen->p4().M());
	theHistograms.fill("WGen_Mass" , "WGen_Mass" , 200, 0, 200, Wgen->p4().M());
	
	//------------Pt--------------
	
	theHistograms.fill("Z0Gen_Pt"  , "Z0Gen_Pt"  , 300, 0, 300, Z0gen->pt()    );
	theHistograms.fill("Z1Gen_Pt"  , "Z1Gen_Pt"  , 300, 0, 300, Z1gen->pt()    );
	theHistograms.fill("WGen_Pt"   , "WGen_Pt"   , 300, 0, 300, Wgen->pt()     );
	
      }
      
 }
 
}




// ------------------------------------------------------------------------------------------ //
// -------------------------------------- Preselection -------------------------------------- //
// ------------------------------------------------------------------------------------------ //


Int_t GenEventAnalyzer::cut() {
  
  bool pass = true;
  
  return pass ? 1 : -1;
}




