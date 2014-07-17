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



void GenEventAnalyzer::Init(TTree *tree)
{
  
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
  
  //  cout << "GenParticlesIn size = " << genParticlesIn->size() << endl;
 theWeight = lumi_*xsec_/200000.;

 if (genParticles->size() >= 3) {

  theHistograms.fill("Number of events", "Number of events", 10, 0, 10, 0, theWeight);

   std::vector<const Particle* > Genj;
   std::vector<const Particle* > Genq;
   std::vector<const Particle* > Genl;
   std::vector<const Particle* > GenZ;
   std::vector<const Particle* > GenW;
   
   
   foreach(const Particle& b, *genParticles) {
     
     int s_id = b.id();
     int id   = abs(s_id);
     
     if ( id < 7 || id == 21 ) Genj.push_back(&b);  //quark and gluons
     
     if ( id < 7 ) Genq.push_back(&b);               //quark
     
     else if ( id >= 11 && id <= 16 ) Genl.push_back(&b);        //leptons
     
     else if (id == 23)  GenZ.push_back(&b);                     // Z
     
     else if ( id == 24 ) GenW.push_back(&b);                    // W

     else cout << "!!!!!!! I've found id = " << id << endl;

     if ( id == 21 ) cout << "I'VE FOUND A GLUON!!!" << endl;  
      
   } 

    cout << "---------- genParticles: history information ----------" << endl; 
    cout << "Number of generated Z= "        << GenZ.size() << endl;
    cout << "Number of generated W= "        << GenW.size() << endl;
    cout << "Number of generated q = "       << Genq.size() << endl;
    cout << "Number of generated q and g = " << Genj.size() << endl;
    cout << "Number of generated l = "       << Genl.size() << endl;
    cout << "Weight = " << theWeight << endl;
    
 
    if ( GenZ.size() >= 2 && GenW.size() >= 1 ) {
	
      const Particle* Z0 = GenZ.at(0);
      const Particle* Z1 = GenZ.at(1);
      const Particle* W  = GenW.at(0);
      
      const Particle* l0 = Genl.at(0);
      const Particle* l1 = Genl.at(1);
      const Particle* l2 = Genl.at(2);
      const Particle* l3 = Genl.at(3);
      
      const Particle* j0 = Genq.at(0);
      const Particle* j1 = Genq.at(1);

      TLorentzVector p_Z0 = Z0->p4();
      TLorentzVector p_Z1 = Z1->p4();
      TLorentzVector p_W  = W->p4();
      
      TLorentzVector p_l0 = l0->p4();
      TLorentzVector p_l1 = l1->p4();
      TLorentzVector p_l2 = l2->p4();
      TLorentzVector p_l3 = l3->p4();
      
      TLorentzVector p_j0 = j0->p4();
      TLorentzVector p_j1 = j1->p4();
      
      TLorentzVector p_4l = p_l0 + p_l1 + p_l2 + p_l3;
      TLorentzVector p_jj = p_j0 + p_j1;
      TLorentzVector p_6f = p_4l + p_jj;

      double Detajj = p_j0.Eta() - p_j1.Eta();
      
      
      //------------Mass--------------

      theHistograms.fill("M_Z0", "M_Z0", 100, 40, 140, p_Z0.M(), theWeight);
      theHistograms.fill("M_Z1", "M_Z1", 100, 40, 140, p_Z1.M(), theWeight);
      theHistograms.fill("M_W" , "M_W" , 100, 40, 140, p_W.M() , theWeight);
      
    
      theHistograms.fill("M_ll0", "M_ll0", 100, 40, 140, (p_l0 + p_l1).M(), theWeight);
      theHistograms.fill("M_ll1", "M_ll1", 100, 40, 140, (p_l2 + p_l3).M(), theWeight);
      theHistograms.fill("M_jj" , "M_jj" , 100, 40, 140, (p_j0 + p_j1).M(), theWeight);
      theHistograms.fill("M_4l" , "M_4l" , 1000, 0, 1000, p_4l.M()        , theWeight);
      theHistograms.fill("M_6f" , "M_6f" , 3000, 0, 3000, p_6f.M()        , theWeight);
      
      
      
      //------------Pt--------------

      theHistograms.fill("Pt_Z0", "Pt_Z0", 500, 0, 500, Z0->pt(), theWeight);
      theHistograms.fill("Pt_Z1", "Pt_Z1", 500, 0, 500, Z1->pt(), theWeight);
      theHistograms.fill("Pt_W" , "Pt_W" , 500, 0, 500, W->pt() , theWeight);
      
      
      theHistograms.fill("Pt_l0", "Pt_l0", 1000, 0, 1000, l0->pt() , theWeight);
      theHistograms.fill("Pt_l1", "Pt_l1", 1000, 0, 1000, l1->pt() , theWeight);
      theHistograms.fill("Pt_l2", "Pt_l2", 1000, 0, 1000, l2->pt() , theWeight);
      theHistograms.fill("Pt_l3", "Pt_l3", 1000, 0, 1000, l3->pt() , theWeight);
      theHistograms.fill("Pt_j0", "Pt_j0", 1000, 0, 1000, j0->pt() , theWeight);
      theHistograms.fill("Pt_j1", "Pt_j1", 1000, 0, 1000, j1->pt() , theWeight);
      theHistograms.fill("Pt_4l", "Pt_4l", 1000, 0, 1000, p_4l.Pt(), theWeight);



      //------------Deta--------------
      
      theHistograms.fill("Deta_jj", "Deta_jj", 100, 0, 5, Detajj , theWeight);
      
    }
    
 }
 
 // else cout << "GenParticles size < 3. It is  = " << genParticles->size() <<endl;
 
}




// ------------------------------------------------------------------------------------------ //
// -------------------------------------- Preselection -------------------------------------- //
// ------------------------------------------------------------------------------------------ //


Int_t GenEventAnalyzer::cut() {
  
  bool pass = true;
  
  return pass ? 1 : -1;
}




