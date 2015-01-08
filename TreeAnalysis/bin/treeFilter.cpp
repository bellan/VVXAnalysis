#include "VVXAnalysis/Commons/interface/Colours.h"
#include "VVXAnalysis/DataFormats/interface/Particle.h"

#include <TChain.h>
#include <TBranch.h>
#include <TFile.h>
#include <iostream>

using std::cout;
using std::endl;

using namespace colour;

#include <boost/foreach.hpp>
#define foreach BOOST_FOREACH

int main (int argc, char ** argv){

  cout<<"Tree Filter."<<endl;
  
  //Get old file, old tree and set top branch address
  TChain *oldTree = new TChain("treePlanter/ElderTree");
  oldTree->Add("samples/WZZJets.root");

  Long64_t nentries = oldTree->GetEntries();

  cout << nentries << endl;
  
  std::vector<phys::Particle> *genParticles = 0;   TBranch *b_genParticles = 0;
  oldTree->SetBranchAddress("genParticles"  , &genParticles  , &b_genParticles);

   //Create a new file + a clone of old tree in new file
   TFile *newFile = new TFile("small.root","recreate");
   TTree *newTree = oldTree->CloneTree(0);

   for (Long64_t i=0;i<nentries; i++) {
      oldTree->GetEntry(i);

      cout << "-------------------------------" << endl;
      foreach(const phys::Particle &genParticle, *genParticles)
	cout << genParticle.id()<<endl;

      if (true) newTree->Fill();
   }
   //newTree->Print();
   newTree->AutoSave();
   
   // Copy ancillary file too.


   delete newFile;
}
