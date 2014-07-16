#include "VVXAnalysis/TreeAnalysis/interface/MCInfo.h"
#include "VVXAnalysis/Commons/interface/Colours.h"

#include "TChain.h"
#include <iostream>

using namespace colour;

MCInfo::MCInfo(const std::string& filename, const double & lumi, const double& externalXSection)
  : luminosity_(lumi)
  , internalCrossSection_(-1)
  , externalCrossSection_(externalXSection)
  , crossSection_(&externalCrossSection_)
  , genEvents_(-1)
  , analyzedEvents_(-1)
  , sampleWeight_(1)
  , mcprocweight_(1)
  , puweight_(1)
  , summcprocweight_(0)
  , sumpuweight_(0)
  , sumpumcprocweight_(0)
  , preSkimCounter_(0) 
  , postSkimCounter_(0)

{

  if(lumi <= 0) return;
  
  TChain *tree = new TChain("treePlanter/HollyTree");
  tree->Add(filename.c_str());

  if (tree == 0) return;
  
  TBranch *b_genEvents            = 0;
  TBranch *b_analyzedEvents       = 0;
  TBranch *b_internalCrossSection = 0;  
  TBranch *b_externalCrossSection = 0;  
  TBranch *b_summcprocweight    = 0;  
  TBranch *b_sumpuweight        = 0;  
  TBranch *b_sumpumcprocweight  = 0;
  TBranch *b_preSkimCounter  = 0;
  TBranch *b_postSkimCounter = 0;
  
  tree->SetBranchAddress("genEvents"             , &genEvents_           , &b_genEvents           );
  tree->SetBranchAddress("analyzedEvents"        , &analyzedEvents_      , &b_analyzedEvents      );
  tree->SetBranchAddress("internalCrossSection"  , &internalCrossSection_, &b_internalCrossSection);
  tree->SetBranchAddress("externalCrossSection"  , &externalCrossSection_, &b_externalCrossSection);
  tree->SetBranchAddress("summcprocweight"       , &summcprocweight_     , &b_summcprocweight     );  
  tree->SetBranchAddress("sumpuweight"           , &sumpuweight_         , &b_sumpuweight         );  
  tree->SetBranchAddress("sumpumcprocweight"     , &sumpumcprocweight_   , &b_sumpumcprocweight   );  
  tree->SetBranchAddress("preSkimCounter"        , &preSkimCounter_      , &b_preSkimCounter      );
  tree->SetBranchAddress("postSkimCounter"       , &postSkimCounter_     , &b_postSkimCounter     );

  Long64_t nentries = tree->GetEntries();  
  
  // temp variables
  double meanIntCrossSection = 0.;
  int    totalAnEvents       = 0 ;
  int    totalGenEvents      = 0 ;
  double totalSumMCProc      = 0.;
  double totalSumPUMCProc    = 0.;

  for (Long64_t jentry=0; jentry<nentries; ++jentry){
    tree->LoadTree(jentry); tree->GetEntry(jentry);
    
    if(genEvents_ != preSkimCounter_)
      std::cout << colour::Warning("WARNING! The number of skimmed events differs from the total generated events. Make sure you are properly weighting the events.") 
		<< " Generated events: " <<  genEvents_ << " Pre-skimmed events: " << preSkimCounter_ 
		<< std::endl;
    
    totalAnEvents       += analyzedEvents_;
    totalGenEvents      += genEvents_;
    totalSumMCProc      += summcprocweight_;
    totalSumPUMCProc    += sumpumcprocweight_;
    meanIntCrossSection += genEvents_*internalCrossSection_;
  }

  // The tree is not needed anymore
  delete tree;
  
  // ... and the variables used to set the branch address can be overwritten too
  genEvents_            = totalGenEvents;
  analyzedEvents_       = totalAnEvents;
  summcprocweight_      = totalSumMCProc;
  sumpumcprocweight_    = totalSumPUMCProc;
  internalCrossSection_ = meanIntCrossSection/genEvents_;

  crossSection_ = &internalCrossSection_;
  std::string xsectype = "internal";
  if(externalCrossSection_ > 0){
    crossSection_ = &externalCrossSection_;
    xsectype = "external";
  }

  sampleWeight_ = luminosity_*crossSection()/genEvents_;
  
  


  std::cout<<"\nThis sample has been made out of a dataset containing " << Green(genEvents()) << " generated events."   << std::endl
	   <<"Out of them, " << Green(analyzedEvents()) << " events have been used to produce the main tree."           << std::endl
	   <<"The cross-section of this sample is " << Green(crossSection()) << Green(" pb ") << "(" << xsectype <<")"  
	   <<" and the integrated luminosity scenario is "<< Green(luminosity_) << Green("/pb.")                        << std::endl
	   <<"The MC process event normalization is " << Green(analyzedEvents_/summcprocweight_)
	   <<" and the sample weight is " << Green(sampleWeight()) 
	   <<". The number of weighted events in the sample is " << Green(analyzedEventsWeighted()) << "."              << std::endl;
	   
}
