#include "VVXAnalysis/TreeAnalysis/interface/MCInfo.h"
#include "VVXAnalysis/Commons/interface/Colours.h"

#include "TChain.h"
#include <iostream>
#include <bitset>

using namespace colour;

MCInfo::MCInfo(const std::string& filename, const double & lumi, const double& externalXSection)
  : luminosity_(lumi)
  , internalCrossSection_(-1)
  , externalCrossSection_(externalXSection)
  , crossSection_(&externalCrossSection_)
  , signalEfficiency_(0.)
  , signalDefinition_(-99)
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
  , signalCounter_(0)
  , postSkimSignalEvents_(0)
  , eventsInEtaAcceptance_(0)
  , eventsInEtaPtAcceptance_(0)
  , eventsIn2P2FCR_(0)
  , eventsIn3P1FCR_(0)
{

  TChain *tree = new TChain("treePlanter/HollyTree");
  tree->Add(filename.c_str());

  if (tree == 0) return;

  TBranch *b_eventsIn2P2FCR   = 0;
  TBranch *b_eventsIn3P1FCR   = 0;

  tree->SetBranchAddress("eventsIn2P2FCR", &eventsIn2P2FCR_, &b_eventsIn2P2FCR);
  tree->SetBranchAddress("eventsIn3P1FCR", &eventsIn3P1FCR_, &b_eventsIn3P1FCR);

  int totalEventsIn2P2FCR = 0;
  int totalEventsIn3P1FCR = 0;

  Long64_t nentries = tree->GetEntries();  

  for (Long64_t jentry=0; jentry<nentries; ++jentry){
    tree->LoadTree(jentry); tree->GetEntry(jentry);

    totalEventsIn2P2FCR += eventsIn2P2FCR_;
    totalEventsIn3P1FCR += eventsIn3P1FCR_;
  }

  eventsIn2P2FCR_ = totalEventsIn2P2FCR;
  eventsIn3P1FCR_ = totalEventsIn3P1FCR;

  std::cout<<"2p2f: "  <<eventsIn2P2FCR_<<" 3p1f: "<< eventsIn3P1FCR_ << std::endl;

  if(lumi <= 0) return; // FIXME: access here CR info
  
  
  TBranch *b_signalDefinition        = 0;
  TBranch *b_genEvents               = 0;
  TBranch *b_analyzedEvents          = 0;
  TBranch *b_internalCrossSection    = 0;  
  TBranch *b_externalCrossSection    = 0;  
  TBranch *b_summcprocweight         = 0;  
  TBranch *b_sumpuweight             = 0;  
  TBranch *b_sumpumcprocweight       = 0;
  TBranch *b_preSkimCounter          = 0;
  TBranch *b_postSkimCounter         = 0;
  TBranch *b_signalCounter           = 0;
  TBranch *b_postSkimSignalEvents    = 0;
  TBranch *b_eventsInEtaAcceptance   = 0;
  TBranch *b_eventsInEtaPtAcceptance = 0;

  tree->SetBranchAddress("signalDefinition"       , &signalDefinition_       , &b_signalDefinition       );
  tree->SetBranchAddress("genEvents"              , &genEvents_              , &b_genEvents              );
  tree->SetBranchAddress("analyzedEvents"         , &analyzedEvents_         , &b_analyzedEvents         );
  tree->SetBranchAddress("internalCrossSection"   , &internalCrossSection_   , &b_internalCrossSection   );
  tree->SetBranchAddress("externalCrossSection"   , &externalCrossSection_   , &b_externalCrossSection   );
  tree->SetBranchAddress("summcprocweight"        , &summcprocweight_        , &b_summcprocweight        );  
  tree->SetBranchAddress("sumpuweight"            , &sumpuweight_            , &b_sumpuweight            );  
  tree->SetBranchAddress("sumpumcprocweight"      , &sumpumcprocweight_      , &b_sumpumcprocweight      );  
  tree->SetBranchAddress("preSkimCounter"         , &preSkimCounter_         , &b_preSkimCounter         );
  tree->SetBranchAddress("postSkimCounter"        , &postSkimCounter_        , &b_postSkimCounter        );
  tree->SetBranchAddress("signalCounter"          , &signalCounter_          , &b_signalCounter          );
  tree->SetBranchAddress("postSkimSignalEvents"   , &postSkimSignalEvents_   , &b_postSkimSignalEvents   );
  tree->SetBranchAddress("eventsInEtaAcceptance"  , &eventsInEtaAcceptance_  , &b_eventsInEtaAcceptance  );
  tree->SetBranchAddress("eventsInEtaPtAcceptance", &eventsInEtaPtAcceptance_, &b_eventsInEtaPtAcceptance);
    


  
  // temp variables
  double meanIntCrossSection = 0.;
  int    totalAnEvents       = 0 ;
  int    totalGenEvents      = 0 ;
  double totalSumMCProc      = 0.;
  double totalSumPUMCProc    = 0.;
  int    totalAccEta         = 0 ;
  int    totalAccEtaPt       = 0 ;
  int    totalSignalEventsPreSel  = 0;
  int    totalSignalEventsPostSel = 0;

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
    totalAccEta         += eventsInEtaAcceptance_;
    totalAccEtaPt       += eventsInEtaPtAcceptance_;
    totalSignalEventsPreSel  += signalCounter_;
    totalSignalEventsPostSel += postSkimSignalEvents_; 
  }

  // The tree is not needed anymore
  delete tree;
  
  // ... and the variables used to set the branch address can be overwritten too
  genEvents_               = totalGenEvents;
  analyzedEvents_          = totalAnEvents;
  summcprocweight_         = totalSumMCProc;
  sumpumcprocweight_       = totalSumPUMCProc;
  internalCrossSection_    = meanIntCrossSection/genEvents_;
  eventsInEtaAcceptance_   = totalAccEta;
  eventsInEtaPtAcceptance_ = totalAccEtaPt;
  signalCounter_           = totalSignalEventsPreSel;
  postSkimSignalEvents_    = totalSignalEventsPostSel; 


  crossSection_ = &internalCrossSection_;
  std::string xsectype = "internal";
  if(externalCrossSection_ > 0){
    crossSection_ = &externalCrossSection_;
    xsectype = "external";
  }

  sampleWeight_ = luminosity_*crossSection()/genEvents_;
  
  double signalFraction     = double(signalCounter_)/genEvents_;
  double signalCrossSection = signalFraction * crossSection();
  signalEfficiency_         = double(postSkimSignalEvents_)/signalCounter_; // FIXME: double check this. That does not account for fakes.

  std::cout<<"\nThis sample has been made out of a dataset containing " << Green(genEvents()) << " generated events."   << std::endl
	   <<"Out of them, " << Green(analyzedEvents()) << " events have been used to produce the main tree."           << std::endl
	   <<"The cross-section of this sample is " << Green(crossSection()) << Green(" pb ") << "(" << xsectype <<")"  
	   <<" and the integrated luminosity scenario is "<< Green(luminosity_) << Green("/pb.")                        << std::endl
	   <<"The MC process event normalization is " << Green(analyzedEvents_/summcprocweight_)
	   <<" and the sample weight is " << Green(sampleWeight()) 
	   <<". The number of weighted events in the sample is " << Green(analyzedEventsWeighted()) << "."              << std::endl
	   << "The signal definition adopted for this analysis is " << Green(signalDefinition()) 
	   << " (" << Green(std::bitset<16>(signalDefinition())) << ")."                                                        << std::endl
	   <<"The fraction of the signal in the sample is " << Green(signalFraction)  
	   <<", that corresponds to a cross section of " <<  Green(signalCrossSection) << Green(" pb.")                 << std::endl
	   <<"The signal efficiency for the baseline selection is " << Green(signalEfficiency()) << "."
	   <<std::endl;
	   
}

