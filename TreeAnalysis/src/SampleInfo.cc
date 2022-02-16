#include "VVXAnalysis/TreeAnalysis/interface/SampleInfo.h"
#include "VVXAnalysis/Commons/interface/Colours.h"

#include "TChain.h"
#include <iostream>
#include <bitset>

using namespace colour;

SampleInfo::SampleInfo(const std::string& filename, const double & lumi, const double& externalXSection, bool blinded, bool applyFRSF, bool forcePosWeight)
  : luminosity_(lumi)
  , blinded_(blinded)
  , applyFRSF_(applyFRSF)
  , forcePosWeight_(forcePosWeight)
  , internalCrossSection_(-1)
  , externalCrossSection_(-1)
  , externalCrossSectionFromCSV_(externalXSection)
  , crossSection_(&externalCrossSection_)
  , signalEfficiency_(0.)
  , signalDefinition_(-99)
  , genEvents_(-1)
  , analyzedEvents_(-1)
  , sampleWeight_(1)
    //  , mcprocweight_(1)
  , summcprocweight_(0)
  , sumpuweight_(0)
  , sumpumcprocweight_(0)
  , preSkimCounter_(0) 
  , postSkimCounter_(0)
  , signalCounter_(0)
  , postSkimSignalEvents_(0)
  , eventsInEtaAcceptance_(0)
  , eventsInEtaPtAcceptance_(0)
{

  genEventWeights_ = new phys::GenEventWeights();
  eventsInRegions_ = new phys::RegionsCounter(); 

  TChain *tree = new TChain("treePlanter/HollyTree");
  tree->Add(filename.c_str());

  if(tree == 0 || lumi <= 0) {isMC_=false; if(tree != 0) extractDataInfo(tree); return;}

  Long64_t nentries = tree->GetEntries();  

  TBranch *b_eventsInRegions         = 0;
  TBranch *b_setup                   = 0;
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

  tree->SetBranchAddress("eventsInRegions", &eventsInRegions_    , &b_eventsInRegions);
  tree->SetBranchAddress("setup"         , &setup_         , &b_setup);
  
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
    
  filename_ = filename;
  filename_.erase(0, filename_.find("/")+1); // erase samples/ 
  filename_.erase(0, 5);  // erase year/
  filename_.erase(filename_.find(".root")); 

  if((filename_.find("Single") != std::string::npos) || (filename_.find("Double") != std::string::npos) || (filename_.find("MuonEG") != std::string::npos)) isMC_ = kFALSE;
  else isMC_=kTRUE;

  // temp variables

  phys::RegionsCounter    *totalEventsInRegions = new phys::RegionsCounter();

  double meanIntCrossSection = 0.;
  int    totalPreSkimCounter = 0 ;
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
    
    *totalEventsInRegions += *eventsInRegions_;

    totalPreSkimCounter += preSkimCounter_;
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
  *eventsInRegions_     = *totalEventsInRegions;
  eventsInRegions_->unblind();
  
  genEvents_               = totalGenEvents;
  preSkimCounter_          = totalPreSkimCounter;

  if(genEvents_ != preSkimCounter_)
    std::cout << "\n" << colour::Warning("WARNING! The number of events before any skim (") << colour::Warning(preSkimCounter_) 
	      << colour::Warning(") differs from the total generated events (") <<colour::Warning(genEvents_) 
	      << colour::Warning("). Make sure you are properly weighting the events.") << std::endl
	      << "If you are running on ZZJetsTo4L without taus in the samples, you can safely ignore this message."
	      << std::endl;
  

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
  if(externalCrossSectionFromCSV_ > 0){
    crossSection_ = &externalCrossSectionFromCSV_;
    xsectype = "externalFromCSV";
  }

  sampleWeight_ = luminosity_*crossSection()/genEvents_;
  
  double signalFraction     = double(signalCounter_)/genEvents_;
  double signalCrossSection = signalFraction * crossSection();
  //signalEfficiency_         = double(eventsInSR_)/signalCounter_; // FIXME: double check this. That does not account for fakes nor trigger selections.

  std::cout<<"\nThis sample has been made out of a dataset containing " << Green(genEvents()) << " generated events."   << std::endl
	   <<"Out of them, " << Green(analyzedEvents()) << " events have been used to produce the main tree."           << std::endl
	   <<"The cross-section of this sample is " << Green(crossSection()) << Green(" pb ") << "(" << xsectype <<")"  
	   <<" and the integrated luminosity scenario is "<< Green(luminosity_) << Green("/pb.")                        << std::endl
	   <<"The MC process event normalization is " << Green(mcWeightNormalization())
	   <<" and the sample weight is " << Green(sampleWeight()) 
	   <<". The number of weighted events in the sample is (approx.) " << Green(analyzedEventsWeighted()) << "."    << std::endl
	   << "The signal definition adopted for this analysis is " << Green(signalDefinition()) 
	   << " (" << Green(std::bitset<16>(signalDefinition())) << ")."                                                << std::endl
	   <<"The fraction of the signal in the sample is " << Green(signalFraction)  
	   <<", that corresponds to a cross section of " <<  Green(signalCrossSection) << Green(" pb.")                 << std::endl

    // FIXME
    //<<"The efficiency for the baseline selection is " << eventsInSR_  <<"/"<< signalCounter_ << " = " << Green(signalEfficiency()) << ".\n"
	   <<"Events in the Regions"    << std::endl
	   <<*eventsInRegions_;
    //	   <<", in the 2P2F CR = " << Green(eventsIn2P2FCR_)
    //	   <<", in the 3P1F CR = " << Green(eventsIn3P1FCR_)
    //	   <<"."
	   
}

void SampleInfo::extractDataInfo(TChain *tree){
  
  Long64_t nentries = tree->GetEntries();  

  TBranch *b_eventsInRegions         = 0;
  TBranch *b_analyzedEvents          = 0;

  tree->SetBranchAddress("eventsInRegions", &eventsInRegions_, &b_eventsInRegions);
  tree->SetBranchAddress("analyzedEvents" , &analyzedEvents_ , &b_analyzedEvents);

  phys::RegionsCounter    *totalEventsInRegions = new phys::RegionsCounter();

  int    totalAnEvents       = 0 ;

  for (Long64_t jentry=0; jentry<nentries; ++jentry){
    tree->LoadTree(jentry); tree->GetEntry(jentry);
    
    *totalEventsInRegions += *eventsInRegions_;
    totalAnEvents       += analyzedEvents_;
  }

  // The tree is not needed anymore
  delete tree;
  
  // ... and the variables used to set the branch address can be overwritten too
  *eventsInRegions_     = *totalEventsInRegions;
  analyzedEvents_       = totalAnEvents;

  if(blinded_) eventsInRegions_->blind();
  else         eventsInRegions_->unblind();

  std::cout<<"\nThe number of analyzed events to produce this tree is " << Green(analyzedEvents()) << "."           << std::endl
	   <<"The selected events distribute in the Control/Search Regions as follow:"    << std::endl
	   <<eventsInRegions();
}

