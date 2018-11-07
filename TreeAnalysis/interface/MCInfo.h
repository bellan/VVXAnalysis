#ifndef VVXAnalysis_TreeAnalysis_MCInfo_H
#define VVXAnalysis_TreeAnalysis_MCInfo_H

#include "VVXAnalysis/DataFormats/interface/DiBoson.h"
#include "VVXAnalysis/DataFormats/interface/Lepton.h"
#include "VVXAnalysis/DataFormats/interface/GenEventWeights.h"

#include <string>

class MCInfo {

 public:
  MCInfo(const std::string &filename="", const double & lumi = -1, const double& externalXSection = -1.);
  ~MCInfo(){}

  int    genEvents()              const {return genEvents_;}
  int    analyzedEvents()         const {return analyzedEvents_;}
  double analyzedEventsWeighted() const {return analyzedEvents() * sampleWeight();}

  double internalCrossSection() const {return internalCrossSection_;}
  double externalCrossSection() const {return externalCrossSection_;}
  double crossSection()         const {return *crossSection_;}
  double sampleWeight()         const {return sampleWeight_*mcWeight();}
  double mcWeight()             const {return genEventWeights_.mcProcWeight()*genEvents_/summcprocweight_;}
  //double mcWeight         const {return 1.;}
  double mcWeightNormalization() const {return genEvents_/summcprocweight_;}
  double mcWeightUnormalized()   const {return genEventWeights_.mcProcWeight();}
  double summcprocweight()      const {return summcprocweight_;}     
  double puWeight()             const {return genEventWeights_.puWeight();}
  double puWeightUncUp()        const {return genEventWeights_.puWeightUncUp();}
  double puWeightUncDn()        const {return genEventWeights_.puWeightUncDn();}

  // Total MC weight of the event. Beware, it does not include DATA/MC correction! See instead below.
  double weight()               const {return luminosity_ >= 0 ? sampleWeight()*puWeight() : 1.;}

  // Total weight of the event, including efficiency scale factors.
  double weight(const phys::DiBoson<phys::Lepton, phys::Lepton> &ZZ) const {return ZZ.fakeRateSF() * (luminosity_ >= 0 ? weight() * ZZ.efficiencySF() : 1.);}
  double signalEfficiency()           const {return signalEfficiency_;}
  double signalEfficiencyCorrection() const {return 1./signalEfficiency_;}
  int    signalDefinition()           const {return signalDefinition_;}

  // Approximate number of events in SR and CRs. It is an approx number because the trigger requirement is not asked.
  int approximateNeventsInSR()     const {return eventsInSR_;}
  int approximateNeventsIn2P2FCR() const {return eventsIn2P2FCR_;}
  int approximateNeventsIn3P1FCR() const {return eventsIn3P1FCR_;}

  bool isMC() const {return isMC_;}

  std::string fileName() const {return filename_; }
  //std::string fileName() const {return filename_;}
  
  
  double kF_ggZZ    () const {return genEventWeights_.kFactor_ggZZ_    ;} 
  double kF_qqZZM   () const {return genEventWeights_.kFactor_qqZZM_   ;}
  double kF_qqZZPt  () const {return genEventWeights_.kFactor_qqZZPt_  ;}
  double kF_qqZZdPhi() const {return genEventWeights_.kFactor_qqZZdPhi_;}
  double kF_EWKqqZZ () const {return genEventWeights_.kFactor_EWKqqZZ_ ;}



 private:
  friend class EventAnalyzer;

  std::string filename_;  
  // integrated luminosity
  double luminosity_;
  double internalCrossSection_;
  double externalCrossSection_;
  double externalCrossSectionFromCSV_;
  double *crossSection_;
  double signalEfficiency_;
  int    signalDefinition_;
  int    genEvents_;
  int    analyzedEvents_;


  // Pure sample weight (from cross section and number of analyzed events)
  double sampleWeight_;
  // Intrinsic weight from MC sample
  //double mcprocweight_;

  phys::GenEventWeights genEventWeights_;
  
  // Weight from PU reweighting. This is a per event weight
  //  double puweight_;
  // double puweightUp_;
  // double puweightDn_;

  // Sum of weights
  double summcprocweight_;
  double sumpuweight_;
  double sumpumcprocweight_;
  
  // Counters for skims
  int preSkimCounter_; 
  int postSkimCounter_;
  int signalCounter_;

  int postSkimSignalEvents_;
  int eventsInEtaAcceptance_;
  int eventsInEtaPtAcceptance_;

  int eventsInSR_;
  int eventsIn2P2FCR_;
  int eventsIn3P1FCR_;

  //others utilities
  bool isMC_;
};


#endif
