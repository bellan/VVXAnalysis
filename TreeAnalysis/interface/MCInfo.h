#ifndef VVXAnalysis_TreeAnalysis_MCInfo_H
#define VVXAnalysis_TreeAnalysis_MCInfo_H

#include "VVXAnalysis/DataFormats/interface/DiBoson.h"
#include "VVXAnalysis/DataFormats/interface/Lepton.h"

#include <string>

class MCInfo {

 public:
  MCInfo(const std::string &filename="", const double & lumi = -1, const double& externalXSection = -1.);
  ~MCInfo(){}

  int    genEvents()              const {return genEvents_;}
  int    analyzedEvents()         const {return analyzedEvents_;}
  double analyzedEventsWeighted() const {return sumpumcprocweight_ * sampleWeight();}

  double internalCrossSection() const {return internalCrossSection_;}
  double externalCrossSection() const {return externalCrossSection_;}
  double crossSection()         const {return *crossSection_;}
  double sampleWeight()         const {return sampleWeight_;}
  //double mcProcWeight()         const {return mcprocweight_*genEvents_/summcprocweight_;}
  double mcProcWeight()         const {return 1.;}
  double mcProcWeightNormalization() const {return genEvents_/summcprocweight_;}
  double mcProcWeightUnormalized()   const {return mcprocweight_;}
  double summcprocweight()      const {return summcprocweight_;}     
  double puWeight()             const {return puweight_;}

  // Total MC weight of the event. Beware, it does not include DATA/MC correction! See instead below.
  double weight()               const {return luminosity_ >= 0 ? sampleWeight_*mcProcWeight()*puweight_ : 1.;}

  // Total weight of the event, including efficiency scale factors.
  double weight(const phys::DiBoson<phys::Lepton, phys::Lepton> &ZZ) const {return ZZ.fakeRateSF() * (luminosity_ >= 0 ? weight() * ZZ.efficiencySF() : 1.);}
  
  double signalEfficiency()           const {return signalEfficiency_;}

  double signalEfficiencyCorrection() const {return 1./signalEfficiency_;}

  int    signalDefinition()           const {return signalDefinition_;}

 private:
  friend class EventAnalyzer;
  
  // integrated luminosity
  double luminosity_;
  double internalCrossSection_;
  double externalCrossSection_;
  double *crossSection_;
  double signalEfficiency_;
  int    signalDefinition_;
  int    genEvents_;
  int    analyzedEvents_;


  // Pure sample weight (from cross section and number of analyzed events)
  double sampleWeight_;
  // Intrinsic weight from MC sample
  double mcprocweight_;
  // Weight from PU reweighting. This is a per event weight
  double puweight_;

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

  int eventsIn2P2FCR_;
  int eventsIn3P1FCR_;
};


#endif
