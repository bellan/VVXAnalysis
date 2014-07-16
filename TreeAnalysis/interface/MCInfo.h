#ifndef VVXAnalysis_TreeAnalysis_MCInfo_H
#define VVXAnalysis_TreeAnalysis_MCInfo_H

#include "VVXAnalysis/DataFormats/interface/DiBoson.h"
#include "VVXAnalysis/Commons/interface/LeptonEfficiency.h"

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
  double mcProcWeight()         const {return mcprocweight_*analyzedEvents_/summcprocweight_;}
  // Total weight of the event, computed on the fly
  double weight()               const {return luminosity_ >= 0 ? sampleWeight_*mcProcWeight()*puweight_ : 1.;}
  // Intermediate step, later the efficiency weight will be taken directly from the event
  double leptonEfficiencyWeight(const phys::DiBoson<phys::Lepton, phys::Lepton> &ZZ) const {return leptonEfficiency_.weight(ZZ);}
  double weight(const phys::DiBoson<phys::Lepton, phys::Lepton> &ZZ) const {return luminosity_ >= 0 ? weight()*leptonEfficiencyWeight(ZZ) : 1.;}
  
 private:
  friend class EventAnalyzer;
  
  // To get Lepton efficiency scale factors. Temporary here!
  LeptonEfficiency leptonEfficiency_;

  // integrated luminosity
  double luminosity_;
  double internalCrossSection_;
  double externalCrossSection_;
  double *crossSection_;
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

};


#endif
