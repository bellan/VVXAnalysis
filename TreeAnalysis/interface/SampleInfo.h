#ifndef VVXAnalysis_TreeAnalysis_SampleInfo_H
#define VVXAnalysis_TreeAnalysis_SampleInfo_H

#include "VVXAnalysis/DataFormats/interface/DiBoson.h"
#include "VVXAnalysis/DataFormats/interface/Lepton.h"
#include "VVXAnalysis/DataFormats/interface/GenEventWeights.h"
#include "VVXAnalysis/DataFormats/interface/RegionsCounter.h"

#include <string>

class TChain;

class SampleInfo {

 public:
  SampleInfo(const std::string &filename="", const double & lumi = -1, const double& externalXSection = -1., bool blinded = true, bool applyFRSF = true, bool forcePosWeight = false);
  ~SampleInfo(){}

  int    genEvents()              const {return genEvents_;}
  int    analyzedEvents()         const {return analyzedEvents_;}
  double analyzedEventsWeighted() const {return analyzedEvents() * sampleWeight();}

  inline phys::RegionsCounter eventsInRegions() const {return *eventsInRegions_;}
  
  double internalCrossSection() const {return internalCrossSection_;}
  double externalCrossSection() const {return externalCrossSection_;}
  double crossSection()         const {return *crossSection_;}
  double sampleWeight()         const {return sampleWeight_*mcWeight();}
  double mcWeight()             const {return genEventWeights_->mcProcWeight()*genEvents_/summcprocweight_;}
  //double mcWeight         const {return 1.;}
  double mcWeightNormalization() const {return genEvents_/summcprocweight_;}
  double mcWeightUnormalized()   const {return genEventWeights_->mcProcWeight();}
  double summcprocweight()      const {return summcprocweight_;}     
  double puWeight()             const {return genEventWeights_->puWeight();}
  double puWeightUncUp()        const {return genEventWeights_->puWeightUncUp();}
  double puWeightUncDn()        const {return genEventWeights_->puWeightUncDn();}
  double L1PrefiringWeight()    const {return genEventWeights_->L1PrefiringWeight();}
  double L1PrefiringWeightUp()  const {return genEventWeights_->L1PrefiringWeightUp();}
  double L1PrefiringWeightDn()  const {return genEventWeights_->L1PrefiringWeightDn();}

  

  // Total MC weight of the event. Beware, it does not include DATA/MC correction! See instead below.
    double weight()               const {return luminosity_ >= 0 ? sampleWeight()*puWeight()*L1PrefiringWeight() : 1.;}

  // Total weight of the event, including efficiency scale factors.
  double weight(const phys::DiBoson<phys::Lepton, phys::Lepton> &VV) const {
    double w = (applyFRSF_ ? VV.fakeRateSF() : 1.) * (luminosity_ >= 0 ? weight() * VV.efficiencySF() : 1.);
    if (forcePosWeight_) return abs(w);
    else return w;
  }
  double weight(const phys::Boson<phys::Lepton> &Z) const {
    double w = (applyFRSF_ ? Z.fakeRateSF() : 1.) * (luminosity_ >= 0 ? weight() * Z.efficiencySF() : 1.);
    if (forcePosWeight_) return abs(w);
    else return w;
  }

  double signalEfficiency()           const {return signalEfficiency_;}
  double signalEfficiencyCorrection() const {return 1./signalEfficiency_;}
  int    signalDefinition()           const {return signalDefinition_;}

  // Approximate number of events in SR and CRs. It is an approx number because the trigger requirement is not asked.
  //int approximateNeventsInSR()     const {return eventsInSR_;} // fixme


  int setup() const {return setup_;}
  
  bool isMC() const {return isMC_;}

  std::string fileName() const {return filename_; }
  //std::string fileName() const {return filename_;}
  
  
  double kF_ggZZ    () const {return genEventWeights_->kF_ggZZ()    ;} 
  double kF_qqZZM   () const {return genEventWeights_->kF_qqZZM()   ;}
  double kF_qqZZPt  () const {return genEventWeights_->kF_qqZZPt()  ;}
  double kF_qqZZdPhi() const {return genEventWeights_->kF_qqZZdPhi();}
  double kF_EWKqqZZ () const {return genEventWeights_->kF_EWKqqZZ() ;}



  float PDFScale           () const {return genEventWeights_->PDFScale           ();}
  float QCDscale_muR1F1    () const {return genEventWeights_->QCDscale_muR1F1    ();} 
  float QCDscale_muR1F2    () const {return genEventWeights_->QCDscale_muR1F2    ();}
  float QCDscale_muR1F0p5  () const {return genEventWeights_->QCDscale_muR1F0p5  ();}
  float QCDscale_muR2F1    () const {return genEventWeights_->QCDscale_muR2F1    ();}
  float QCDscale_muR2F2    () const {return genEventWeights_->QCDscale_muR2F2    ();}
  float QCDscale_muR2F0p5  () const {return genEventWeights_->QCDscale_muR2F0p5  ();}
  float QCDscale_muR0p5F1  () const {return genEventWeights_->QCDscale_muR0p5F1  ();}
  float QCDscale_muR0p5F2  () const {return genEventWeights_->QCDscale_muR0p5F2  ();}
  float QCDscale_muR0p5F0p5() const {return genEventWeights_->QCDscale_muR0p5F0p5();}
  float PDFVar_Up          () const {return genEventWeights_->PDFVar_Up          ();}
  float PDFVar_Down        () const {return genEventWeights_->PDFVar_Down        ();}
  float alphas_MZ_Up       () const {return genEventWeights_->alphas_MZ_Up       ();}
  float alphas_MZ_Down     () const {return genEventWeights_->alphas_MZ_Down     ();}

  // Not so elegant, as a class named SampleInfo is giving info about data
  void extractDataInfo(TChain *tree);

 private:
  friend class EventAnalyzer;

  std::string filename_;  
  double luminosity_;
  bool   blinded_;
  bool   applyFRSF_;
  bool   forcePosWeight_;
  
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

  phys::GenEventWeights *genEventWeights_;
  

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

  phys::RegionsCounter *eventsInRegions_;


  int setup_;
  
  //others utilities
  bool isMC_;
};


#endif
