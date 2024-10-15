#ifndef VZGAnalyzer_h
#define VZGAnalyzer_h

/** \class WZZAnalyzer
 *  Concrete class for VZG analysis
 *
 *  $Date: 2013/03/15 13:37:31 $
 *  $Revision: 1.6 $
 *  \author R. Bellan - UNITO <riccardo.bellan@cern.ch>
 */


#include "EventAnalyzer.h"
#include "RegistrableAnalysis.h"
#include "VVXAnalysis/Commons/interface/Constants.h"

#include <TString.h>
#include <TTree.h>

class VZGAnalyzer: public EventAnalyzer, RegistrableAnalysis<VZGAnalyzer>{

public:

  //, const std::string& filename, const double& lumi = 1., const double& externalXSection = -1., bool doBasicPlots = false

//  VZGAnalyzer(const AnalysisConfiguration& configuration)
//    : EventAnalyzer(*(new Selector<VZGAnalyzer>(*this)), 
// 		   configuration){
//     //theHistograms.profile(genCategory);
//   }

VZGAnalyzer(const AnalysisConfiguration& configuration)
    : EventAnalyzer(*(new Selector<VZGAnalyzer>(*this)),
		    configuration){
    //theHistograms.profile(genCategory);
    // Memory allocation
    doFeats_ = true;//CT new flag for handling the features 
    leptons_      = new std::vector<phys::Lepton>;
    genQuarks_    = new std::vector<phys::Particle>;
    genChLeptons_ = new std::vector<phys::Particle>;
    genNeutrinos_ = new std::vector<phys::Particle>;
    genPhotons_   = new std::vector<phys::Particle>;
    genPhotonsPrompt_.reset(new std::vector<phys::Particle>);
    
    genZlepCandidates_ = new std::vector<phys::Boson<phys::Particle>>;
    genWlepCandidates_ = new std::vector<phys::Boson<phys::Particle>>;
    genZhadCandidates_ = new std::vector<phys::Boson<phys::Particle>>;
    genWhadCandidates_ = new std::vector<phys::Boson<phys::Particle>>;
  }

  virtual ~VZGAnalyzer(){
    delete leptons_;
    delete genQuarks_;
    delete genChLeptons_;
    delete genNeutrinos_;
    delete genPhotons_;
    
    delete genZlepCandidates_;
    delete genWlepCandidates_;
    delete genZhadCandidates_;
    delete genWhadCandidates_;
  }
  //virtual ~VZGAnalyzer(){}
  
  void PlotJet(const phys::Particle &, std::string , const float , std::string );
  
  void PlotJets(const phys::Particle &, const phys::Particle &, std::string , const float , std::string );
  
  void ResolutionPlots(const phys::Particle &, const phys::Particle &, std::string , const float , std::string );

  void genEventSetup();

  virtual void analyze();

  virtual void  fillFeatTree(FeatList&, bool&);

  virtual void genAnalyze();
  
  virtual void QuarksToJets();
  
  //virtual void AlternativegenAnalyze();
  
  virtual void genVBAnalyzer();
  
  virtual void genPhotonsAnalyzer();
  
  //virtual Int_t GENsignalConstraint();
    
  virtual bool baselineRequirements();

  virtual bool   IN_GENsignalDef();
  
  virtual bool   LeptonicSignalConstraint();
  
  virtual bool   HadronicSignalConstraint();
  
  virtual bool   PhotonSignalConstraint();

  virtual Bool_t cut(Int_t, phys::Boson<phys::Jet>,phys::Jet,std::vector<phys::Photon>,int);

  //  virtual void Reconstruct(phys::Boson<phys::Jet>*,phys::Jet*,bool*,bool*);
  int Reconstruct(phys::Boson<phys::Jet>*,phys::Jet*,bool*,bool*,phys::Photon*);
  
  virtual void PhotonSelection(std::vector<phys::Photon> *);

  virtual void CompatibilityTest(phys::Boson<phys::Jet>, phys::Boson<phys::Particle>, std::string, std::string);

  virtual void printHistos(uint, std::string, phys::Boson<phys::Jet>,phys::Jet,std::vector<phys::Photon>,int);

 private:
  std::vector<phys::Lepton>* leptons_;
	
  // Systematics: photons {EScale, ESigma} x {Up, Down} + {central}
  std::map<const char*, std::unique_ptr<std::vector<phys::Photon>>> kinPhotons_;    // Only kinematic selection
  std::map<const char*, std::unique_ptr<std::vector<phys::Photon>>> loosePhotons_;  // Loose ID: currently 3/5 cuts of ID
  std::map<const char*, std::unique_ptr<std::vector<phys::Photon>>> goodPhotons_;   // Tight ID: currently Loose WP of POG cut-based ID
  // Vectors of gen particles
  std::vector<phys::Particle>* genQuarks_;
  std::vector<phys::Particle>* genChLeptons_;
  std::vector<phys::Particle>* genNeutrinos_;
  std::vector<phys::Particle>* genPhotons_;
  std::unique_ptr<std::vector<phys::Particle>> genPhotonsPrompt_;
  // Vectors of gen Bosons
  std::vector<phys::Boson<phys::Particle>>* genZlepCandidates_;
  std::vector<phys::Boson<phys::Particle>>* genWlepCandidates_;
  std::vector<phys::Boson<phys::Particle>>* genZhadCandidates_;
  std::vector<phys::Boson<phys::Particle>>* genWhadCandidates_;
 // Gen objects
  phys::DiBoson<phys::Particle, phys::Particle> genZZ_;
  phys::DiBoson<phys::Particle, phys::Particle> genWZ_;
  // V --> j (j)
  phys::Boson<phys::Jet> candVTojj_;
  phys::Jet              candVToJ_ ;
  // KinPhoton that passes the largest number of cuts of the Loose ID
  phys::Photon* bestKinPh_;

  std::unique_ptr<TH2F> hPhotonFR_;
  std::unique_ptr<TH2F> hPhotonFR_KtoVL_;
  std::unique_ptr<TH2F> hPhotonFR_KtoVLexcl_;
  std::unique_ptr<TH2F> hPhotonFRSF_LtoT_;
  std::string channelReco_;



  template<class T, class V>
  bool haveCommonDaughter(const phys::Boson<T>& a, const phys::Boson<V>& b, const float tol=0.001){
    return (
	    (a.daughter(0).p4() - b.daughter(0).p4()).P() < tol ||
	    (a.daughter(0).p4() - b.daughter(1).p4()).P() < tol ||
	    (a.daughter(1).p4() - b.daughter(0).p4()).P() < tol ||
	    (a.daughter(1).p4() - b.daughter(1).p4()).P() < tol   );
  }
	
  template<class T1, class T2, class V1, class V2>
  bool haveCommonDaughter(const phys::DiBoson<T1, T2>& a, const phys::DiBoson<V1, V2>& b, const float tol=0.001){
    return (
	    haveCommonDaughter(a.daughter(0), b.daughter(0), tol) ||
	    haveCommonDaughter(a.daughter(0), b.daughter(1), tol) ||
	    haveCommonDaughter(a.daughter(1), b.daughter(0), tol) ||
	    haveCommonDaughter(a.daughter(1), b.daughter(1), tol)   );
  }

  
  friend class Selector<VZGAnalyzer>;

template< class PAR >
  bool ZBosonDefinition(phys::Boson<PAR> cand) const{
    bool checkCharge = cand.daughter(0).charge() + cand.daughter(1).charge() == 0;
    return checkCharge && fabs(cand.p4().M() - phys::ZMASS) < 30;
  }


  template< class PAR >
  bool WBosonDefinition(phys::Boson<PAR> cand) {
    bool gooddaughters = (fabs(cand.daughter(0).eta()) < 2.5 && cand.daughter(0).pt() > 30 &&
			  cand.daughter(0).passPUID() && cand.daughter(0).passLooseJetID() &&
			  fabs(cand.daughter(1).eta()) < 2.5 && cand.daughter(1).pt() > 30 &&
			  cand.daughter(1).passPUID() && cand.daughter(1).passLooseJetID());
    bool goodmass = fabs(cand.p4().M() - phys::WMASS) < 50;
    return (goodmass && gooddaughters);
  }
  
  bool GenWBosonDefinition(phys::Boson<phys::Particle> cand) {
    bool gooddaughters = (fabs(cand.daughter(0).eta()) < 2.5 && cand.daughter(0).pt() > 30 &&
			  fabs(cand.daughter(1).eta()) < 2.5 && cand.daughter(1).pt() > 30);
    bool goodmass = fabs(cand.p4().M() - phys::WMASS) < 50;
    return (goodmass && gooddaughters);
  }

};
#endif

