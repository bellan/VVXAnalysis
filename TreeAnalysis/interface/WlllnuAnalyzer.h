#ifndef WlllnuAnalyzer_h
#define WlllnuAnalyzer_h

/** \class WlllnuAnalyzer
 *  Concrete class for Wlllnu analysis
 *
 *  $Date: 2013/03/15 13:37:31 $
 *  $Revision: 1.4 $
 *  \author A. Corrado - UNITO <arianna.corrado@cern.ch>
 */


#include "EventAnalyzer.h"
#include "RegistrableAnalysis.h"
#include "VVXAnalysis/Commons/interface/Constants.h"
#include "VVXAnalysis/DataFormats/interface/Boson.h"
#include "VVXAnalysis/DataFormats/interface/Particle.h"

using std::string;

class WlllnuAnalyzer: public EventAnalyzer, RegistrableAnalysis<WlllnuAnalyzer>{

public:

  //, const std::string& filename, const double& lumi = 1., const double& externalXSection = -1., bool doBasicPlots = false
  
  WlllnuAnalyzer(const AnalysisConfiguration& configuration)
   : EventAnalyzer(*(new Selector<WlllnuAnalyzer>(*this)), 
		   configuration){
    //theHistograms.profile(genCategory);
    // Memory allocation
    genQuarks_       .reset(new std::vector<phys::Particle>);
    genChLeptons_    .reset(new std::vector<phys::Particle>);
    //genElectrons_    .reset(new std::vector<phys::Particle>);
    //genMuons_        .reset(new std::vector<phys::Particle>);
    genNeutrinos_    .reset(new std::vector<phys::Particle>);
    genPhotons_      .reset(new std::vector<phys::Particle>);
    genPhotonsPrompt_.reset(new std::vector<phys::Particle>);

    genZlepCandidates_.reset(new std::vector<phys::Boson<phys::Particle>>);
    genWlepCandidates_.reset(new std::vector<phys::Boson<phys::Particle>>);
    genZhadCandidates_.reset(new std::vector<phys::Boson<phys::Particle>>);
    genWhadCandidates_.reset(new std::vector<phys::Boson<phys::Particle>>);
  }

  virtual ~WlllnuAnalyzer(){}

  virtual void begin();
  virtual void end(TFile &);
  
  virtual void analyze();
  
  virtual Int_t cut();
  typedef std::pair<bool,int> boolInt;
  
  void genEventSetup();
  void reconstructionLepCompatibility(std::vector<phys::Particle>*, std::vector<phys::Lepton>*, string, string);

  template<class T, class V>
  static bool haveCommonDaughter(const phys::Boson<T>& a, const phys::Boson<V>& b, const float tol=0.001){
    return (
	    (a.daughter(0).p4() - b.daughter(0).p4()).P() < tol ||
	    (a.daughter(0).p4() - b.daughter(1).p4()).P() < tol ||
	    (a.daughter(1).p4() - b.daughter(0).p4()).P() < tol ||
	    (a.daughter(1).p4() - b.daughter(1).p4()).P() < tol   );
  }
	
  template<class T1, class T2, class V1, class V2>
  static bool haveCommonDaughter(const phys::DiBoson<T1, T2>& a, const phys::DiBoson<V1, V2>& b, const float tol=0.001){
    return (
	    haveCommonDaughter(a.daughter(0), b.daughter(0), tol) ||
	    haveCommonDaughter(a.daughter(0), b.daughter(1), tol) ||
	    haveCommonDaughter(a.daughter(1), b.daughter(0), tol) ||
	    haveCommonDaughter(a.daughter(1), b.daughter(1), tol)   );
  }

  bool GenWtoLNuDefinition(phys::Boson<phys::Particle> cand) const {
    bool gooddaughters = (fabs(cand.daughter(0).eta()) < 2.5 && cand.daughter(0).pt() > 20);
    return gooddaughters;
  }

  bool GenZtoLLDefinition(phys::Boson<phys::Particle> cand) const {
    bool gooddaughters = (fabs(cand.daughter(0).eta()) < 2.5 && cand.daughter(0).pt() > 5 &&
			  fabs(cand.daughter(1).eta()) < 2.5 && cand.daughter(1).pt() > 5);
    bool goodmass = 60 < cand.p4().M() && cand.p4().M() < 120;
    return goodmass && gooddaughters;
  }

  bool GenVtoQQDefinition(phys::Boson<phys::Particle> cand) const {
    bool gooddaughters = (fabs(cand.daughter(0).eta()) < 4.7 && cand.daughter(0).pt() > 30 &&
			  fabs(cand.daughter(1).eta()) < 4.7 && cand.daughter(1).pt() > 30);
    bool goodmass = 60 < cand.p4().M() && cand.p4().M() < 120;
    return goodmass && gooddaughters;
  }

  bool isPhotonPrompt(const phys::Photon& ph, double tolerance=0.2) const {
    return genPhotonsPrompt_->size() > 0 && physmath::deltaR( *closestDeltaR(ph, *genPhotonsPrompt_), ph ) < tolerance;
  }

private:
  
  friend class Selector<WlllnuAnalyzer>;
  template< class PAR >
    bool ZBosonDefinition(phys::Boson<PAR> cand) const{
    bool checkCharge = cand.daughter(0).charge() + cand.daughter(1).charge() == 0;
    return checkCharge && fabs(cand.p4().M() - phys::ZMASS) < 30;
  }

  
  template< class PAR >
    bool WBosonDefinition(phys::Boson<PAR> cand) {
    
    bool gooddaughters = false;
    if(fabs(cand.daughter(0).eta()) < 2.5 && cand.daughter(0).pt() > 30 &&
       cand.daughter(0).passPUID() && cand.daughter(0).passLooseJetID() &&
       fabs(cand.daughter(1).eta()) < 2.5 && cand.daughter(1).pt() > 30 &&
       cand.daughter(1).passPUID() && cand.daughter(1).passLooseJetID())
      gooddaughters = true;
    
    if(fabs(cand.p4().M() - phys::WMASS) < 150 && gooddaughters)
      return true;
    return false;
    
  }
  
  // Sample variables
  double genElPairInvMass;
  double genMuPairInvMass;
  double genFourLepInvMass;
  double genFourLepTransverseMass;
  double genLepPairTransverseMass;
  
  phys::Particle genMu1;
  phys::Particle genMu2;
  phys::Particle genEl;
  phys::Particle genNu_2;
  
  // Check minimum invariant mass of the possible l+l- pair -- //
  std::vector<phys::Boson<phys::Particle> > possibleLepPair(std::vector<phys::Particle>*);
  phys::Boson<phys::Particle> minInvMassChLepPair(std::vector<phys::Boson<phys::Particle> >);
  phys::Boson<phys::Particle> maxInvMassChLepPair(std::vector<phys::Boson<phys::Particle> >);
  // Check W to 3 leptons decay configuration
  int eventMode(std::vector<phys::Particle>*, std::vector<phys::Particle>*, std::vector<phys::Lepton>*, std::vector<phys::Lepton>*);
  // Check if there are Gen or/and Rec events
  bool genEvents;
  bool recEvents;
  // Check three leptons of same flavour charge
  bool checkLeptonsCharge(phys::Particle, phys::Particle, phys::Particle);
	
  // Efficiency parameters
  bool isGen_mode1(double);
  bool isRec_mode1(double);
  bool isRec_mode11(double);
  bool isRec_mode111(double);
  bool isRec_mode1111(double);
  
  bool isGen_mode2(double);
  bool isRec_mode2(double);
  bool isRec_mode22(double);
  bool isRec_mode222(double);
  bool isRec_mode2222(double);
  
  bool isGen_mode3(double);
  
  bool isGen_mode4(double);
  
  // Vectors of gen particles
  std::unique_ptr<std::vector<phys::Particle>> genQuarks_;
  std::unique_ptr<std::vector<phys::Particle>> genChLeptons_;
  std::unique_ptr<std::vector<phys::Particle>> genNeutrinos_;
  std::unique_ptr<std::vector<phys::Particle>> genPhotons_;
  std::unique_ptr<std::vector<phys::Particle>> genPhotonsPrompt_;
  std::vector<phys::Particle>* genElectrons_;
  std::vector<phys::Particle>* genMuons_;
  // Vectors of gen Bosons
  std::unique_ptr<std::vector<phys::Boson<phys::Particle> > > genZlepCandidates_;
  std::unique_ptr<std::vector<phys::Boson<phys::Particle> > > genWlepCandidates_;
  std::unique_ptr<std::vector<phys::Boson<phys::Particle> > > genZhadCandidates_;
  std::unique_ptr<std::vector<phys::Boson<phys::Particle> > > genWhadCandidates_;
  // Gen objects
  phys::DiBoson<phys::Particle, phys::Particle> genZZ_;
  phys::DiBoson<phys::Particle, phys::Particle> genZW_;
  
  // Vector of leptons
  std::vector<phys::Lepton>* leptons;
  
};
#endif

