#ifndef VVGammaAnalyzer_h
#define VVGammaAnalyzer_h

/** \class VVGammaAnalyzer
 *  ZZG --> 4l G,   WZG --> 3l nu G,   VZG --> 2l 2j G
 *
 *  $Date: 2013/03/15 13:37:31 $
 *  $Revision: 1.4 $
 *  \author A. Mecca - UNITO <alberto.mecca@cern.ch>
 */

#include <set>
#include <map>
#include <unordered_map>
#include <fstream>

#include "EventAnalyzer.h"
#include "RegistrableAnalysis.h"
#include "VVXAnalysis/Commons/interface/Constants.h"

class VVGammaAnalyzer: public EventAnalyzer, RegistrableAnalysis<VVGammaAnalyzer>{
  class SignalDefinitionHelper {
    /*
      Aggregate the information about the signal definition status for the current event.
      Provides a common method pass() to check the appropriate logic
     */
  public:
    SignalDefinitionHelper(VVGammaAnalyzer* a) : analyzer_(a) {}
    bool pass() const;
    bool pass_photon() const { return photon_; };
    void eval(/*const VVGammaAnalyzer* analyzer*/);

    // This enum should probably be declared in Commons/interface/SignalDefinition.h
    enum TopologyBit : unsigned int {
      ZZto4L      = 0,
      WZto3LNu    = 1,
      ZtoLL       = 2,
      Fiducial    = 4,
      Detector    = 5,
      TrigPlateau = 6,
      VVmass100   = 7
    };

  protected:
    bool eval_photon();
    bool eval_ZZ4L  ();
    bool eval_WZ3L  ();
    bool eval_ZV2L2j();

  private:
    VVGammaAnalyzer* analyzer_;

    bool photon_;  // There is a gen photon that passes the signal definition
    bool passllLowMass_;
    bool ZZ4L_;
    bool WZ3L_;
    bool ZV2L2j_;
    bool Vhad_;
    bool Zll_;
  };

public:

  //, const std::string& filename, const double& lumi = 1., const double& externalXSection = -1., bool doBasicPlots = false

  VVGammaAnalyzer(const AnalysisConfiguration& configuration)
    : EventAnalyzer(*(new Selector<VVGammaAnalyzer>(*this)),
		    configuration)
    , sigdefHelper(this) {
    //theHistograms.profile(genCategory);
    // Memory allocation
    leptons_      = new std::vector<phys::Lepton>;
    jets_noph_ .reset(new std::vector<phys::Jet>);
    fsrPhotons_.reset(new std::vector<phys::Particle>);
    genQuarks_    = new std::vector<phys::Particle>;
    genChLeptons_ = new std::vector<phys::Particle>;
    genNeutrinos_ = new std::vector<phys::Particle>;
    genPhotons_   = new std::vector<phys::Particle>;
    genPhotonsPrompt_.reset(new std::vector<phys::Particle>);

    genZtoLepWithTau_.reset(new std::vector<phys::Boson<phys::Particle> >);
    genWtoLepWithTau_.reset(new std::vector<phys::Boson<phys::Particle> >);
  }

  virtual ~VVGammaAnalyzer(){
    delete leptons_;
    delete genQuarks_;
    delete genChLeptons_;
    delete genNeutrinos_;
    delete genPhotons_;
  }
	
  virtual void begin();
  	
  void initEvent();
  virtual Int_t cut();

  virtual void analyze();

  virtual void end(TFile &);
	
  virtual void finish();


  friend class Selector<VVGammaAnalyzer>;
  template< class PAR >
  bool ZBosonDefinition(phys::Boson<PAR> cand) const{
    bool checkCharge = cand.daughter(0).charge() + cand.daughter(1).charge() == 0;
    return checkCharge && fabs(cand.p4().M() - phys::ZMASS) < 30;
  }


  template< class PAR >
  bool WBosonDefinition(phys::Boson<PAR> cand) const {
    bool gooddaughters = (fabs(cand.daughter(0).eta()) < 2.5 && cand.daughter(0).pt() > 30 &&
			  cand.daughter(0).passPUID() && cand.daughter(0).passLooseJetID() &&
			  fabs(cand.daughter(1).eta()) < 2.5 && cand.daughter(1).pt() > 30 &&
			  cand.daughter(1).passPUID() && cand.daughter(1).passLooseJetID());
    bool goodmass = fabs(cand.p4().M() - phys::WMASS) < 50;
    return (goodmass && gooddaughters);
  }

  static bool GenWtoLNuDefinition(const phys::Boson<phys::Particle>& cand) {
    bool gooddaughters = (fabs(cand.daughter(0).eta()) < 2.5 && cand.daughter(0).pt() > 20);
    return gooddaughters;
  }

  static bool GenZtoLLDefinition(const phys::Boson<phys::Particle>& cand) {
    bool gooddaughters = (fabs(cand.daughter(0).eta()) < 2.5 && cand.daughter(0).pt() > 5 &&
			  fabs(cand.daughter(1).eta()) < 2.5 && cand.daughter(1).pt() > 5);
    bool goodmass = 60 < cand.p4().M() && cand.p4().M() < 120;
    return goodmass && gooddaughters;
  }

  static bool GenVtoQQDefinition(const phys::Boson<phys::Particle>& cand) {
    bool gooddaughters = (fabs(cand.daughter(0).eta()) < 4.7 && cand.daughter(0).pt() > 30 &&
			  fabs(cand.daughter(1).eta()) < 4.7 && cand.daughter(1).pt() > 30);
    bool goodmass = 60 < cand.p4().M() && cand.p4().M() < 120;
    return goodmass && gooddaughters;
  }

  bool isPhotonPrompt(const phys::Photon& ph, double tolerance=0.2) const {
    return genPhotonsPrompt_->size() > 0 && physmath::deltaR( *closestDeltaR(ph, *genPhotonsPrompt_), ph ) < tolerance;
  }
  
  
private:
	
  std::vector<phys::Lepton>* leptons_;
  std::unique_ptr<std::vector<phys::Jet>> jets_noph_;
	
  // Systematics: photons {EScale, ESigma} x {Up, Down} + {central}
  const std::vector<const char*> photonSystKeys_ = {"central", "EScale_Up", "EScale_Down", "ESigma_Up", "ESigma_Down"};
  std::unordered_map<std::string, std::unique_ptr<std::vector<phys::Photon>>> kinPhotons_;    // Only kinematic selection
  std::unordered_map<std::string, std::unique_ptr<std::vector<phys::Photon>>> loosePhotons_;  // Loose ID: currently 3/5 cuts of ID
  std::unordered_map<std::string, std::unique_ptr<std::vector<phys::Photon>>> goodPhotons_;   // Tight ID: currently Loose WP of POG cut-based ID
  std::unique_ptr<std::vector<phys::Particle>> fsrPhotons_;
 	
  // Vectors of gen particles
  std::vector<phys::Particle>* genQuarks_;
  std::vector<phys::Particle>* genChLeptons_;
  std::vector<phys::Particle>* genNeutrinos_;
  std::vector<phys::Particle>* genPhotons_;
  std::unique_ptr<std::vector<phys::Particle>> genPhotonsPrompt_;
  // Vectors of gen Bosons
  std::unique_ptr<std::vector<phys::Boson<phys::Particle> > > genZtoLepWithTau_;
  std::unique_ptr<std::vector<phys::Boson<phys::Particle> > > genWtoLepWithTau_;
  // Gen objects
  phys::DiBoson<phys::Particle, phys::Particle> genZZ_;
  phys::DiBoson<phys::Particle, phys::Particle> genZW_;
  // V --> j (j)
  phys::Boson<phys::Jet> candVTojj_;
  phys::Jet              candVToJ_ ;
  // KinPhoton that passes the largest number of cuts of the Loose ID
  phys::Photon* bestKinPh_;
  phys::Photon* bestMVAPh_;

  std::unique_ptr<TH2F> hPhotonFR_VLtoL_data_;
  std::unique_ptr<TH2F> hPhotonFR_VLtoL_dataZG_;
  std::unique_ptr<TH2F> hPhotonFR_KtoVLexcl_;
  std::unique_ptr<TH2F> hPhotonFRSF_VLtoL_;

  std::unique_ptr<TH2F> hPhotonEffSF_;
  double hPhotonEffSF_maxPt_;

  std::map<phys::Photon::MVAwp, std::unique_ptr<TH2F>> mapPhotonMVASF_;
  std::map<phys::Photon::MVAwp, float                > mapPhotonMVASF_maxPt_;

  std::string channelReco_;

  bool passFSRcut_;

  // Objects reconstruction for each event
  void makeChannelReco();  // sets channelReco_
  void genEventSetup();
  /* const phys::Jet* candAK8(const std::vector<phys::Jet>*) const; */
  void hadronicObjectsReconstruction();
  static std::vector<phys::Boson<phys::Jet>> candidatesVTojj(const std::vector<phys::Jet>&);
 	
  // Basic histograms
  void genEventHistos();
  void baseHistos_cut();
  void baseHistos_analyze();
  void fillPhotonPlots(const phys::Photon& ph, const char* name, const char* title);
  void photonHistos();
  void jetHistos();
  void PKU_comparison();
 	
  // Sub analyses
  void plotsVVGstatus(const char* name, const char* title, const TLorentzVector& p4_VV, const char* mType="mass");
  void leptonFakeRate();
  void photonFakeRate_ABCD();
  void photonFakeRate_LtoT(const char* method, const phys::Photon& thePh, bool isPass, double effSF);
  void photonFRClosure(const char* method, const phys::Photon& thePh, bool isPass, double f_FR);
  void studyFSRregion(const std::vector<phys::Photon>&);
  template<class PAR>
  void efficiency(const std::vector<PAR>& vRec, const std::vector<phys::Particle>& vGen, const char* recLabel, const char* genLabel, double tolerance=0.4);
  void photonIsolation(const std::vector<phys::Photon>&, const char*);
  void photonIsolation_bestKin();
  void orphanPhotonStudy();  // study reco photons that are not matched to gen
  void systematicsStudy(  const char* sys_label);
  void SYSplots_inclusive(const char* sys_label, const char* syst, double weight);
  void SYSplots_photon(   const char* sys_label, const char* syst, double weight, const phys::Photon& ph, const char* ph_selection);
  void SYSplots_phCut(    const char* sys_label, const char* syst, double weight, const phys::Photon& phCut);
  void SYSplots_phMVA(    const char* sys_label, const char* syst, double weight, const phys::Photon& phMVA);
  void SYSplots(          const char* sys_label, const char* syst, double weight, const phys::Photon* phCut, const phys::Photon* phMVA);
  void debug3Lregion();
  void photonGenStudy();
  void ZllVsZllGstudy(const std::vector<phys::Photon>&, const char*);
  void ZllVsZllGstudyGEN(const std::vector<phys::Particle>&, const char*);

  // void printCSVheader(std::ofstream&);
  // void printCSV(std::ofstream&);

  void endNameHistos();
 	
  // Utilities
  static std::unique_ptr<TH2F> getHistfromFile(const char* fname, const char* hname="PhFR", const char* info="");
  bool canBeFSR(const phys::Photon&, const std::vector<phys::Lepton>&) const;

  template <class T, class UnaryPredicate>
  static std::vector<phys::Boson<T>> makeBosons(const std::vector<T>&, UnaryPredicate); //predicate acts on the newly made bosons
  void initCherryPick();
  bool cherrypickEvt() const;
  std::map<unsigned long,                    // run
	   std::map<unsigned long,           // lumi block
		    std::set<unsigned long>  // event
		    >
	   >
  cherryEvents;

  void fillCutsNm1(const std::string& name, const std::string& title, const std::vector<std::pair<std::string, bool>>& cuts, const double& weight);
  void fillCutFlow(const std::string& name, const std::string& title, const std::vector<std::pair<std::string, bool>>& cuts, const double& weight);

  double getPhotonFR_VLtoL       (const phys::Photon& ph) const;
  double getPhotonFRUnc_VLtoL    (const phys::Photon& ph) const;
  double getPhotonFR_VLtoL_data  (const phys::Photon& ph) const;
  double getPhotonFR_VLtoL_dataZG(const phys::Photon& ph) const;
  double getPhotonFR_KtoVLexcl   (const phys::Photon& ph) const;
  double getPhotonFRSF_VLtoL     (const phys::Photon& ph) const;

  double getPhotonEffSF_MVA(const phys::Photon&, phys::Photon::MVAwp) const;
  double getPhotonEffSFUnc_MVA(const phys::Photon&, phys::Photon::MVAwp) const;

  int photonEffSF_getBin(const phys::Photon& ph) const{
    double pt = ph.pt() < hPhotonEffSF_maxPt_ ? ph.pt() : hPhotonEffSF_maxPt_ - 0.1;
    return hPhotonEffSF_->FindFixBin(ph.eta(), pt);
  }

  inline float getPhotonEffSF(   const phys::Photon& ph) const{
    return hPhotonEffSF_->GetBinContent(photonEffSF_getBin(ph));  // ph->efficiencySF();     // Note: to be restored when using new ntuples
  }

  inline float getPhotonEffSFUnc(const phys::Photon& ph) const{
    return hPhotonEffSF_->GetBinError(  photonEffSF_getBin(ph));  // ph->efficiencySF();     // Note: to be restored when using new ntuples
  }

  static bool is4Lregion(const phys::RegionTypes reg){
    return (reg == phys::SR4P || reg == phys::CR3P1F || reg == phys::CR2P2F ||
	    reg == phys::SR4P_1L || reg == phys::SR4P_1P || reg == phys::CR4P_1F);
  }
  static bool is3Lregion(const phys::RegionTypes reg){
    return (reg == phys::SR3P || reg == phys::CR110 || reg == phys::CR101 || reg == phys::CR011 ||
	    reg == phys::CR100 || reg == phys::CR010 || reg == phys::CR001 || reg == phys::CR000 ||
	    reg == phys::SR3P_1L || reg == phys::SR3P_1P || reg == phys::CR3P_1F);
  }
  static bool is2Lregion(const phys::RegionTypes reg){
    return (reg == phys::SR2P || reg == phys::SR2P_1L || reg == phys::SR2P_1P || reg == phys::CR2P_1F);
  }
	
  static bool passVeryLoose(const phys::Photon& ph);
	
  static char phABCD(const phys::Photon&, const phys::Photon::IdWp);
  static char phABCD_study(const phys::Photon&, const double& barrel_thr, const double& endcap_thr);
 	
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

  std::pair<double, double> getZllAndZllgMasses(const phys::Photon&);
  std::pair<double, double> getZllAndZllgMasses(const std::vector<phys::Photon>&);
		
  static const std::vector<double> pt_bins;
  static const std::vector<double> pt_bins_LFR;
  static const std::vector<double> eta_bins;
  static const std::vector<double> aeta_bins;
  static const std::vector<double> ph_pt_bins;
  static const std::vector<double> ph_ptExtended_bins;
  static const std::vector<double> ph_aeta_bins;
  static const std::vector<double> ph_aetaExtended_bins;
  static const std::vector<double> mVV_bins;
  static const std::vector<double> mVVG_bins;
  static const std::vector<double> mZ_bins;
  static const std::vector<double> mZG_bins;
  
protected:
  unsigned long evtN_;  // Used in the print in cut()
  std::map<phys::RegionTypes, unsigned long> evtNInReg_, analyzedNInReg_;  // Used to count processed events.
  clock_t startTime_;  // Used to calculate elapsed time
  std::map<phys::RegionTypes, float> evtWInReg_, analyzedWInReg_;  // Weighted events passing cut and total
  SignalDefinitionHelper sigdefHelper;
};
#endif

