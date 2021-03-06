#ifndef ZVAnalyzer_h
#define ZVAnalyzer_h

/** \class ZVAnalyzer
 *  Concrete class for VVX analysis
 *
 *  $Date: 2013/03/15 13:37:31 $
 *  $Revision: 1.4 $
 *  \author D. Usseglio - UNITO <davide.usseglio@edu.unito.it>
 */


#include "EventAnalyzer.h"
#include "RegistrableAnalysis.h"
#include "VVXAnalysis/Commons/interface/Constants.h"
#include "VVXAnalysis/Commons/interface/Comparators.h"
#include "TGraphAsymmErrors.h"

class ZVAnalyzer: public EventAnalyzer, RegistrableAnalysis<ZVAnalyzer>{

public:

  //, const std::string& filename, const double& lumi = 1., const double& externalXSection = -1., bool doBasicPlots = false

 ZVAnalyzer(const AnalysisConfiguration& configuration)
   : EventAnalyzer(*(new Selector<ZVAnalyzer>(*this)),
		   configuration){
    //theHistograms.profile(genCategory);
  }

  virtual ~ZVAnalyzer(){}
    
  virtual void begin();

  virtual void analyze();
  
  virtual Int_t cut();
    
  virtual void end(TFile &);


 private:
  friend class Selector<ZVAnalyzer>;
  template< class PAR >
    bool ZBosonDefinition(phys::Boson<PAR> cand) const{
    bool checkCharge = cand.daughter(0).charge() + cand.daughter(1).charge() == 0;
    bool checkMass = fabs(cand.p4().M() - phys::ZMASS) < 30;
    return checkCharge && checkMass;
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
    
    template <class P, class T>
    bool checkMatch(const P&, T&, const float&);
   
   template<class P>
    phys::Particle* findMatchingParticle(const phys::Particle&,std::vector<P>*);
    
    double minDeltaR;
    double minPos;
    
    void fillBasicPlots();
    template<class P>
    void fillParticlePlots(const std::string &, const P &);
    void normalizeHistograms(std::string);
    TGraphAsymmErrors* HistogramsErrors(TH1 *, TH1 *);
    
    void doSomeFits(std::string);
    void getFitInfo(TF1*);
    
    void initStatistics();
    void tempStatisticParticles(const phys::Particle&);
    void tempStatisticEvents();
    
    float counter;
    int eCounter;
    int mCounter;
    
    int promptCounter;
    int peCounter;
    int pmCounter;
    
    long electronEvents;
    long muonEvents;
    long passingSelection;
    float totalEvents;
    float badevents;
    
    float matchedElectrons;
    float matchedMuons;
    float matchedJets;
    float totalElectrons;
    float totalMuons;
    float totalJets;
    float totalJets2;
    
    
    float Counter;
    long checkElectrons;
    long checkMuons;
    
    template <class T, class P, typename C>
    void analyzeEfficiency(std::vector<T>*, std::vector<P>*, std::string, C&, const float&);
    
    template <class T, class P>
    void analyzeResolutionpt(const T&, const P&, std::string);
    
    template <class T, class P>
    void analyzeResolutionEnergy(const T&, const P&, std::string);
    
    void GenAnalysis();
    void endGenAnalysis();
    
    
    clock_t startTime;
};
#endif

