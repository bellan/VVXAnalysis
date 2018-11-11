#ifndef WZAnalyzer_h
#define WZAnalyzer_h

/** \class WZAnalyzer
 *  Concrete class for WZ analysis
 *
 *  $Date: 2017/05/24 $
 *  $Revision: 0.5 $
 *  \author E. Racca - UNITO <eleonora.racca@cern.ch>
 */


#include "EventAnalyzer.h"
#include "RegistrableAnalysis.h"
#include "VVXAnalysis/Commons/interface/AriEle.h"
#include "VVXAnalysis/Commons/interface/Comparators.h"
#include "VVXAnalysis/Commons/interface/Constants.h"
#include "VVXAnalysis/Commons/interface/LeptonScaleFactors.h"
#include "VVXAnalysis/DataFormats/interface/Boson.h"
#include "VVXAnalysis/DataFormats/interface/DiBoson.h"
#include "VVXAnalysis/DataFormats/interface/Particle.h"
#include <time.h>

class WZAnalyzer: public EventAnalyzer, RegistrableAnalysis<WZAnalyzer>{

public:

  //, const std::string& filename, const double& lumi = 1., const double& externalXSection = -1., bool doBasicPlots = false

 WZAnalyzer(const AnalysisConfiguration& configuration)
   : EventAnalyzer(*(new Selector<WZAnalyzer>(*this)), 
		   configuration){
    //theHistograms.profile(genCategory);
  }

  virtual ~WZAnalyzer(){}

  virtual void analyze();

  virtual void begin();

  virtual void end(TFile &);
  
  virtual Int_t cut();

  virtual void GenAnalysis(ZZtype &, Particle &, Particle &);

  virtual void RecoAnalysis(DiBosonLepton &, Particle &, Particle &);

  virtual void GenRecoAnalysis(const ZZtype, const Particle, const Particle, const DiBosonLepton, const Particle, const Particle);


 private:
  Int_t eventGen;
  Int_t eventReco;
  Int_t eventGenReco;
  Int_t eventSample;
  Int_t recoAfterCut;
  Int_t recoJetless2;
  Int_t recoZlempty;

  Int_t gen3e;
  Int_t gen3m;
  Int_t gen2e1m;
  Int_t gen2m1e;
  Int_t reco3e;
  Int_t reco3m;
  Int_t reco2e1m;
  Int_t reco2m1e;
  
  Int_t counter1;
  Int_t counter2;
  Int_t counter3;
  Int_t counter4;

  Float_t begintime;
  Float_t endtime;
  
  friend class Selector<WZAnalyzer>;
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

  void printTime(float btime, float etime){
    cout << "\nExecution time: " << (int)((etime - btime)/3600) << " h " << (((int)(etime - btime)%3600)/60) << " m " << etime - btime - (int)((etime - btime)/3600)*3600 - (((int)(etime - btime)%3600)/60)*60 << " s." << endl;
  }

};
#endif

