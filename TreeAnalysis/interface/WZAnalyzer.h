#ifndef WZAnalyzer_h
#define WZAnalyzer_h

/** \class WZAnalyzer
 *  Concrete class for WZ analysis
 *
 *  $Date: 2021/06/17 $
 *  $Revision: 0.5 $
 *  \author E. Racca - UNITO <eleonora.racca@cern.ch>
 */


#include "EventAnalyzer.h"
#include "RegistrableAnalysis.h"
#include "VVXAnalysis/DataFormats/interface/TypeDefs.h"
#include "VVXAnalysis/Commons/interface/Constants.h"
#include "VVXAnalysis/DataFormats/interface/Boson.h"
#include "VVXAnalysis/DataFormats/interface/DiBoson.h"
#include "VVXAnalysis/DataFormats/interface/Particle.h"
#include "VVXAnalysis/TreeAnalysis/interface/VVjjHelper.h"

using namespace phys;
using namespace std;

class WZAnalyzer: public EventAnalyzer, RegistrableAnalysis<WZAnalyzer>{

public:
  
 WZAnalyzer(const AnalysisConfiguration& configuration)
   : EventAnalyzer(*(new Selector<WZAnalyzer>(*this)), 
		   configuration){
    helper_ = new VVjjHelper(&theHistograms);
  }
  
  virtual ~WZAnalyzer(){
    delete helper_;
  }
  
  virtual void analyze();
  
  virtual void begin();
  
  virtual void end(TFile &);
  
  virtual Int_t cut();
  

private:  
  VVjjHelper* helper_;
  
  Int_t eventGen;
  Int_t eventReco;
  Int_t eventGenReco;
  Int_t eventSample;
  Int_t eventGenaftercut;
  Int_t eventRecoaftercut;
  Int_t eventGenNOReco;
  Int_t eventRecoNOGen;
  
  Float_t weightGen;
  Float_t weightReco;
  Float_t weightGenReco;
  Float_t weightGenaftercut;
  Float_t weightRecoaftercut;
  Float_t weightGenNOReco;
  Float_t weightRecoNOGen;

  Float_t begintime;
  Float_t endtime;  

  void GenAnalysis(DiBosonParticle &, Particle &, Particle &);
  void RecoAnalysis(DiBosonLepton &, Particle &, Particle &);
  void GenRecoAnalysis(const DiBosonParticle, const Particle, const Particle, const DiBosonLepton, const Particle, const Particle);
  
  void BuildingWZ();
  void Article();
  
  
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

};
#endif
