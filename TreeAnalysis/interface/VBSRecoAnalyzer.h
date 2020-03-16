#ifndef VBSRecoAnalyzer_h
#define VBSRecoAnalyzer_h

/** \class VBSRecoAnalyzer
 *  Concrete class for VBS analysis
 *
 *  $Date: 2013/03/15 13:37:31 $
 *  $Revision: 1.4 $
 *  \author R. Bellan - UNITO <riccardo.bellan@cern.ch>
 */

#include "TRandom.h"
#include "EventAnalyzer.h"
#include "RegistrableAnalysis.h"
#include "VVXAnalysis/Commons/interface/Constants.h"
class VBSRecoAnalyzer: public EventAnalyzer, RegistrableAnalysis<VBSRecoAnalyzer>{

public:
 VBSRecoAnalyzer(const AnalysisConfiguration& configuration)
   : EventAnalyzer(*(new Selector<VBSRecoAnalyzer>(*this)),
		   configuration){}
  
  virtual ~VBSRecoAnalyzer(){}

  void FillHistosJets(std::string decay, float Wh, std::vector<phys::Jet> *jetsVec,std::string type);
  virtual void analyze();
  virtual void end( TFile &);
  virtual void begin();
  


  /* // Jets obtained by gaussian JER smearing */
  std::vector<phys::Jet> *UpJER_jets;
  std::vector<phys::Jet> *DownJER_jets;
  
  // Jets obtained correcting up and down for the JES uncertainty
  std::vector<phys::Jet> *UpJES_jets;
  std::vector<phys::Jet> *DownJES_jets;
  
  //  Jets obtained correcting up and down for the JES uncertainty (data distributions = no JER correction applied)
  std::vector<phys::Jet> *UpJESData_jets;
  std::vector<phys::Jet> *DownJESData_jets;
  

 private:

  friend class Selector<VBSRecoAnalyzer>;
  template< class PAR >
    bool ZBosonDefinition(phys::Boson<PAR> cand) const{
    return fabs(cand.p4().M() - phys::ZMASS) < 20;
  }
  template< class PAR >
    bool WBosonDefinition(phys::Boson<PAR> cand) const{
    return fabs(cand.p4().M() - phys::WMASS) < 40;
  }

  Int_t nJets;
  Float_t mjj;
  Float_t deltaEtajj;
  std::vector<double> Xbins_mass; 
};
#endif

