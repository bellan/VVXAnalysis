#ifndef ZZSAnalyzer_h
#define ZZSAnalyzer_h

/** \class ZZSAnalyzer
 *  Concrete class for ZZS analysis
 *
 *  $Date: 2013/03/15 13:37:31 $
 *  $Revision: 1.4 $
 *  \author R. Bellan - UNITO <riccardo.bellan@cern.ch>
 */


#include "EventAnalyzer.h"
#include "RegistrableAnalysis.h"
#include "VVXAnalysis/Commons/interface/Constants.h"
#include <vector>
#include <boost/assign/std/vector.hpp> // for 'operator+=()'
#include <boost/assert.hpp>
using namespace boost::assign;

class ZZSAnalyzer: public EventAnalyzer, RegistrableAnalysis<ZZSAnalyzer>{

public:
 ZZSAnalyzer(const AnalysisConfiguration& configuration)
   : EventAnalyzer(*(new Selector<ZZSAnalyzer>(*this)),
		   configuration){}
  
  virtual ~ZZSAnalyzer(){}

  virtual void analyze();
  virtual void end(TFile &);
  /* virtual Int_t cut(); */

  /* void ZZplots(int id =-1); */
  int sig;
  int bkg;
  int e;
  Long64_t nentries;
  int tot_gen;
  int sig_def;
  /* int sig_100 = 0; */
  /* int bkg_100= 0; */
  /* int sig_200 = 0; */
  /* int bkg_200 = 0; */

  /* int sig_250 = 0; */
  /* int bkg_250 = 0; */
  /* int sig_300 = 0; */
  /* int bkg_300 = 0; */
  /* int sig_350 = 0; */
  /* int bkg_350 = 0; */
  /* int sig_400 = 0; */
  /* int bkg_400 = 0; */
  /* int sig_500 = 0; */
  /* int bkg_500 = 0; */
  /* int sig_600 = 0; */
  /* int bkg_600 = 0;  */
  /* int sig_800 = 0; */
  /* int bkg_800 = 0; */

  float m4L_r;
  float m4L_g;
  float m4L_g_tot;

  std::vector<double> bins;
  //std::vector<double> ybins;

 private:
  friend class Selector<ZZSAnalyzer>;
  template< class PAR >
    bool ZBosonDefinition(phys::Boson<PAR> cand) const{
    return fabs(cand.p4().M() - phys::ZMASS) < 20;
  }
  template< class PAR >
    bool WBosonDefinition(phys::Boson<PAR> cand) const{
    return fabs(cand.p4().M() - phys::WMASS) < 40;
  }
};
#endif

