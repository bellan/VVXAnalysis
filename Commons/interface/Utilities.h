#ifndef VVXAnalysis_Commons_Utilities_h
#define VVXAnalysis_Commons_Utilities_h

#include "VVXAnalysis/DataFormats/interface/Particle.h"
#include "VVXAnalysis/DataFormats/interface/Photon.h"

namespace physmath{
 
  // ----- Delta Phi
  inline double deltaPhi(const double& phi1, const double &phi2) { 
    double result = phi1 - phi2;
    while (result > M_PI) result -= 2*M_PI;
    while (result <= -M_PI) result += 2*M_PI;
    return result;
  }

  inline double deltaPhi (const TLorentzVector &a, const TLorentzVector &b) {
    return deltaPhi(a.Phi(),b.Phi());
  }
  
  template<typename T1, typename T2> inline double deltaPhi(const T1& p1, const T2& p2){
  	return deltaPhi(p1.phi(), p2.phi());
  }

  
  // ----- Delta R
  template<typename T1, typename T2> double deltaR(const T1& p1, const T2& p2){
    return sqrt( deltaPhi(p1.phi(),p2.phi())*deltaPhi(p1.phi(),p2.phi()) +
		 (p1.eta()-p2.eta())*(p1.eta()-p2.eta()) );
  }
	
  inline double deltaR(const TLorentzVector &a, const TLorentzVector &b){
    //Overloads the template deltaR<T1,T2>(...). Created for deltaRComparator in "Comparators.h"
    double dPhi = deltaPhi(a.Phi(), b.Phi());
    return sqrt( dPhi*dPhi + (a.Eta() - b.Eta())*(a.Eta() - b.Eta()) );
  }
  
  
  // ----- Check if particles are the same
  template<typename T> bool isAlmostEqual(const T& a, const T& b, const double &tollerance = 0.0001){
    if      (a != 0) return fabs(a-b)/a < tollerance;
    else if (b != 0) return fabs(a-b)/b < tollerance;
    else             return fabs(a-b)   < tollerance*1e-5;
  }
  
  enum MatchingType{ratio,DR};
  template<typename T1, typename T2> bool isAlmostEqual(const T1& a, const T2& b, MatchingType mtype, const double &tollerance = 0.01){

    if(mtype == MatchingType::DR) return deltaR(a,b) < tollerance;
    else return false;
  }
  
  template<typename T> bool isTheSame(const T& p1, const T& p2){
    return abs(p1.pt()- p2.pt()) < 0.01;
  }
  
  
  // ----- Transverse Mass
  template<typename T> double mT(const T& p1, const T& p2){
    return sqrt( 2*p1.pt()*p2.pt()*(1-TMath::Cos(deltaPhi(p1.phi(), p2.phi()))) );
  }

  template<typename T> double mT(const T& p1, const T& p2, const T& p3){
    return sqrt(mT(p1, p2)*mT(p1, p2) + mT(p1, p3)*mT(p1, p3) + mT(p2, p3)*mT(p2, p3));
  }
  
}


// ----- Concatenating two vectors
template <typename T> std::vector<T> concatenate(std::vector<T> &a, std::vector<T> &b) {
  std::vector<T> ab = a;
    ab.insert(ab.end(), b.begin(), b.end());
    return ab;
}


// Test how many cuts in the ID does the photon pass. See https://twiki.cern.ch/twiki/bin/view/CMS/CutBasedPhotonIdentificationRun2#Working_points_for_94X_and_later
static int nCutsIDLoose(const phys::Photon& ph){
  bool isBarrel = fabs(ph.eta()) < 1.4442;
  auto pho_pt = ph.pt();
  // bools can be implicitly cast to ints (false -> 0, true -> 1)
  int HoverE = ph.HoverE()          < (isBarrel ? 0.04596 : 0.0590);
  int sIeIe = ph.sigmaIetaIeta()    < (isBarrel ? 0.0106  : 0.0272);
  int chIso = ph.chargedIsolation() < (isBarrel ? 1.694   : 2.089);
  int neIso = ph.neutralHadronIsolation() < (isBarrel ? 
              24.032 + 0.01512*pho_pt + 2.259e-05*pho_pt*pho_pt : 
              19.722 + 0.0117 *pho_pt + 2.3e-05  *pho_pt*pho_pt);
  int phIso = ph.photonIsolation() < (isBarrel ? 
              2.876  + 0.004017*pho_pt :
              4.162  + 0.0037*pho_pt);
  return HoverE + sIeIe + chIso + neIso + phIso;
}

static int nCutsIDMedium(const phys::Photon& ph){
  bool isBarrel = fabs(ph.eta()) < 1.4442;
  auto pho_pt = ph.pt();
  // bools can be implicitly cast to ints (false -> 0, true -> 1)
  int HoverE = ph.HoverE()          < (isBarrel ? 0.02197 : 0.0326);
  int sIeIe = ph.sigmaIetaIeta()    < (isBarrel ? 0.01015 : 0.0272);
  int chIso = ph.chargedIsolation() < (isBarrel ? 1.141   : 1.051);
  int neIso = ph.neutralHadronIsolation() < (isBarrel ? 
              1.189 + 0.01512*pho_pt + 2.259e-05*pho_pt*pho_pt : 
              2.718 + 0.0117 *pho_pt + 2.3e-05  *pho_pt*pho_pt);
  int phIso = ph.photonIsolation() < (isBarrel ? 
              2.08  + 0.004017*pho_pt :
              3.867 + 0.0037  *pho_pt);
  return HoverE + sIeIe + chIso + neIso + phIso;
}


#endif
