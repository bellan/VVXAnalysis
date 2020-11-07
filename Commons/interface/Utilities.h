#ifndef VVXAnalysis_Commons_Utilities_h
#define VVXAnalysis_Commons_Utilities_h

#include "VVXAnalysis/DataFormats/interface/Particle.h"

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


#endif
