#ifndef VVXAnalysis_Commons_Utilities_h
#define VVXAnalysis_Commons_Utilities_h

#include "VVXAnalysis/DataFormats/interface/Particle.h"
#include "VVXAnalysis/Commons/interface/Constants.h"

#include <vector>

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
  
  // ----- Minimum of |m - mZ| and |m - mW|
  inline double minDM(const double& m, const double& m1 = phys::WMASS, const double& m2 = phys::ZMASS){
    return std::min( fabs(m - m1), fabs(m - m2) );
  }
}


// ----- Concatenating two vectors
template <typename T> std::vector<T> concatenate(std::vector<T> &a, std::vector<T> &b) {
  std::vector<T> ab = a;
    ab.insert(ab.end(), b.begin(), b.end());
    return ab;
}


template <class T, class C>
typename C::const_iterator closestDeltaR(const T& p, const C& container){
  return std::min_element(container.cbegin(), container.cend(),
			  [p](typename C::const_reference a, typename C::const_reference b){
			    return physmath::deltaR(p, a) < physmath::deltaR(p, b);
			  });
}

template <class T, class C>
typename C::const_iterator closestDeltaR_p(const T& p, const C& container){  // same as previous, but for containers of pointers
  return std::min_element(container.cbegin(), container.cend(),
			  [p](typename C::const_reference a, typename C::const_reference b){
			    return physmath::deltaR(p, *a) < physmath::deltaR(p, *b);
			  });
}


template <class F, class T>
std::vector<std::pair<const F*, const T*>> matchDeltaR(const std::vector<F>& vFrom, const std::vector<T>& vTo, const double& tol = 0.2){
  // Performs associations in deltaR of objects from vFrom to vTo
  // Ensures that no more than one object in vFrom is associated to the same object in vTo
  // Precedence is given following the order in which they appear
  // If a candidate match for an object in vFrom exceeds the dR threshold, it is discarded and no object from vTo is associated to it
  
  std::vector<std::pair<const F*, const T*>> out;
  out.reserve(vFrom.size());
  
  typedef typename std::vector< typename std::vector<T>::const_iterator > indices_container;  // "vector of iterators to vector<T> elements"
  indices_container indicesTo;
  for(auto it = vTo.cbegin(); it != vTo.cend(); ++it) indicesTo.push_back(it);
  
  for(const F& from : vFrom){
    typename indices_container::const_iterator bestTo = closestDeltaR_p(from, indicesTo);
    
    if(bestTo != indicesTo.end() && physmath::deltaR(from, **bestTo) < tol){
      out.push_back( std::make_pair(&from, &(**bestTo)) );
      indicesTo.erase(bestTo);
    }
    else
      out.push_back( std::make_pair(&from, nullptr  ) ); // If no corresponding "to" is found, pair::second is left as a nullptr
  }
  
  return out;
}


#endif
