template<typename T> bool isAlmostEqual(const T& a, const T& b, const double &tollerance = 0.0001){
  if      (a != 0) return abs(a-b)/a < tollerance;
  else if (b != 0) return abs(a-b)/b < tollerance;
  else             return abs(a-b)   < tollerance*1e-5;
  
}

namespace phys{
  template<typename T1, typename T2> double deltaR(const T1& p1, const T2& p2){
    return sqrt( (p1.phi()-p2.phi())*(p1.phi()-p2.phi()) +
		 (p1.eta()-p2.eta())*(p1.eta()-p2.eta()) );
  }
}
