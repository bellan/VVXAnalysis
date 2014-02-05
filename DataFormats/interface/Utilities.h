template<typename T> bool isAlmostEqual(const T& a, const T& b, const double &tollerance = 0.0001){
  if      (a != 0) return abs(a-b)/a < tollerance;
  else if (b != 0) return abs(a-b)/b < tollerance;
  else             return abs(a-b)   < tollerance*1e-5;
  
}
