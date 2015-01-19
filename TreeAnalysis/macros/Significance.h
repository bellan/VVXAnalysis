#include <iostream>

enum Significance {SoverB, SoverSB, Punzi, PunziSimple, PunziImproved, Stop};

double punziSmin(const double &B, const double &a = -1, const double &b = -1){
  if(a <= 0  || b <= 0){
    std::cout << "Error! Punzi's formula not correctly implemented! " << B << " " << a << " " << b << std::endl;
    return -1;
  }
  return 0.5*(b*b + 2*a*sqrt(B) + b*sqrt(b*b + 4*a*sqrt(B)+4*B)); // here S means signal efficency!!!
}

double significance(Significance type, const double &S, const double &B, const double &a = -1, const double &b = -1){
  
  switch (type){
  case SoverB:
    return S/sqrt(B);

  case SoverSB:
    return S/sqrt(S+B);

  case Stop:
    if(a <= 0  || b <= 0){
      std::cout << "Error! Stop formula not correctly implemented! " << S << " " << a << " " << b << std::endl;
      return -1;
    }
    return S/(S + B + a*a*S*S + b*b*B*B + 3./2); // here a and b means systematic uncertainty!
    
  case Punzi:
    if(a <= 0  || b <= 0){
      std::cout << "Error! Punzi's formula not correctly implemented! " << S << " " << a << " " << b << std::endl;
      return -1;
    }
    return S/punziSmin(B, a , b); // here S means signal efficency!!!

  case PunziSimple:
    if(a <= 0){
      std::cout << "Error! Punzi's formula not correctly implemented! " << S << " " << a << std::endl;
      return -1;
    }
    return S/(a/2 + sqrt(B)); // here S means signal efficency!!!
    
  case PunziImproved:
    if(a <= 0  || b <= 0){
      std::cout << "Error! Punzi's formula not correctly implemented! " << S << " " << a << " " << b << std::endl;
      return -1;
    }
    return S/(a*a/8 + 9*b*b/13 + a*sqrt(B) + b*sqrt(b*b + 4*a*sqrt(B)+4*B)/2); // here S means signal efficency!!!

       
  default:
    return -1;
    }
  
}
