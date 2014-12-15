#include <iostream>
#include <TH1F.h>

#include "Significance.h"

void testSignificance(){

  double S = 15;
  double B = 9; 

  std::cout << "Signal: " << S << " background: " << B << std::endl;

  std::cout << "Test S over sqrt B: " << significance(SoverB, S, B) << std::endl;

  std::cout << "Test S over sqrt S+B: " << significance(SoverSB, S, B) << std::endl;

  std::cout << "Test Punzi: " << significance(Punzi, S/100., B, 1.645, 1.645) << std::endl;
  std::cout << "Test Punzi: " << significance(Punzi, S/100., B, 3, 1.645) << std::endl;
  std::cout << "Test Punzi: " << significance(Punzi, S/100., B, 5, 1.645) << std::endl;

  std::cout << "Test PunziSimple: " << significance(PunziSimple, S/100., B, 1.645) << std::endl;

  std::cout << "Test PunziImproved: " << significance(PunziImproved, S/100., B, 1.645, 1.645) << std::endl;
  std::cout << "Test PunziImproved: " << significance(PunziImproved, S/100., B, 3    , 1.645) << std::endl;
  std::cout << "Test PunziImproved: " << significance(PunziImproved, S/100., B, 5    , 1.645) << std::endl;

  std::cout << "Test Stop: " << significance(Stop, S, B, 0.20, 0.70) << std::endl;


}
