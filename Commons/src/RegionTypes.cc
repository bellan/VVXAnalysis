#include "VVXAnalysis/Commons/interface/RegionTypes.h"

#include <iostream>

phys::RegionTypes phys::regionType(const std::string& input){
  if      (input == "SR")     return SR;
  else if (input == "CR")     return CR;
  else if (input == "CR2P2F") return CR2P2F;
  else if (input == "CR3P1F") return CR3P1F;
  else if (input == "CR2P2F_HZZ") return CR2P2F_HZZ;
  else if (input == "CR3P1F_HZZ") return CR3P1F_HZZ;
  else{
    std::cout << "Unknown region: " << input << std::endl;
    abort();
  }
}
