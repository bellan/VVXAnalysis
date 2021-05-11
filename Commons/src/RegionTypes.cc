#include "VVXAnalysis/Commons/interface/RegionTypes.h"

#include <iostream>

phys::RegionTypes phys::regionType(const std::string& input){
  if      (input == "SR")     return SR;
  else if (input == "CR")     return CR;
  else if (input == "CR2P2F") return CR2P2F;
  else if (input == "CR3P1F") return CR3P1F;
  else if (input == "SR_HZZ") return SR_HZZ;
  else if (input == "CR2P2F_HZZ") return CR2P2F_HZZ;
  else if (input == "CR3P1F_HZZ") return CR3P1F_HZZ;
  else if (input == "MC") return MC;
  else if (input == "MC_HZZ") return MC_HZZ;
  else if (input == "SR3L") return SR3L;
  else{
    std::cout << "Unknown region: " << input << std::endl;
    abort();
  }
}

std::string phys::regionType(phys::RegionTypes input){
  if      (input == SR)     return "SR";
  else if (input == CR)     return "CR";
  else if (input == CR2P2F) return "CR2P2F";
  else if (input == CR3P1F) return "CR3P1F";
  else if (input == SR_HZZ) return "SR_HZZ";
  else if (input == CR2P2F_HZZ) return "CR2P2F_HZZ";
  else if (input == CR3P1F_HZZ) return "CR3P1F_HZZ";
  else if (input == MC) return "MC";
  else if (input == MC_HZZ) return "MC_HZZ";
  else if (input == SR3L) return "SR3L";
  else{
    std::cout << "Unknown region: " << input << std::endl;
    abort();
  }
}
