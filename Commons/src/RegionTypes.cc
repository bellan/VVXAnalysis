#include "VVXAnalysis/Commons/interface/RegionTypes.h"

#include <iostream>

phys::RegionTypes phys::regionType(const std::string& input){
  if      (input == "SR")     return SR;
  else if (input == "CR")     return CR;
  else if (input == "CR2P2F") return CR2P2F;
  else if (input == "CR3P1F") return CR3P1F;
  else if (input == "MC")            return MC;
  else if (input == "SR_HZZ")     return SR_HZZ;
  else if (input == "CR_HZZ")     return CR_HZZ;
  else if (input == "CR2P2F_HZZ") return CR2P2F_HZZ;
  else if (input == "CR3P1F_HZZ") return CR3P1F_HZZ;
  else if (input == "MC_HZZ")        return MC_HZZ;
  else if (input == "SR_ZZFull")     return SR_ZZFull;
  else if (input == "CR_ZZFull")     return CR_ZZFull;
  else if (input == "CR2P2F_ZZFull") return CR2P2F_ZZFull;
  else if (input == "CR3P1F_ZZFull") return CR3P1F_ZZFull;
  else if (input == "MC_ZZFull")     return MC_ZZFull;
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
  else if (input == MC)     return "MC";
  else if (input == SR_HZZ)     return "SR_HZZ";
  else if (input == CR_HZZ)     return "CR_HZZ";
  else if (input == CR2P2F_HZZ) return "CR2P2F_HZZ";
  else if (input == CR3P1F_HZZ) return "CR3P1F_HZZ";
  else if (input == MC_HZZ)     return "MC_HZZ";  
  else if (input == SR_ZZFull)     return "SR_ZZFull";
  else if (input == CR_ZZFull)     return "CR_ZZFull";
  else if (input == CR2P2F_ZZFull) return "CR2P2F_ZZFull";
  else if (input == CR3P1F_ZZFull) return "CR3P1F_ZZFull";
  else if (input == MC_ZZFull)     return "MC_ZZFull";
  else{
    std::cout << "Unknown region: " << input << std::endl;
    abort();
  }
}
