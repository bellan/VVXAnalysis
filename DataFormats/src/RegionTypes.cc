#include "VVXAnalysis/DataFormats/interface/RegionTypes.h"

#include <iostream>

phys::RegionTypes phys::regionType(const std::string& input){
  if      (input == "SR4P")       return SR4P  ;
  else if (input == "CR3P1F")     return CR3P1F;
  else if (input == "CR2P2F")     return CR2P2F;
  else if (input == "SR4P_1L")    return SR4P_1L;
  else if (input == "SR4P_1P")    return SR4P_1P;
  else if (input == "CR4P_1F")    return CR4P_1F;

  else if (input == "SR3P")       return SR3P  ;
  else if (input == "CR110")      return CR110 ;
  else if (input == "CR101")      return CR101 ;
  else if (input == "CR011")      return CR011 ;
  else if (input == "CR100")      return CR100 ;
  else if (input == "CR001")      return CR001 ;
  else if (input == "CR010")      return CR010 ;
  else if (input == "CR000")      return CR000 ;
  else if (input == "SR3P_1L")    return SR3P_1L;
  else if (input == "SR3P_1P")    return SR3P_1P;
  else if (input == "CR3P_1F")    return CR3P_1F;
  else if (input == "CRLFR")      return CRLFR;

  else if (input == "SR2P")       return SR2P  ;
  else if (input == "SR2P_1L")    return SR2P_1L;
  else if (input == "SR2P_1P")    return SR2P_1P;
  else if (input == "CR2P_1F")    return CR2P_1F;
  else if (input == "CR")         return CR;
  else if (input == "SR_HZZ")     return SR_HZZ;
  else if (input == "CR2P2F_HZZ") return CR2P2F_HZZ;
  else if (input == "CR3P1F_HZZ") return CR3P1F_HZZ;
  else if (input == "MC_HZZ")     return MC_HZZ;
  else if (input == "MC")         return MC;

  else{
    std::cout << "Unknown region: " << input << std::endl;
    abort();
  }
}

std::string phys::regionType(phys::RegionTypes input){

  if      (input == SR4P)       return "SR4P"  ;
  else if (input == CR3P1F)     return "CR3P1F";
  else if (input == CR2P2F)     return "CR2P2F";
  else if (input == SR4P_1L)    return "SR4P_1L";
  else if (input == SR4P_1P)    return "SR4P_1P";
  else if (input == CR4P_1F)    return "CR4P_1F";

  else if (input == SR3P)       return "SR3P"  ;
  else if (input == CR110)      return "CR110" ;
  else if (input == CR101)      return "CR101" ;
  else if (input == CR011)      return "CR011" ;
  else if (input == CR100)      return "CR100" ;
  else if (input == CR001)      return "CR001" ;
  else if (input == CR010)      return "CR010" ;
  else if (input == CR000)      return "CR000" ;
  else if (input == SR3P_1L)    return "SR3P_1L";
  else if (input == SR3P_1P)    return "SR3P_1P";
  else if (input == CR3P_1F)    return "CR3P_1F";
  else if (input == CRLFR)      return "CRLFR";

  else if (input == SR2P)       return "SR2P"  ;
  else if (input == SR2P_1L)    return "SR2P_1L";
  else if (input == SR2P_1P)    return "SR2P_1P";
  else if (input == CR2P_1F)    return "CR2P_1F";
  else if (input == CR)         return "CR";
  else if (input == SR_HZZ)     return "SR_HZZ";
  else if (input == CR2P2F_HZZ) return "CR2P2F_HZZ";
  else if (input == CR3P1F_HZZ) return "CR3P1F_HZZ";
  else if (input == MC_HZZ)     return "MC_HZZ";
  else if (input == MC)         return "MC";

  else{
    std::cout << "Unknown region: " << input << std::endl;
    abort();
  }
}


phys::Channel phys::channelType(const std::string& input){
  if      (input == "ZZ")    return ZZ;
  else if (input == "ZW")    return ZW;
  else if (input == "WZ")    return ZW;
  else if (input == "ZL")    return ZL;
  else if (input == "ZV")    return ZV;
  else if (input == "UNDEF") return UNDEF;

  else{
    std::cout << "Unknown channel: " << input << std::endl;
    abort();
  }
}

std::string phys::channelType(phys::Channel input){

  if      (input == ZZ)    return "ZZ";
  else if (input == ZW)    return "ZW";
  else if (input == ZL)    return "ZL";
  else if (input == ZV)    return "ZV";
  else if (input == UNDEF) return "UNDEF";


  else{
    std::cout << "Unknown channel: " << input << std::endl;
    abort();
  }
}
