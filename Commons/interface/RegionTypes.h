#ifndef VVXAnalysis_Commons_RegionTypes_h
#define VVXAnalysis_Commons_RegionTypes_h

#include <string>

namespace phys{
  enum RegionTypes {SR, CR, CR2P2F, CR3P1F, SR_HZZ, CR_HZZ, CR2P2F_HZZ, CR3P1F_HZZ, MC, MC_HZZ};
 
  RegionTypes regionType(const std::string& input);
  std::string regionType(RegionTypes input);

};

#endif
