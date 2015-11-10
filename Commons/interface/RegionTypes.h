#ifndef VVXAnalysis_Commons_RegionTypes_h
#define VVXAnalysis_Commons_RegionTypes_h

#include <string>

namespace phys{
  enum RegionTypes {SR, CR, CR2P2F, CR3P1F, CR2P2F_HZZ, CR3P1F_HZZ, MC,SR_ZZOnShell};
 
  RegionTypes regionType(const std::string& input);
  std::string regionType(RegionTypes input);

};

#endif
