#ifndef VVXAnalysis_Commons_RegionTypes_h
#define VVXAnalysis_Commons_RegionTypes_h

#include <string>

namespace phys{
  enum RegionTypes {SR, CR2P2F, CR3P1F};
  
  RegionTypes regionType(const std::string& input);

};

#endif
