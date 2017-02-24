#ifndef VVXAnalysis_Commons_StringTools_H
#define VVXAnalysis_Commons_StringTools_H
#include <string>
#include <iostream> 

namespace strtool{
  bool sortEvents(std::string i, std::string j);
  
  std::string sRound(float val, std::string digit = ".2");
 
}
#endif
