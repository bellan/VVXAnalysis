#include "VVXAnalysis/Commons/interface/StringTools.h"
#include <boost/algorithm/string.hpp>
#include <boost/format.hpp>

bool strtool::sortEvents(std::string i, std::string j){

  std::vector<std::string> strs1;
  std::vector<std::string> strs2;

  boost::split(strs1, i, boost::is_any_of(":"));
  boost::split(strs2, j, boost::is_any_of(":"));

  int run1  = std::stoi(strs1[0]);
  int run2  = std::stoi(strs2[0]);
 
  if (run1==run2){
    int lumi1  = std::stoi(strs1[1]);
    int lumi2  = std::stoi(strs2[1]);
    if(lumi1==lumi2){
      unsigned long ev1  = std::stoul(strs1[2]);
      unsigned long ev2  = std::stoul(strs2[2]);
      return ev1<ev2;
    }
    else return lumi1<lumi2;
  }
  else  return run1<run2;
}


std::string strtool::sRound(float val,std::string digit)
{
  //  return  str( boost::format("%.2f") % val );
  return  str( boost::format("%"+digit+"f") % val );
}
