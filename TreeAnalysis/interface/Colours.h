#ifndef VVXAnalysis_TreeAnalysis_Colours_H
#define VVXAnalysis_TreeAnalysis_Colours_H

/** \class Colours
 *  Not exactly a data format... but text data format
 *
 *  $Date: 2013/03/15 13:37:31 $
 *  $Revision: 1.3 $
 *  \author R. Bellan - UCSB <riccardo.bellan@cern.ch>
 */

#include <sstream> 

class Colour{
 public:
 Colour(const std::string& trailer)
    : trailer_(trailer)
    , tail_("\033[00m"){}

  template<typename T>
  std::string operator()(T i){

    std::ostringstream s1;
    s1 << trailer_ << i << tail_;
    return s1.str();
  }
  
 protected:
  std::string trailer_;
  std::string tail_;

};

#endif
