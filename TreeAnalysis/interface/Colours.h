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

namespace colour{
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

  static Colour White("\033[0;07m");
  
  static Colour Red("\033[0;31m");

  static Colour Green("\033[0;32m");
  
  static Colour Yellow("\033[0;33m");
  
  static Colour Blue("\033[0;34m");

  static Colour Violet("\033[0;35m");
  
  static Colour Important("\033[0;41m");

  static Colour OK("\033[0;42m");
  
  static Colour Warning("\033[0;43m");
}
#endif
