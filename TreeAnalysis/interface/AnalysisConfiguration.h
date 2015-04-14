#ifndef VVXAnalysis_TreeAnalysis_AnalysisConfiguration_h
#define VVXAnalysis_TreeAnalysis_AnalysisConfiguration_h

/** \class AnalysisConfiguration
 *  Fold the configuration of the analysis in a unique object, so new additions are more easier to be propagated to the concrete classes
 *
 *  $Date: 2013/03/15 13:37:31 $
 *  $Revision: 1.4 $
 *  \author R. Bellan - UNITO <riccardo.bellan@cern.ch>
 */

#include <iostream>
#include <map>
#include <boost/any.hpp>

class AnalysisConfiguration{
 public:

  typedef std::map<std::string, boost::any> Parameters;

 public:
  
  AnalysisConfiguration(){}

  void addParameter(const std::string& name, const boost::any & value){
    parameters_[name] = value;
  }
  
  template<typename T>
    T getParameter(const std::string& name)const {
    auto search = parameters_.find(name);
    if(search != parameters_.end()) {
      return boost::any_cast<T>(search->second);
    }
    else{
      std::cout << "AnalysisConfiguration: parameter " << name << " not found" << std::endl;
      abort();
    }
  }

 private:
  Parameters parameters_;
};


#endif
