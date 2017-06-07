#include "VVXAnalysis/TreeAnalysis/interface/AnalysisFactory.h"
#include "VVXAnalysis/Commons/interface/Colours.h"
#include "VVXAnalysis/Commons/interface/RegionTypes.h"

#include <iostream>
#include <boost/lexical_cast.hpp>

using std::cout;
using std::endl;
using namespace colour;

int main (int argc, char ** argv){

  if(argc < 4){ 
    cout<<Important("Too few arguments")<<endl;
    cout<<"vvxAnalysis <input file name> <output filename> <luminosity (for MC)> <external cross section (for MC)>"<<endl;
    return 1;
  }
  
  AnalysisConfiguration analysisConfig;
  analysisConfig.addParameter("analysisName"    , std::string(argv[1]));
  analysisConfig.addParameter("region"          , phys::regionType(std::string(argv[2])));
  analysisConfig.addParameter("filename"        , std::string(argv[3]));
  analysisConfig.addParameter("outputfile"      , std::string(argv[4]));
  analysisConfig.addParameter("lumi"            , atof(argv[5]));
  analysisConfig.addParameter("externalXSection", atof(argv[6]));
  analysisConfig.addParameter("maxNumEvents"    , atoi(argv[7]));
  analysisConfig.addParameter("doBasicPlots"    , false);
  analysisConfig.addParameter("test"    , phys::RegionTypes::CR3P1F);
  

  if(atof(argv[8]))   analysisConfig.addParameter("doSF"    , true);
  else   analysisConfig.addParameter("doSF"    , false);

  cout<<Yellow("Analyzing "+analysisConfig.getParameter<std::string>("filename")+" ... please wait... ")<<endl ;

  EventAnalyzer *analysis = AnalysisFactory::get()->createAnalysis(analysisConfig);

  analysis->loop(analysisConfig.getParameter<std::string>("outputfile"));

  cout<<"Output saved in --> "<<Green(analysisConfig.getParameter<std::string>("outputfile"))<<endl;
  cout<<"\nAnalysis status: "<<OK("DONE")<<"\n"<<endl;

  return 0;
}


