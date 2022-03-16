#include "VVXAnalysis/TreeAnalysis/interface/AnalysisFactory.h"
#include "VVXAnalysis/Commons/interface/Colours.h"
#include "VVXAnalysis/DataFormats/interface/RegionTypes.h"

#include <iostream>
#include <boost/lexical_cast.hpp>

using std::cout;
using std::endl;
using namespace colour;

std::vector<phys::RegionTypes> regions_from_str(const std::string&);

int main (int argc, char ** argv){

  if(argc < 4){ 
    cout<<Important("Too few arguments")<<endl;
    cout<<"vvxAnalysis <input file name> <output filename> <luminosity (for MC)> <external cross section (for MC)>"<<endl;
    return 1;
  }
  
  AnalysisConfiguration analysisConfig;
  analysisConfig.addParameter("analysisName"    , std::string(argv[1]));
  analysisConfig.addParameter("regions"         , regions_from_str(std::string(argv[2])));
  analysisConfig.addParameter("filename"        , std::string(argv[3]));
  analysisConfig.addParameter("outputfile"      , std::string(argv[4]));

  analysisConfig.addParameter("year"            , atoi(argv[5]));

  analysisConfig.addParameter("lumi"            , atof(argv[6]));
  analysisConfig.addParameter("externalXSection", atof(argv[7]));
  analysisConfig.addParameter("maxNumEvents"    , atoi(argv[8]));
  analysisConfig.addParameter("doBasicPlots"    , false);
  analysisConfig.addParameter("test"    , phys::RegionTypes::CR3P1F);
  

  if(atof(argv[9]))   analysisConfig.addParameter("doSF"    , true);
  else   analysisConfig.addParameter("doSF"    , false);

  if(atoi(argv[10])) analysisConfig.addParameter("blinded"  , false);
  else analysisConfig.addParameter("blinded"  , true);

  if(atoi(argv[11])) analysisConfig.addParameter("applyFRSF", false);
  else analysisConfig.addParameter("applyFRSF"  , true);

  if(atoi(argv[12])) analysisConfig.addParameter("forcePosWeight", true);
  else analysisConfig.addParameter("forcePosWeight"  , false);

  
  EventAnalyzer *analysis = AnalysisFactory::get()->createAnalysis(analysisConfig);

  analysis->loop(analysisConfig.getParameter<std::string>("outputfile"));

  cout<<"Output saved in --> "<<Green(analysisConfig.getParameter<std::string>("outputfile"))<<endl;
  cout<<"\nAnalysis status: "<<OK("DONE")<<"\n"<<endl;

  return 0;
}

std::vector<phys::RegionTypes> regions_from_str(const std::string& in){
  // Parses a string like "SR4P;CR3P1F;CR2P2F" and returns a vector of phys::RegionTypes
  std::istringstream iss(in);
  std::string item;
  std::vector<phys::RegionTypes> regions;
  
  while(std::getline(iss, item, ';'))
    regions.push_back(phys::regionType(item));
  
  return regions;
}

