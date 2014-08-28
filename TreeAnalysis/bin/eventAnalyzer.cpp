#include "VVXAnalysis/TreeAnalysis/interface/AnalysisFactory.h"
#include "VVXAnalysis/Commons/interface/Colours.h"

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

  // Region:
  std::string region(argv[2]);

  // input filename
  std::string filename(argv[3]);

  float lumi  = atof(argv[5]); 
  float externalXsec = atof(argv[6]);
    
  std::cout<<Yellow("Analyzing "+filename+" ... please wait... ")<<endl ;
    
  std::string analysisName = argv[1];

  EventAnalyzer *analysis = AnalysisFactory::get()->createAnalysis(analysisName, region, filename, lumi, externalXsec);
  analysis->loop(argv[4]);

  cout<<"Output saved in --> "<<Green(argv[4])<<endl;
  cout<<"\nAnalysis status: "<<OK("DONE")<<"\n"<<endl;

  return 0;
}


