#include "VVXAnalysis/TreeAnalysis/interface/ZZSAnalyzer.h"
#include "VVXAnalysis/TreeAnalysis/interface/Colours.h"

#include <iostream>
#include <boost/lexical_cast.hpp>

using std::cout;
using std::endl;
using namespace colour;

int main (int argc, char ** argv){

  if(argc < 4){ 
    cout<<Important("Too few arguments")<<endl;
    cout<<"zzsAnalysis <input file name> <output filename> <luminosity (for MC)> <external cross section (for MC)>"<<endl;
    return 1;
  }

  // input filename
  std::string filename(argv[1]);

  float lumi  = atof(argv[3]); 
  float externalXsec = atof(argv[4]);
    
  std::cout<<Red("Analyzing "+filename+" ... please wait... ")<<endl ;
  
  ZZSAnalyzer analysis(filename, lumi, externalXsec);
  analysis.loop(argv[2]);

  cout<<"Output saved in --> "<<Green(argv[2])<<endl;
  cout<<"\nAnalysis status: "<<OK("done!")<<"\n"<<endl;

  return 0;
}


