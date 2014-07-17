#include "VVXAnalysis/TreeAnalysis/interface/AnalysisFactory.h"
#include "VVXAnalysis/Commons/interface/Colours.h"
#include "VVXAnalysis/TreeAnalysis/interface/GenEventAnalyzer.h"

#include <iostream>
#include <string>
#include <boost/lexical_cast.hpp>

using std::cout;
using std::endl;
using namespace colour;

int main (int argc, char ** argv){

  if(argc < 4){ 
    cout<<Important("Too few arguments")<<endl;
    cout<<"MadGraph Analysis <input file name> <output filename> <luminosity (for MC)> <cross section>"<<endl;
    return 1;
  }


  // input filename
  std::string filename(argv[1]);

  float lumi  = atof(argv[3]); 
    
  std::cout<<Yellow("Analyzing "+filename+" ... please wait... ")<<endl ;
    
  float xsec = std::stof(argv[4]);

  GenEventAnalyzer *analysis = new GenEventAnalyzer(filename, lumi, xsec);
  analysis->loop(argv[2]);

  cout<<"Output saved in --> "<<Green(argv[2])<<endl;
  cout<<"\nAnalysis status: "<<OK("DONE")<<"\n"<<endl;

  return 0;


}


