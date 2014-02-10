#include <iostream>
#include <boost/lexical_cast.hpp>

#include "ZZWAnalyzer.h"

using std::cout;
using std::endl;

int main (int argc, char ** argv){

  cout<<"Ciao"<<endl;
  
  if(argc < 4){ 
    cout<<"Too few arguments"<<endl;
    cout<<"eventAnalyzer <input file name> <output filename> <luminosity (for MC)> <external cross section (for MC)>"<<endl;
    return 1;
  }

  // input filename
  std::string filename(argv[1]);

  float lumi  = atof(argv[3]); 
  float externalXsec = atof(argv[4]);
    
  std::cout<<"Analyzing "<<filename<<" ... please wait... "<<endl ;
  
  ZZWAnalyzer analysis(filename, lumi, externalXsec);
  analysis.loop(argv[2]);

  std::cout<<filename<<" --> "<<argv[2]<<" done!"<<endl ;

  return 0;
}


