#if !defined(__CINT__) || defined(__MAKECINT__)
#include <Riostream.h>
#include <TSystem.h>
#endif


void compilemyclass(TString myopt="fast"){
  char opt[4];
  if(myopt.Contains("force")){
    sprintf(opt, "kfg");
  }
  else {
    sprintf(opt ,"kg");
  }
  
  gSystem->CompileMacro("ResponseMatrix.cxx",opt);
  gSystem->CompileMacro("DataToUnfold.cxx",opt);
  gSystem->CompileMacro("PurityAndStability.cxx",opt);
  

}

