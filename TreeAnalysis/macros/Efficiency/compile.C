#include "TSystem.h"

void compile(){
	gSystem->SetBuildDir("bin");
  gSystem->CompileMacro("Efficiency.C","kg-");
  gSystem->CompileMacro("EffCalls.C","kg-");
}
