{
  char opt[4] ="kfg";

  gSystem->CompileMacro("DataToUnfold.cxx",opt);
  gSystem->CompileMacro("PurityAndStability.cxx",opt);
  gSystem->CompileMacro("ResponseMatrix.cxx",opt);

  gROOT->ProcessLine(".L ToUnfold.C++");
}
