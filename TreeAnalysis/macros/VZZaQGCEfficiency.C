#include <TCanvas.h>
#include <TFile.h>
#include <TH1F.h>
#include <TLegend.h>

#include <boost/foreach.hpp>
#define foreach BOOST_FOREACH

#define TEST_MODE

using namespace std;

void WZZEfficiencyAnalysis(){

  TFile* result = TFile::Open("~/VVXAnalysis/TreeAnalysis/results/2018/WZZAnalyzer_MC/WZZ.root");
  result->Close("R");
}
