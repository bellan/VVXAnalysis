#include <TCanvas.h>
#include <TFile.h>
#include <TH1F.h>
#include <TLegend.h>

#include <boost/foreach.hpp>
#define foreach BOOST_FOREACH

#define TEST_MODE

using namespace std;

void Draw(){

  vector<string> samples = {"DYJetsToLL_M50", "TTTo2L2Nu", "WZGTo2L2jG", "ZZGTo2L2jG"};
  vector<string> variables = {"mj (MATCHED AK8)", "mjj (MATCHED AK4)"};
  vector<int> colors = {kBlue, kRed, kGreen, kBlack}; // Define the line colors for each sample

  for (int j = 0; j < variables.size(); j++) {

    string variable = variables[j];
    TCanvas* canvas = new TCanvas(variable.c_str(), variable.c_str());
    canvas->cd();

    for (int i = 0; i < samples.size(); i++) {

      string sample = samples[i];
      string filename = "~/VVXAnalysis/TreeAnalysis/results/2017/ZProva_SR2P/" + sample + ".root";
      TFile* file = TFile::Open(filename.c_str());

      TH1F* h = (TH1F*)file->Get(variable.c_str());
      //h->SetLineColor(colors[i]); // Set the line color from the index i of the sample
      h->Draw((i == 0) ? "" : "same");

      file->Close();
    }

    canvas->Update();
    canvas->SaveAs((variable + ".png").c_str());
    delete canvas;
  }

}
