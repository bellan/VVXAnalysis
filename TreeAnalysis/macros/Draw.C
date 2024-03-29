#include <TCanvas.h>
#include <TFile.h>
#include <TH1F.h>
#include <TLegend.h>

#include <boost/foreach.hpp>
#define foreach BOOST_FOREACH

#define TEST_MODE

using namespace std;

void Draw(){

  vector<string> samples = { "WZGTo2L2jG", "ZZGTo2L2jG"};
  vector<string> variables = {"1.0size_selectedGENjets"};//,"1.0size_selectedRECOjets","1.1size_selectedRECOjets","3size_RECOtruematchedjets","4size_RECOmatchedjets"
  vector<int> colors = {kBlue, kRed, kGreen,kOrange}; // Define the line colors for each sample

  for (int j = 0; j < variables.size(); j++) {

    string variable = variables[j];
    TCanvas* canvas = new TCanvas(variable.c_str(), variable.c_str(), 200,10,600,480);
    canvas->cd();
    TLegend * legend =new TLegend(0.7,0.5,0.9,0.7);

    for (int i = 0; i < samples.size(); i++) {

      string sample = samples[i];
      string filename = "~/VVXAnalysis/TreeAnalysis/results/2017/ZProva_SR2P/" + sample + ".root";
      cout << filename;
      TFile *file = TFile::Open(filename.c_str());

      TH1F *h = (TH1F *)file->Get(variable.c_str());
      h->SetLineColor(colors[i]); // Set the line color from the index i of the sample
      h->Draw((i == 0) ? "" : "same");
      legend->AddEntry(h,sample.c_str(),"l");//options:lpf
      legend->Draw();

      // h->Draw("same");
      // file->Close();
      // auto legend = new TLegend(0.1, 0.7, 0.48, 0.9);
      // legend->SetHeader("Sample", "C"); // option "C" allows to center the header
      // legend->AddEntry(h,filename.c_str(), "f");
      // legend->Draw("same");
      }


//    canvas->Update();
    // canvas->BuildLegend();
    
    canvas->SaveAs((variable + ".png").c_str());
    delete canvas;
    delete legend;

  }

}
