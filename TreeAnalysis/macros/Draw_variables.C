#include <TCanvas.h>
#include <TFile.h>
#include <TH1F.h>
#include <TLegend.h>

#include <boost/foreach.hpp>
#define foreach BOOST_FOREACH

#define TEST_MODE

using namespace std;

void Draw_variables(){

  vector<string> methods={"mZ","mW","mWZ","maxVPt","minTotPt"};
  for(auto suffix:methods)
  {
  vector<string> samples = { "WZGTo2L2jG", "ZZGTo2L2jG"};
  vector<string> variables = {"mjj_(TRUEMATCHED_AK4)","mjj_(PASSED)_"+suffix,"mjj_(FAILED)_"+suffix};//"0size_jets","1.0size_selectedRECOjets","1.1size_selectedRECOjets",
  vector<int> colors = { kRed, kBlue,kGreen,kMagenta,kOrange}; // Define the line colors for each variable

  for (int j = 0; j < samples.size(); j++) {

    string sample = samples[j];
    TCanvas* canvas = new TCanvas(sample.c_str(), sample.c_str(), 200,10,600,480);
    canvas->cd();
    TLegend * legend =new TLegend(0.7,0.5,0.9,0.7);

    for (int i = 0; i < variables.size(); i++) {

      string variable = variables[i];
      string filename = "~/VVXAnalysis/TreeAnalysis/results/2017/ZProva_SR2P/" + sample + ".root";
      cout << filename;
      TFile *file = TFile::Open(filename.c_str());

      TH1F *h = (TH1F *)file->Get(variable.c_str());
      h->SetLineColor(colors[i]); // Set the line color from the index i of the variables
      h->Draw((i == 0) ? "" : "same");//
      legend->AddEntry(h,variable.c_str(),"l");//options:lpf
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
    canvas->SaveAs((sample + "mjj_Passed_vs_Failed"+ suffix +".png").c_str());
    delete canvas;
    delete legend;

  }

  }


}
