
#include <TFile.h>
#include <TROOT.h>
#include <TStyle.h>
#include <TH1.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <iostream>
#include <string>


void SignalHisto_data(string finalstate = "4m", string unc = "Up")
{
  gROOT->Reset();  
  gROOT->SetStyle("Plain");   
  gStyle->SetOptStat(0);
  
  TFile *output = new TFile("Jets_test/DataToUnfold_JES.root", "UPDATE");
  //Data (ZZRecoAnalyzer not yet done FIXME)
  TFile *data = new TFile("../results/ZZRecoAnalyzer_SR/data.root");
  //reducible background
  TFile *red = new TFile("../results/ZZRecoAnalyzer_CR/data.root"); 
  //TFile *red = new TFile("../results/ZZRecoAnalyzer_CR_150623_data/data.root"); 
  //irreducible backgrounds
  TFile *ttZ = new TFile("../results/ZZRecoAnalyzer_SR/TTZJets.root");
  TFile *ttWW = new TFile("../results/ZZRecoAnalyzer_SR/TTWWJets.root");
  TFile *WWZ = new TFile("../results/ZZRecoAnalyzer_SR/WWZJets.root");
  
  TH1 * h_data;
  TH1 * h_totdata;
  TH1 * h_red;
  TH1 * h_ttZ;
  TH1 * h_ttWW;
  TH1 * h_WWZ;
  
  string histoName = "ZZTo" + finalstate + "_Jets_JESData"+unc+"Smear_01"; 
  // string histoName_red = "ZZTo" + finalstate + "_Jets_01";
  string histoMCName =  "ZZTo" + finalstate + "_Jets_JERCentralSmear_01";

  h_totdata = (TH1*) data->Get(histoName.c_str()); 
  h_data = (TH1*)h_totdata->Clone("h_totdata");
  h_red = (TH1*) red->Get(histoName.c_str()); 
  h_ttZ = (TH1*) ttZ->Get(histoMCName.c_str()); 
  h_ttWW = (TH1*) ttWW->Get(histoMCName.c_str()); 
  h_WWZ = (TH1*) WWZ->Get(histoMCName.c_str()); 

  // float err_tot = 0;
  // float err_sist = 0; 
  // float err_stat = 0;
  // float err_red = 0;
  // float err_irr = 0;
  float dataminusbkg = 0;

  std::cout << " =======================================================================================================================================================================" << std::endl;
    std::cout << "                                          " << finalstate.c_str() << " " << unc.c_str() << std::endl;
 std::cout << " ========================================================================================================================================================================" << std::endl;
  for(int i =1; i<5; i++){
    dataminusbkg = 0;
    dataminusbkg = h_totdata->GetBinContent(i)- h_red->GetBinContent(i)- h_ttZ->GetBinContent(i)-h_ttWW->GetBinContent(i)-h_WWZ->GetBinContent(i);
    if(dataminusbkg>0.) h_data-> SetBinContent(i,dataminusbkg);
    else h_data-> SetBinContent(i,0.);
    
    //adding background systematic uncertainties to the statistical uncertainty FIXME
    // err_tot =0;
    // err_red =0;
    // err_irr =0;
    // err_sist = 0; 
    // err_stat = 0;
    // err_red = h_red->GetBinError(i);
    // err_irr =sqrt(h_ttZ->GetBinError(i)*h_ttZ->GetBinError(i)+h_ttWW->GetBinError(i)*h_ttWW->GetBinError(i)+h_WWZ->GetBinError(i)*h_WWZ->GetBinError(i));
    // err_sist = sqrt(err_red*err_red+err_irr*err_irr);
    // err_stat = h_totdata->GetBinError(i);
    // err_tot = sqrt(err_stat*err_stat+err_sist*err_sist);
    // h_data->SetBinError(i,err_tot);

    //some information
  
    std::cout << "bin " << i << " data = "<< h_totdata->GetBinContent(i) << " +- " << h_totdata->GetBinError(i) << " red = " << h_red->GetBinContent(i) << " +- " << h_red->GetBinError(i) <<   " ttZ = " << h_ttZ->GetBinContent(i) << " +- " << h_ttZ->GetBinError(i) <<  " ttWW = " << h_ttWW->GetBinContent(i) << " +- " << h_ttWW->GetBinError(i) <<   " WWZ = " << h_WWZ->GetBinContent(i)   << " +- " << h_WWZ->GetBinError(i) << " data-bkg = " << dataminusbkg <<" +- " << h_data->GetBinError(i) << " " << h_data->GetBinContent(i)<< std::endl;
  }
  
  TLegend *leg = new TLegend(0.65,0.65,0.6,0.85);
  leg->SetFillColor(kWhite);
  leg->SetBorderSize(0);
  leg->SetTextSize(0.03); 
 
  TCanvas *c = new TCanvas ("c","c");
  h_totdata->Draw(); 
  h_data->SetLineColor(1);
  h_data->Draw("SAME");
  h_data->SetLineColor(2);

   leg->AddEntry(h_totdata,"full dataset","l"); 
  leg->Draw(); 
  leg->AddEntry(h_data,"data-background","l"); 
  leg->Draw("SAME");
   
  string dataName = "DataminusBkg_Jets_ZZTo"+finalstate+"_JES"+unc;
  string TotdataName = "TotData_Jets_ZZTo"+finalstate+"_JES"+unc;
 
  output->cd();   
  h_data->Write(dataName.c_str());
  h_totdata->Write(TotdataName.c_str());
  output->Close();
 
}

void MakeAllFinalStates(){
  SignalHisto_data("4m","Up"); 
  SignalHisto_data("4m","Down");
  SignalHisto_data("4e","Up"); 
  SignalHisto_data("4e","Down");
  SignalHisto_data("2e2m","Up"); 
  SignalHisto_data("2e2m","Down");
 }

void PlotVariations(string finalstate = "4e"){

  gROOT->Reset();  
  gROOT->SetStyle("Plain");   
  gStyle->SetOptStat(0);
  gStyle->SetPaintTextFormat("5.3f");

  TFile *file_JES = new TFile("DataToUnfold_JES.root");
  TFile *file = new TFile("DataToUnfold.root");
  string data =  "DataminusBkg_Jets_ZZTo"+finalstate;
  string data_up =  "DataminusBkg_Jets_ZZTo"+finalstate+"_JESUp";
  string data_down =  "DataminusBkg_Jets_ZZTo"+finalstate+"_JESDown";

  //string title = "Response Matrix N jets - "+finalstate+"final state (dataset: " + dataset+" unc: "+ unc+")";
  string xAxis = "N jets ("+ finalstate+" final state)";
  string yAxis = "Events";


  TH1 * h_data = (TH1*) file->Get(data.c_str());
  TH1 * h_data_up = (TH1*) file_JES->Get(data_up.c_str());
  TH1 * h_data_down = (TH1*) file_JES->Get(data_down.c_str());
  TCanvas *c = new TCanvas("c","c"); 
 
  // double max = h_data->GetBinContent(1)+10;

  h_data_down->SetTitle("");
  h_data_down->GetXaxis()->SetTitle(xAxis.c_str());
  h_data_down->GetYaxis()->SetTitle(yAxis.c_str());
  //  h_data_down->SetMaximum(max);
  h_data->SetMarkerColor(1);
  h_data->SetMarkerStyle(8);
  h_data->SetLineColor(1);
  //h_data->Draw("E");

  h_data_down->Draw("");
  h_data_down->SetLineColor(2);
  h_data_down->SetLineWidth(2);
  h_data_up->Draw("SAME");
  h_data_up->SetLineColor(2); 
  h_data_up->SetLineWidth(2);
  h_data->Draw("E SAME");
  TLegend *leg = new TLegend(0.65,0.45,0.55,0.55); 
  leg->SetFillColor(kWhite);
  leg->SetBorderSize(0);
  leg->SetTextSize(0.03);
  leg->SetTextFont(42); 
 
  leg->AddEntry(h_data,"Signal yield from data","lep"); 
  leg->Draw(); 
  leg->AddEntry(h_data_up,"Variations due to #pm d#sigma_{JES} ","l"); 
  leg->Draw("SAME");

  string pdf ="DataMinusBkg_Jets_JES_ZZTo" + finalstate +".pdf";
  string png ="DataMinusBkg_Jets_JES_ZZTo" + finalstate +".png";
 
  c->Print(png.c_str());
}
