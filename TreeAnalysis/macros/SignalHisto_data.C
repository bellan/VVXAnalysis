#include <TFile.h>
#include <TROOT.h>
#include <TStyle.h>
#include <TH1.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <iostream>
#include <string>


void SignalHisto_data(string finalstate = "4m")
{
  gROOT->Reset();  
  gROOT->SetStyle("Plain");   
  gStyle->SetOptStat(0);
  
  TFile *output = new TFile("DataToUnfold.root", "UPDATE");
  //Data (ZZRecoAnalyzer not yet done FIXME)
  TFile *data = new TFile("/afs/cern.ch/user/l/lfinco/work/VVScattering/CMSSW_5_3_11/src/VVXAnalysis/TreeAnalysis/results/ZZRecoAnalyzer_SR/data.root");
  //reducible background
  TFile *red = new TFile("/afs/cern.ch/user/l/lfinco/work/VVScattering/CMSSW_5_3_11/src/VVXAnalysis/TreeAnalysis/results/ZZRecoAnalyzer_CR/data.root"); 
  //irreducible backgrounds
  TFile *ttZ = new TFile("/afs/cern.ch/user/l/lfinco/work/VVScattering/CMSSW_5_3_11/src/VVXAnalysis/TreeAnalysis/results/ZZRecoAnalyzer_SR/TTZJets.root");
  TFile *ttWW = new TFile("/afs/cern.ch/user/l/lfinco/work/VVScattering/CMSSW_5_3_11/src/VVXAnalysis/TreeAnalysis/results/ZZRecoAnalyzer_SR/TTWWJets.root");
  TFile *WWZ = new TFile("/afs/cern.ch/user/l/lfinco/work/VVScattering/CMSSW_5_3_11/src/VVXAnalysis/TreeAnalysis/results/ZZRecoAnalyzer_SR/WWZJets.root");
  
  TH1 * h_data;
  TH1 * h_totdata;
  TH1 * h_red;
  TH1 * h_ttZ;
  TH1 * h_ttWW;
  TH1 * h_WWZ;
  
  string histoName = "ZZTo" + finalstate + "_Mass_01";

  h_totdata = (TH1*) data->Get(histoName.c_str()); 
  h_data = (TH1*)h_totdata->Clone("h_totdata");
  h_red = (TH1*) red->Get(histoName.c_str()); 
  h_ttZ = (TH1*) ttZ->Get(histoName.c_str()); 
  h_ttWW = (TH1*) ttWW->Get(histoName.c_str()); 
  h_WWZ = (TH1*) WWZ->Get(histoName.c_str()); 

  // float err_tot = 0;
  // float err_sist = 0; 
  // float err_stat = 0;
  // float err_red = 0;
  // float err_irr = 0;
  float dataminusbkg = 0;
  for(int i =1; i<9; i++){
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
    // std::cout << "bin " << i << " data = "<< h_totdata->GetBinContent(i) << " +- " << h_totdata->GetBinError(i) << " red = " << h_red->GetBinContent(i) << " +- " << h_red->GetBinError(i) <<   " ttZ = " << h_ttZ->GetBinContent(i) << " +- " << h_ttZ->GetBinError(i) <<  " ttWW = " << h_ttWW->GetBinContent(i) << " +- " << h_ttWW->GetBinError(i) <<   " WWZ = " << h_WWZ->GetBinContent(i)   << " +- " << h_WWZ->GetBinError(i) << " data-bkg = " << dataminusbkg <<" +- " << h_data->GetBinError(i) << " " << h_data->GetBinContent(i)<< std::endl;
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
   
  string dataName = "DataminusBkg_Mass_ZZTo"+finalstate;
  string TotdataName = "TotData_Mass_ZZTo"+finalstate;
 
  output->cd();   
  h_data->Write(dataName.c_str());
  h_totdata->Write(TotdataName.c_str());
  output->Close();
 
}

void MakeAllFinalState(){
  SignalHisto_data("4m");
  SignalHisto_data("4e");
  SignalHisto_data("2e2m"); 
 
}
