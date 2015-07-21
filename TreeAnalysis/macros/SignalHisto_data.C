#include <TFile.h>
#include <TROOT.h>
#include <TStyle.h>
#include <TH1.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <iostream>
#include <string>


void SignalHisto_data(string var = "Mass",string finalstate = "4m", bool syst = 0)
{
  gROOT->Reset();  
  gROOT->SetStyle("Plain");   
  gStyle->SetOptStat(0); 

  cout << "=================================================================================================================" << endl; 
  cout << "                                       " << finalstate.c_str() << endl; 
  cout << "=================================================================================================================" << endl; 
  TFile *output;
  string fileName = var+"_test/DataToUnfold.root";
  string fileNameSyst =  var+"_test/DataToUnfold_syst.root";

  if(syst == 0) {output = new TFile(fileName.c_str(), "UPDATE");}
  else {output = new TFile(fileNameSyst.c_str(), "UPDATE");}

  //Data (ZZRecoAnalyzer not yet done FIXME)
  TFile *data = new TFile("../results/ZZRecoAnalyzer_SR/data.root");
  //reducible background
  TFile *red = new TFile("../results/ZZRecoAnalyzer_CR/data.root"); 
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
  TH1 * h_data_irrp;
  TH1 * h_data_irrm;
  TH1 * h_data_redp;
  TH1 * h_data_redm;
  string variable;
  string XaxisTitle;
  string fs;  
  int b =0;

  if(finalstate == "4m") fs = "4#mu";
  else if(finalstate == "2e2m") fs = "2e2#mu";
  else fs = finalstate;
  
  if(var == "Mass") {
    variable = var;
    XaxisTitle = "m_{" + fs + "}";
    b = 9;
   }
  else if(var == "Jets"){
    variable = "Jets_JERCentralSmear"; 
    XaxisTitle = "N Jets (" + finalstate + " final state)";
    b = 5;
  }

  string histoName = "ZZTo" + finalstate + "_" + var + "_01";
  string histoMCName = "ZZTo" + finalstate + "_" + variable + "_01";

  h_totdata = (TH1*) data->Get(histoName.c_str()); 
  h_data = (TH1*)h_totdata->Clone("h_totdata");
  h_red = (TH1*) red->Get(histoName.c_str()); 
  h_ttZ = (TH1*) ttZ->Get(histoMCName.c_str()); 
  h_ttWW = (TH1*) ttWW->Get(histoMCName.c_str()); 
  h_WWZ = (TH1*) WWZ->Get(histoMCName.c_str()); 
  
  if(syst == 1){
    h_data_irrp = (TH1*)h_totdata->Clone("h_totdata");
    h_data_irrm = (TH1*)h_totdata->Clone("h_totdata");
    h_data_redp = (TH1*)h_totdata->Clone("h_totdata");
    h_data_redm = (TH1*)h_totdata->Clone("h_totdata");
  }
  
  // float err_tot = 0;
  // float err_sist = 0; 
  // float err_stat = 0;
  float err_red = 0;
  float err_irr = 0;
  float dataminusbkg = 0;
  float dmb_irrp =0;
  float dmb_irrm =0;
  float dmb_redp =0;
  float dmb_redm =0;

  for(int i =1; i<b; i++){
    dataminusbkg = 0;
    dataminusbkg = h_totdata->GetBinContent(i)- h_red->GetBinContent(i)- h_ttZ->GetBinContent(i)-h_ttWW->GetBinContent(i)-h_WWZ->GetBinContent(i);
    if(dataminusbkg>0.) h_data-> SetBinContent(i,dataminusbkg);
    else h_data-> SetBinContent(i,0.);
    
    err_irr = 0;
    err_red = 0; 
    err_irr = sqrt(h_ttZ->GetBinError(i)*h_ttZ->GetBinError(i)+h_ttWW->GetBinError(i)*h_ttWW->GetBinError(i)+h_WWZ->GetBinError(i)*h_WWZ->GetBinError(i));
    err_red = h_red->GetBinError(i);
    
    if(syst == 1){ 
      dmb_irrp = dataminusbkg + err_irr;
      dmb_irrm = dataminusbkg - err_irr;
      dmb_redp = dataminusbkg + err_red;
      dmb_redm = dataminusbkg - err_red;
      
      if(dmb_irrp>0.) h_data_irrp-> SetBinContent(i,dmb_irrp);
      else h_data_irrp-> SetBinContent(i,0.);
      if(dmb_irrm>0.) h_data_irrm-> SetBinContent(i,dmb_irrm);
      else h_data_irrm-> SetBinContent(i,0.);
      if(dmb_redp>0.) h_data_redp-> SetBinContent(i,dmb_redp);
      else h_data_redp-> SetBinContent(i,0.);
      if(dmb_redm>0.) h_data_redm-> SetBinContent(i,dmb_redm);
      else h_data_redm-> SetBinContent(i,0.);
    }

    cout << "=================================================================================================================" << endl; 
    cout <<  "                                   bin " << i << endl; 
    cout << "=================================================================================================================" << endl; 
    //some information
    std::cout << "bin " << i << " data = "<< h_totdata->GetBinContent(i) << " +- " << h_totdata->GetBinError(i) << " red = " << h_red->GetBinContent(i) << " +- " << h_red->GetBinError(i) <<   " ttZ = " << h_ttZ->GetBinContent(i) << " +- " << h_ttZ->GetBinError(i) <<  " ttWW = " << h_ttWW->GetBinContent(i) << " +- " << h_ttWW->GetBinError(i) <<   " WWZ = " << h_WWZ->GetBinContent(i)   << " +- " << h_WWZ->GetBinError(i) << " data-bkg = " << dataminusbkg <<" +- " << h_data->GetBinError(i) << " " << h_data->GetBinContent(i)<< std::endl;
  }
  
  TLegend *leg1 = new TLegend(0.65,0.65,0.6,0.85);
  leg1->SetFillColor(kWhite);
  leg1->SetBorderSize(0);
  leg1->SetTextSize(0.03); 
 
  TCanvas *c1 = new TCanvas ("c1","c1");
  h_totdata->Draw(); 
  h_data->SetLineColor(1);
  h_data->Draw("SAME");
  h_data->SetLineColor(2);

  leg1->AddEntry(h_totdata,"full dataset","l"); 
  leg1->Draw(); 
  leg1->AddEntry(h_data,"data-background","l"); 
  leg1->Draw("SAME");
   
  double err_data =0;
  double err_finaldata =0;
  double err_red_tot = 0;
  double err_ttZ = 0;
  double err_ttWW = 0;
  double err_WWZ = 0;
  double tot_data =  h_totdata->IntegralAndError(0,9,err_data);
  double tot_red =  h_red->IntegralAndError(0,9,err_red_tot);
  double tot_ttZ =  h_ttZ->IntegralAndError(0,9,err_ttZ);
  double tot_ttWW =  h_ttWW->IntegralAndError(0,9,err_ttWW);
  double tot_WWZ =  h_WWZ->IntegralAndError(0,9,err_WWZ);
  double tot_irr = tot_ttZ + tot_ttWW + tot_WWZ ;
  double err_irr_tot = sqrt(err_ttZ*err_ttZ + err_ttWW*err_ttWW + err_WWZ*err_WWZ);
  double finaldata = h_data->IntegralAndError(0,9,err_finaldata);
 
  std::cout << finalstate.c_str() << std::endl;
  std::cout <<  "    tot Data Yield = " << tot_data << " +- " << err_data <<std::endl;
  std::cout <<  "    red_tot = " << tot_red << " +- " << err_red  << std::endl; 
  std::cout <<  "    ttZ_tot = " << tot_ttZ << " +- " << err_ttZ << std::endl; 
  std::cout <<  "    ttWW_tot = " << tot_ttWW <<" +- " << err_ttWW << std::endl;
  std::cout <<  "    WWZ_tot = " << tot_WWZ << " +- " << err_WWZ  << std::endl; 
  std::cout <<  "    tot_irr = " << tot_irr << " +- " << err_irr_tot << std::endl; 
  std::cout <<  "    final data = " << finaldata << " +- " << err_finaldata << std::endl;
  std::cout <<   "==============================================" << std::endl;

  
  if(syst ==1){
    h_data_redp->GetXaxis()->SetTitle(XaxisTitle.c_str());
    
    TLegend *leg = new TLegend(0.65,0.65,0.6,0.85);
    leg->SetFillColor(kWhite);
    leg->SetBorderSize(0);
    leg->SetTextSize(0.03); 
  
    TCanvas *c = new TCanvas ("c","c");
  
    h_data_redp->Draw();
    h_data_redp->SetLineColor(4);
    h_data_redp->SetMarkerStyle(7);
    h_data_redp->SetMarkerColor(4); 
    //h_data_redp->SetMarkerSize(10);
    h_data_redm->Draw("SAME");
    h_data_redm->SetLineColor(kMagenta); 
    h_data_redm->SetMarkerStyle(7); 
    h_data_redm->SetMarkerColor(kMagenta); 
    //h_data_redm->SetMarkerSize(10);
    h_data_irrp->Draw("SAME");
    h_data_irrp->SetLineColor(2); 
    h_data_irrp->SetMarkerStyle(7); 
    h_data_irrp->SetMarkerColor(2);
    //h_data_irrp->SetMarkerSize(10);
    h_data_irrm->Draw("SAME");
    h_data_irrm->SetLineColor(3); 
    h_data_irrm->SetMarkerStyle(7); 
    h_data_irrm->SetMarkerColor(3);
    //h_data_irrm->SetMarkerSize(10);
    h_data->Draw("SAME");
    h_data->SetLineColor(1); 
    h_data->SetMarkerStyle(7);
    h_data->SetMarkerColor(1);
    
    leg->AddEntry(h_data,"data-background","l"); 
    leg->Draw();
    leg->AddEntry(h_data_irrp,"data-background (+ err_irr)","l"); 
    leg->Draw("SAME");
    leg->AddEntry(h_data_irrm,"data-background (- err_irr)","l"); 
    leg->Draw("SAME");
    leg->AddEntry(h_data_redp,"data-background (+ err_red)","l"); 
    leg->Draw("SAME");
    leg->AddEntry(h_data_redm,"data-background (- err_red)","l"); 
    leg->Draw("SAME");
  
  }  

  string dataName = "DataminusBkg_"+var+"_ZZTo"+finalstate;
  string TotdataName = "TotData_"+var+"_ZZTo"+finalstate;
  string dataIrrpName = "DataminusBkg_irrp_"+var+"_ZZTo"+finalstate; 
  string dataIrrmName = "DataminusBkg_irrm_"+var+"_ZZTo"+finalstate; 
  string dataRedpName = "DataminusBkg_redp_"+var+"_ZZTo"+finalstate; 
  string dataRedmName = "DataminusBkg_redm_"+var+"_ZZTo"+finalstate;
  
 
    output->cd();    
    if(syst ==0){
      h_data->Write(dataName.c_str());
      h_totdata->Write(TotdataName.c_str());
    }
    else{
      h_data_irrp->Write(dataIrrpName.c_str());
      h_data_irrm->Write(dataIrrmName.c_str());
      h_data_redp->Write(dataRedpName.c_str());
      h_data_redm->Write(dataRedmName.c_str());
    }
    output->Close();

 
}

void MakeAllFinalStates(string var="Mass",bool syst =0){
  SignalHisto_data(var.c_str(),"4m",syst);
  SignalHisto_data(var.c_str(),"4e",syst);
  SignalHisto_data(var.c_str(),"2e2m",syst); 
 
}
