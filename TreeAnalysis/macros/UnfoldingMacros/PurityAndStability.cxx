#if !defined(__CINT__) || defined(__MAKECINT__)
#include "PurityAndStability.h"
#include "tdrstyle.C"
#include "CMS_lumi.C"
#endif

using namespace std;

ClassImp(PurityAndStability)

PurityAndStability::PurityAndStability(): TObject()
{}

PurityAndStability::PurityAndStability(bool mad): TObject() 
{ 
  if(mad ==1) mc = "_Mad";
  else mc = "_Pow";
}

PurityAndStability::~PurityAndStability(){}

//Build purity and stability distributions
void PurityAndStability::Build(string var, string finalstate)
{
  fileName = var+"_test/MatrixPurityStability"+mc+".root";
  matrixFileName =  var+"_test/matrices"+mc+".root";
  output = new TFile(fileName.c_str(), "UPDATE");
  matrixFile = new TFile(matrixFileName.c_str());
 
  matrixName = "ResMat_qqggJJ_"+var+"_ZZTo"+finalstate+"_st_01";
  histoName = var+"_qqggJJ_ZZTo"+finalstate+"_st_01";
  h_Resmat = (TH2*) matrixFile->Get(matrixName.c_str()); 
  h_purity = (TH1*) matrixFile->Get(histoName.c_str()); 
  h_stability = (TH1*) matrixFile->Get(histoName.c_str());  
  
  int b = 0;

  if(var == "Mass")  b = 9;
  else if(var == "nJets" || var == "nJets_Central") b = 5;
  else if(var == "Mjj" || var == "CentralMjj" || var == "Deta" || var == "CentralDeta") b = 3;
  else if(var == "PtJet1") b = 6;
  else if(var == "PtJet2") b = 5;
  else if(var == "EtaJet1" || var == "EtaJet2" ) b = 6;
  else if(var == "dRZZ") b = 6;

  float binContent = 0; 
  float binError = 0;
  float all_reco = 0;
  float all_gen = 0;
  float p_err =0;  
  float s_err =0;  
  double reco_err =0;  
  double gen_err =0;  
  
  for(int i =1; i<b; i++){
    reco_err=0;
    gen_err=0;
    binContent=0;  
    binError=0;    
    all_reco=0;
    all_gen=0; 
    p_err=0;
    s_err =0;
    
    binError = h_Resmat->GetBinError(i,i);
    binContent = h_Resmat->GetBinContent(i,i);
    all_reco = h_Resmat->IntegralAndError(i,i,1,b,reco_err);
    all_gen = h_Resmat->IntegralAndError(1,b,i,i,gen_err);
    p_err = sqrt((binError/binContent)*(binError/binContent) + (reco_err/all_reco)*(reco_err/all_reco))*(binContent/all_reco);    
    s_err = sqrt((binError/binContent)*(binError/binContent) + (gen_err/all_gen)*(gen_err/all_gen))*(binContent/all_gen);    


    // std::cout <<"bin " << i <<  " all_reco  = " << all_reco << " +- " << reco_err << " p= "<<  binContent/all_reco << " +- " << p_err << std::endl;
    // std::cout <<"bin " << i <<" all_gen  = " << all_gen <<" +- " << gen_err << " s= " << binContent/all_gen <<" +- " << s_err  << std::endl;
    
    if(all_reco!=0){
      h_purity->SetBinContent(i,binContent/all_reco);
      h_purity->SetBinError(i,p_err);
    }
    else h_purity->SetBinContent(i,0.);
    if(all_gen!=0) {
      h_stability->SetBinContent(i,binContent/all_gen);
      h_stability->SetBinError(i,s_err);
    }
    else h_stability->SetBinContent(i,0.);
  }

  h_purity->SetMaximum(1.1);
  h_purity->SetMinimum(0.); 
  h_stability->SetMaximum(1.1);
  h_stability->SetMinimum(0.);

  string purityName = "hPurity_ZZTo" + finalstate + "_01";
  string stabilityName = "hStability_ZZTo" + finalstate + "_01";

  output->cd(); 
  h_purity->Write(purityName.c_str()); 
  h_stability->Write(stabilityName.c_str()); 
  output->Close();
  
}

//Plot Purity and Stability distributions
void PurityAndStability::Plot(string var,string finalstate, string path) 
{
  gROOT->Reset();  
  gROOT->SetStyle("Plain");   
  gStyle->SetOptStat(0);
  
  string xAxis; 
  string fs;

  if(finalstate == "4m") fs = "4#mu";
  else if(finalstate == "2e2m") fs = "2e2#mu";
  else fs = finalstate;
   
  fileNameMad = var+"_test/MatrixPurityStability_Mad.root";
  fileNamePow =  var+"_test/MatrixPurityStability_Pow.root";
  
  madgraph = new TFile(fileNameMad.c_str()); 
  powheg = new TFile(fileNamePow.c_str()); 
  
  pName = "hPurity_ZZTo"+finalstate+"_01";
  sName = "hStability_ZZTo"+finalstate+"_01";
 
  p_mad = (TH1*) madgraph->Get(pName.c_str());
  s_mad = (TH1*) madgraph->Get(sName.c_str());
  p_pow = (TH1*) powheg->Get(pName.c_str());
  s_pow = (TH1*) powheg->Get(sName.c_str());

  if(var =="Mass"){
    xAxis = "reco m_{"+fs+"}";
   }
  else if(var =="nJets"){
    xAxis = "reco Njets";
  }
  else if(var =="Mjj"){
    xAxis = "reco m_{jj}";
  }
  else if(var =="Deta"){
    xAxis = "reco #Delta#eta_{jj}";
  }
  else if(var =="nJets_Central"){
    xAxis = "reco Ncentraljets";
  }
  else if(var =="Mjj_Central"){
    xAxis = "reco m_{jj}";
  }
  else if(var =="Deta_Central"){
    xAxis = "reco #Delta#eta_{jj}";
  }
  else if(var =="PtJet1"){
    xAxis = "reco p_{T}^{jet1}";
  }
  else if(var =="PtJet2"){
    xAxis = "reco p_{T}^{jet2}";
  }
 else if(var =="EtaJet1"){
   xAxis = "reco |#eta^{jet1}|";
 }
 else if(var =="EtaJet2"){
   xAxis = "reco |#eta^{jet2}|";
 } 
 else if(var =="dRZZ"){
   xAxis = "reco #DeltaR(Z_1,Z_2)";
 }
  
  TCanvas *c_p = new TCanvas ("c_p","c_p");
  TLegend *leg_p = new TLegend(0.65,0.35,0.6,0.43); 
  leg_p->SetFillColor(kWhite);
  leg_p->SetBorderSize(0);
  leg_p->SetTextSize(0.03);
 
  string sTitle = "";  
  string pTitle = "";
 
  p_mad->Draw("E");  
  p_mad->SetLineColor(1);
  p_mad->SetTitle(pTitle.c_str());
  p_mad->GetXaxis()->SetTitle(xAxis.c_str());
 
  p_pow->Draw("E SAME");
  p_pow->SetLineColor(2);
  
  leg_p->AddEntry(p_mad,"Madgraph set","l"); 
  leg_p->Draw(); 
  leg_p->AddEntry(p_pow,"Powheg set","l"); 
  leg_p->Draw("SAME");
  
  TCanvas *c_s = new TCanvas ("c_s","c_s");
  TLegend *leg_s = new TLegend(0.65,0.35,0.6,0.43); 
  leg_s->SetFillColor(kWhite);
  leg_s->SetBorderSize(0);
  leg_s->SetTextSize(0.03);

  s_mad->Draw("E");  
  s_mad->SetLineColor(1);
  s_mad->SetTitle(sTitle.c_str());
  s_mad->GetXaxis()->SetTitle(xAxis.c_str());
 
  s_pow->Draw("E SAME");
  s_pow->SetLineColor(2);
  
  leg_s->AddEntry(s_mad,"Madgraph set","l"); 
  leg_s->Draw(); 
  leg_s->AddEntry(s_pow,"Powheg set","l"); 
  leg_s->Draw("SAME");
  
 
  string p_png = "~/www/VBS/"+path+"/"+var+"/PurityAndStability/"+"purity_MadPow_"+finalstate+"_"+var+".png";
  string s_png = "~/www/VBS/"+path+"/"+var+"/PurityAndStability/"+"stability_MadPow_"+finalstate+"_"+var+".png";
  string p_pdf = "~/www/VBS/"+path+"/"+var+"/PurityAndStability/"+"purity_MadPow_"+finalstate+"_"+var+".pdf";
  string s_pdf = "~/www/VBS/"+path+"/"+var+"/PurityAndStability/"+"stability_MadPow_"+finalstate+"_"+var+".pdf";
  c_p->Print(p_png.c_str());
  c_s->Print(s_png.c_str());
  c_p->Print(p_pdf.c_str());
  c_s->Print(s_pdf.c_str());
}


void PurityAndStability::Plot_PAS(string var,string finalstate, string path) 
{ 
  gROOT->Reset();  
  gROOT->SetStyle("Plain");   
  gStyle->SetOptStat(0); 
  
  setTDRStyle();
  int iPeriod = 2; 
  int iPos = 11; 
  writeExtraText = true;    
  extraText  = "Simulation";
  extraText2 = "";

  string xAxis;;
  string fs;

  fileName = var+"_test/MatrixPurityStability"+mc+".root";
  file = new TFile(fileName.c_str(), "UPDATE");
 
  pName = "hPurity_ZZTo"+finalstate+"_01";
  sName = "hStability_ZZTo"+finalstate+"_01";
  
  p_mad = (TH1*) file->Get(pName.c_str());
  s_mad = (TH1*) file->Get(sName.c_str());
 
  if(finalstate == "4m") fs = "4#mu";
  else if(finalstate == "2e2m") fs = "2e2#mu";
  else fs = finalstate;
  
  if(var == "nJets" || var == "nJets_Central"){
    p_mad->GetXaxis()->SetBinLabel(1,"0");
    p_mad->GetXaxis()->SetBinLabel(2,"1");
    p_mad->GetXaxis()->SetBinLabel(3,"2");
    p_mad->GetXaxis()->SetBinLabel(4,">2");  
    //p_mad->GetXaxis()->SetLabelFont(42);
    //p_mad->GetXaxis()->SetLabelOffset(0.02);
    p_mad->GetXaxis()->SetLabelSize(0.05);
  }

  if(var =="Mass"){
    xAxis = "reco m_{"+fs+"} [GeV]";
   }
  else if(var =="nJets"){
    xAxis = "reco N jets (|#eta^{jet}|<4.7)";
  }
  else if(var =="Mjj"){
    xAxis = "reco m_{jj} (|#eta^{jet}|<4.7) [GeV]";
  }
  else if(var =="Deta"){
    xAxis = "reco #Delta#eta_{jj} (|#eta^{jet}|<4.7)";
  }
  else if(var =="nJets_Central"){
    xAxis = "reco N jets (|#eta^{jet}|<2.4)";
  }
  else if(var =="Mjj_Central"){
    xAxis = "reco m_{jj} (|#eta^{jet}|<2.4) [GeV]";
  }
  else if(var =="Deta_Central"){
    xAxis = "reco #Delta#eta_{jj} (|#eta^{jet}|<2.4)";
  }
  else if(var =="PtJet1"){
    xAxis = "reco p_{T}^{jet1} [GeV]";
  }
  else if(var =="PtJet2"){
    xAxis = "reco p_{T}^{jet2} [GeV]";
  }
 else if(var =="EtaJet1"){
   xAxis = "reco |#eta^{jet1}|";
 }
 else if(var =="EtaJet2"){
   xAxis = "reco |#eta^{jet2}|";
 }
 else if(var =="dRZZ"){
   xAxis = "reco #DeltaR(Z_1,Z_2)";
 }
  TCanvas *c_p = new TCanvas ("c_p","c_p");
  TLegend *leg_p = new TLegend(0.70,0.35,0.6,0.45); 
  leg_p->SetFillColor(kWhite);
  leg_p->SetBorderSize(0);
  leg_p->SetTextSize(0.04);
  leg_p->SetTextFont(42);   

  p_mad->Draw("E");  
  p_mad->SetMaximum(1.2);  
  p_mad->SetLineColor(1);
  p_mad->SetMarkerStyle(8);
  //p_mad->SetTitle(pTitle.c_str());
  p_mad->SetTitle("");
  //hReco->GetXaxis()->SetRange(0,25);
  p_mad->GetXaxis()->SetTitle(xAxis.c_str());
  p_mad->GetXaxis()->SetTitleOffset(1.2);

  s_mad->Draw("E SAME");
  s_mad->SetLineColor(2);
  s_mad->SetMarkerStyle(8);
  s_mad->SetMarkerColor(2);

  leg_p->AddEntry(p_mad,"Purity","p"); 
  leg_p->Draw(); 
  leg_p->AddEntry(s_mad,"Stability","p"); 
  leg_p->Draw("SAME");
  
  //lumiTextSize     = 0.4;
  // cmsTextSize      = 0.48;
  // extraOverCmsTextSize  = 0.80;//0.63;
  CMS_lumi(c_p,iPeriod,iPos);
  // pad1->Update();
  // pad1->RedrawAxis();
  //pad1->GetFrame()->Draw();

  string p_pdf ="~/www/VBS/"+path+"/"+var+"/PurityAndStability/"+ "PurityStability_"+finalstate+"_"+var+mc+".pdf";
  string p_png = "~/www/VBS/"+path+"/"+var+"/PurityAndStability/"+"PurityStability_"+finalstate+"_"+var+mc+".png";

  c_p->Print(p_pdf.c_str());
  c_p->Print(p_png.c_str());
}
