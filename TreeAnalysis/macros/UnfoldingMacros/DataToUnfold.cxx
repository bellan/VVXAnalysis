#if !defined(__CINT__) || defined(__MAKECINT__)

#include "DataToUnfold.h"
#include "tdrstyle.h"
#include "CMS_lumi.h"
// #include "CMS_lumi.C"

#endif

using namespace std;

ClassImp(DataToUnfold)


DataToUnfold::DataToUnfold(): TObject() 
{  
  data = new TFile("../../results/ZZRecoAnalyzer_SR/data.root");
  //reducible background
  red3P1F     = new TFile("../../results/ZZRecoAnalyzer_CR3P1F/data.root"); 
  red3P1FqqZZ = new TFile("../../results/ZZRecoAnalyzer_CR3P1F/ZZTo4lamcatnlo.root");
  red3P1FggZZ = new TFile("../../results/ZZRecoAnalyzer_CR3P1F/gg_4l.root"); 
  red2P2F     = new TFile("../../results/ZZRecoAnalyzer_CR2P2F/data.root"); 
  //irreducible backgrounds
  ttZ = new TFile("../../results/ZZRecoAnalyzer_SR/TTZToLL.root");
  WWZ = new TFile("../../results/ZZRecoAnalyzer_SR/WWZ.root");
  //   ZZbkg = new TFile("../../results/ZZRecoAnalyzer_SR/sig_pow.root"); //new

}


DataToUnfold::~DataToUnfold(){}


//Build data and data-minus-bkg distributions for the 4e, 4mu and 2e2mu final states. 
//Build data distributions to estimate the reducible and irreducible background systematic uncertainties 
void DataToUnfold::Build(string var, string finalstate)
{
  output = new TFile((var+"_test/DataToUnfold.root").c_str(), "UPDATE");
  output_syst =   new TFile((var+"_test/DataToUnfold_syst.root").c_str(), "UPDATE");
  
  variable = var;

  histoName          = "ZZTo" + finalstate + "_" + var + "_01";
  histoMCName        = "ZZTo" + finalstate + "_" + variable + "_01";
  dataName           = "DataminusBkg_"+var+"_ZZTo"+finalstate;
  TotdataName        = "TotData_"+var+"_ZZTo"+finalstate;
  dataIrrpName       = "DataminusBkg_irrp_"+var+"_ZZTo"+finalstate; 
  dataIrrmName       = "DataminusBkg_irrm_"+var+"_ZZTo"+finalstate; 
  dataRedpName       = "DataminusBkg_redp_"+var+"_ZZTo"+finalstate; 
  dataRedmName       = "DataminusBkg_redm_"+var+"_ZZTo"+finalstate;

  h_totdata  = (TH1*) data->Get(histoName.c_str()); 
  h_data     = (TH1*) h_totdata->Clone("h_totdata");

  h_red          = (TH1*) red3P1F->Get(histoName.c_str()); 
  h_red3P1FqqZZ  = (TH1*) red3P1FqqZZ->Get(histoName.c_str()); 
  h_red3P1FggZZ  = (TH1*) red3P1FggZZ->Get(histoName.c_str()); 
  h_red2P2F      = (TH1*) red2P2F->Get(histoName.c_str()); 
  
  
  if(h_red          ==NULL) {cout <<"histo "<<histoName.c_str()<<" in "<<  red3P1F->GetName()    <<" is null. Abort"<<endl; abort();} 
  if(h_red3P1FqqZZ  ==NULL) {cout <<"histo "<<histoName.c_str()<<" in "<<  red3P1FqqZZ->GetName()<<" is null. Abort"<<endl; abort();} 
  if(h_red3P1FggZZ  ==NULL) {cout <<"histo "<<histoName.c_str()<<" in "<<  red3P1FggZZ->GetName()<<" is null. Abort"<<endl; abort();} 
  if(h_red2P2F      ==NULL) {cout <<"histo "<<histoName.c_str()<<" in "<<  red2P2F->GetName()    <<" is null. Abort"<<endl; abort();} 
  

  h_red->Add(h_red3P1FqqZZ,-1);
  h_red->Add(h_red3P1FggZZ,-1);
  h_red->Add(h_red2P2F);

  h_red_up          = (TH1*) red3P1F->Get(("ZZTo" + finalstate + "_" + var + "_RedUp_01").c_str()); 
  h_red3P1FqqZZ_up  = (TH1*) red3P1FqqZZ->Get(("ZZTo" + finalstate + "_" + var + "_RedUp_01").c_str()); 
  h_red3P1FggZZ_up  = (TH1*) red3P1FggZZ->Get(("ZZTo" + finalstate + "_" + var + "_RedUp_01").c_str()); 
  h_red2P2F_up      = (TH1*) red2P2F->Get(("ZZTo" + finalstate + "_" + var + "_RedUp_01").c_str()); 


  if(h_red_up          ==NULL) {cout <<"histo "<<histoName.c_str()<<" in "<<  red3P1F->GetName()    <<" is null. Abort"<<endl; abort();} 
  if(h_red3P1FqqZZ_up  ==NULL) {cout <<"histo "<<histoName.c_str()<<" in "<<  red3P1FqqZZ->GetName()<<" is null. Abort"<<endl; abort();} 
  if(h_red3P1FggZZ_up  ==NULL) {cout <<"histo "<<histoName.c_str()<<" in "<<  red3P1FggZZ->GetName()<<" is null. Abort"<<endl; abort();} 
  if(h_red2P2F_up      ==NULL) {cout <<"histo "<<histoName.c_str()<<" in "<<  red2P2F->GetName()    <<" is null. Abort"<<endl; abort();} 
  

  h_red_up->Add(h_red3P1FqqZZ_up,-1);
  h_red_up->Add(h_red3P1FggZZ_up,-1);
  h_red_up->Add(h_red2P2F_up);

  h_red_down          = (TH1*) red3P1F->Get(("ZZTo" + finalstate + "_" + var + "_RedDn_01").c_str()); 
  h_red3P1FqqZZ_down  = (TH1*) red3P1FqqZZ->Get(("ZZTo" + finalstate + "_" + var + "_RedDn_01").c_str()); 
  h_red3P1FggZZ_down  = (TH1*) red3P1FggZZ->Get(("ZZTo" + finalstate + "_" + var + "_RedDn_01").c_str()); 
  h_red2P2F_down      = (TH1*) red2P2F->Get(("ZZTo" + finalstate + "_" + var + "_RedDn_01").c_str()); 

  if(h_red_down          ==NULL) {cout <<"histo "<<histoName.c_str()<<" in "<<  red3P1F->GetName()    <<" is null. Abort"<<endl; abort();} 
  if(h_red3P1FqqZZ_down  ==NULL) {cout <<"histo "<<histoName.c_str()<<" in "<<  red3P1FqqZZ->GetName()<<" is null. Abort"<<endl; abort();} 
  if(h_red3P1FggZZ_down  ==NULL) {cout <<"histo "<<histoName.c_str()<<" in "<<  red3P1FggZZ->GetName()<<" is null. Abort"<<endl; abort();} 
  if(h_red2P2F_down      ==NULL) {cout <<"histo "<<histoName.c_str()<<" in "<<  red2P2F->GetName()    <<" is null. Abort"<<endl; abort();} 


  h_red_down->Add(h_red3P1FqqZZ_down,-1);
  h_red_down->Add(h_red3P1FggZZ_down,-1);
  h_red_down->Add(h_red2P2F_down);

  h_ttZ      = (TH1*) ttZ->Get(histoMCName.c_str()); 
  h_WWZ      = (TH1*) WWZ->Get(histoMCName.c_str()); 
  //  h_ZZBkg    = (TH1*) ZZbkg->Get(("ZZTo" + finalstate + "_" + var + "_01_nofr").c_str());  //new

  if(h_totdata  ==NULL) {cout <<"histo "<<histoName.c_str()<<" in "<<h_totdata->GetName()<<" is null. Abort"<<endl; abort();} 
  if(h_red      ==NULL) {cout <<"histo "<<histoName.c_str()<<" in "<<h_red->GetName()<<" is null. Abort"<<endl; abort();}  
  if(h_ttZ      ==NULL) {cout <<"histo "<<histoName.c_str()<<" in "<<h_ttZ->GetName()<<" is null. Abort"<<endl; abort();} 
  if(h_WWZ      ==NULL) {cout <<"histo "<<histoName.c_str()<<" in "<<h_WWZ->GetName()<<" is null. Abort"<<endl; abort();} 
  
  h_data_irrp = (TH1*)h_totdata->Clone("h_totdata");
  h_data_irrm = (TH1*)h_totdata->Clone("h_totdata");
  h_data_redp = (TH1*)h_totdata->Clone("h_totdata");
  h_data_redm = (TH1*)h_totdata->Clone("h_totdata");
   
  float err_red = 0;
  float err_irr = 0;
  float dataminusbkg = 0;
  float dataminusbkg_redUp = 0;
  float dataminusbkg_redDn = 0;
  float dmb_irrp =0;
  float dmb_irrm =0;
  float dmb_redp =0;
  float dmb_redm =0;

  int b = h_totdata->GetNbinsX();
  for(int i =1; i<=b; i++){
    dataminusbkg = 0;
    err_irr = 0;
    err_red = 0;
    err_red_tot = 0;
    dmb_irrp =0;
    dmb_irrm =0;
    dmb_redp =0;
    dmb_redm =0;
    
    //    cout<<"all "<< h_totdata->GetBinContent(i) <<" red "<<h_red->GetBinContent(i)<<" irr "<<h_ttZ->GetBinContent(i)<<" irr 2 "<<h_WWZ->GetBinContent(i)<<endl; 
    //    dataminusbkg = h_totdata->GetBinContent(i)- h_red->GetBinContent(i)- h_ttZ->GetBinContent(i)-h_ttWW->GetBinContent(i)-h_WWZ->GetBinContent(i); //FIXME

    // dataminusbkg       = h_totdata->GetBinContent(i)- h_red->GetBinContent(i)- h_ttZ->GetBinContent(i)-h_WWZ->GetBinContent(i) -  h_ZZBkg->GetBinContent(i);  //new
    // dataminusbkg_redUp = h_totdata->GetBinContent(i)- h_red_up->GetBinContent(i)- h_ttZ->GetBinContent(i)-h_WWZ->GetBinContent(i) - h_ZZBkg->GetBinContent(i); 
    // dataminusbkg_redDn = h_totdata->GetBinContent(i)- h_red_down->GetBinContent(i)- h_ttZ->GetBinContent(i)-h_WWZ->GetBinContent(i) - h_ZZBkg->GetBinContent(i); 

    //    cout<< h_ZZBkg->GetBinContent(i)<<endl;
    dataminusbkg = h_totdata->GetBinContent(i)- h_red->GetBinContent(i)- h_ttZ->GetBinContent(i)-h_WWZ->GetBinContent(i); 
    dataminusbkg_redUp = h_totdata->GetBinContent(i)- h_red_up->GetBinContent(i)- h_ttZ->GetBinContent(i)-h_WWZ->GetBinContent(i); 
    dataminusbkg_redDn = h_totdata->GetBinContent(i)- h_red_down->GetBinContent(i)- h_ttZ->GetBinContent(i)-h_WWZ->GetBinContent(i); 

    //    dataminusbkg = h_totdata->GetBinContent(i)- h_red->GetBinContent(i)- h_Irr->GetBinContent(i);
    if(dataminusbkg>0.) h_data-> SetBinContent(i,dataminusbkg);
    else h_data-> SetBinContent(i,0.);


    err_irr = sqrt(h_ttZ->GetBinError(i)*h_ttZ->GetBinError(i)+h_WWZ->GetBinError(i)*h_WWZ->GetBinError(i));

    //    err_irr = sqrt(h_Irr->GetBinError(i)*h_Irr->GetBinError(i));
    //    err_red_tot=+ err_red*err_red;

    err_red  = h_red->GetBinContent(i)*0.3; //To be changed
    dmb_irrp = dataminusbkg + err_irr;
    dmb_irrm = dataminusbkg - err_irr;

    //    dmb_redp = dataminusbkg + err_red;
    // dmb_redm = dataminusbkg - err_red;
    
    if(dmb_irrp>0.) h_data_irrp-> SetBinContent(i,dmb_irrp);
    else h_data_irrp-> SetBinContent(i,0.);
    if(dmb_irrm>0.) h_data_irrm-> SetBinContent(i,dmb_irrm);
    else h_data_irrm-> SetBinContent(i,0.);
    if(dataminusbkg_redUp>0.) h_data_redp-> SetBinContent(i, dataminusbkg_redUp);
    else h_data_redp-> SetBinContent(i,0.);
    if(dataminusbkg_redDn >0.) h_data_redm-> SetBinContent(i,dataminusbkg_redDn);
    else h_data_redm-> SetBinContent(i,0.);
    
    // cout << "=================================================================================================================" << endl; 
    // cout <<  "                                   bin " << i << endl; 
    // cout << "=================================================================================================================" << endl; 
    // //some information
    // std::cout << "bin " << i << " data = "<< h_totdata->GetBinContent(i) << " +- " << h_totdata->GetBinError(i) << " red = " << h_red->GetBinContent(i) << " +- " << h_red->GetBinError(i) <<   " ttZ = " << h_ttZ->GetBinContent(i) << " +- " << h_ttZ->GetBinError(i) <<  " ttWW = " << h_ttWW->GetBinContent(i) << " +- " << h_ttWW->GetBinError(i) <<   " WWZ = " << h_WWZ->GetBinContent(i)   << " +- " << h_WWZ->GetBinError(i) << " data-bkg = " << dataminusbkg <<" +- " << h_data->GetBinError(i) << " " << h_data->GetBinContent(i)<< std::endl;
  }
  
  // //Print some values:
  // double err_data =0;
  // double err_finaldata =0;
  // double sqrt_err_red_tot = sqrt(err_red_tot);
  // double err_ttZ = 0;
  // double err_ttWW = 0;
  // double err_WWZ = 0;
  // double tot_data =  h_totdata->IntegralAndError(0,9,err_data);
  // double tot_red =  h_red->IntegralAndError(0,9,sqrt_err_red_tot);
  // double tot_ttZ =  h_ttZ->IntegralAndError(0,9,err_ttZ);
  // double tot_ttWW =  h_ttWW->IntegralAndError(0,9,err_ttWW);
  // double tot_WWZ =  h_WWZ->IntegralAndError(0,9,err_WWZ);
  // double tot_irr = tot_ttZ + tot_ttWW + tot_WWZ ;
  // double err_irr_tot = sqrt(err_ttZ*err_ttZ + err_ttWW*err_ttWW + err_WWZ*err_WWZ);
  // double finaldata = h_data->IntegralAndError(0,9,err_finaldata);
  
  // std::cout << finalstate.c_str() << std::endl;
  // std::cout <<  "    tot Data Yield = " << tot_data << " +- " << err_data <<std::endl;
  // std::cout <<  "    red_tot = " << tot_red << " +- " << sqrt_err_red_tot  << "(" <<sqrt_err_red_tot/tot_data*100 << "%) " << sqrt_err_red_tot/tot_red*100<< std::endl; 
  // std::cout <<  "    ttZ_tot = " << tot_ttZ << " +- " << err_ttZ << std::endl; 
  // std::cout <<  "    ttWW_tot = " << tot_ttWW <<" +- " << err_ttWW << std::endl;
  // std::cout <<  "    WWZ_tot = " << tot_WWZ << " +- " << err_WWZ  << std::endl; 
  // std::cout <<  "    tot_irr = " << tot_irr << " +- " << err_irr_tot << "(" <<err_irr_tot/tot_data*100 << "%) " << err_irr_tot/tot_irr*100 <<std::endl; 
  // std::cout <<  "    final data = " << finaldata << " +- " << err_finaldata << std::endl;
  // std::cout <<   "==============================================" << std::endl;  

 output->cd();
 h_data->Write(dataName.c_str(),TObject::kOverwrite);
 h_totdata->Write(TotdataName.c_str(),TObject::kOverwrite);
 output->Close();
 
 output_syst->cd();
 h_data_irrp->Write(dataIrrpName.c_str(),TObject::kOverwrite);
 h_data_irrm->Write(dataIrrmName.c_str(),TObject::kOverwrite);
 h_data_redp->Write(dataRedpName.c_str(),TObject::kOverwrite);
 h_data_redm->Write(dataRedmName.c_str(),TObject::kOverwrite);
 output_syst->Close();

}

// //Build data and data-minus-bkg distributions(for the 4e, 4mu and 2e2mu final states) for the JES systematic uncertainty, where data distributions are shifted up and down by the JES uncertainty.
// void DataToUnfold::Build_JE(string var, string finalstate)
// {
//   output = new TFile((var+"_test/DataToUnfold_JES.root").c_str(), "UPDATE");

//   histoName_up     = "ZZTo" + finalstate + "_" + var+"_JESDataUp_01"; 
//   histoName_down   = "ZZTo" + finalstate + "_" + var+"_JESDataDn_01"; 
//   histoMCName      = "ZZTo" + finalstate + "_" + var+"_01";

//   dataName_up      = "DataminusBkg_"+var+"_ZZTo"+finalstate+"_JESUp";
//   dataName_down    = "DataminusBkg_"+var+"_ZZTo"+finalstate+"_JESDn";
//   TotdataName_up   = "TotData_"+var+"_ZZTo"+finalstate+"_JESUp";
//   TotdataName_down = "TotData_"+var+"_ZZTo"+finalstate+"_JESDn";

//   h_totdata_up     = (TH1*) data->Get(histoName_up.c_str());    
//   h_totdata_down   = (TH1*) data->Get(histoName_down.c_str()); 
//   h_data_up        = (TH1*) h_totdata_up->Clone("h_totdata_up");  
//   h_data_down      = (TH1*) h_totdata_down->Clone("h_totdata_down"); 

//   h_red_up         = (TH1*) red3P1F->Get(histoName_up.c_str());
//   h_red_down       = (TH1*) red3P1F->Get(histoName_down.c_str());


//   h_ttZ = (TH1*) ttZ->Get(histoMCName.c_str());
//   h_WWZ = (TH1*) WWZ->Get(histoMCName.c_str());

 
//   if(h_red_up == NULL)  {cout<<"histo "<<h_red_up->GetName()<<" is Null"<<endl; abort();}
//   if(h_red_down == NULL){cout<<"histo "<<h_red_down->GetName()<<" is Null"<<endl; abort();}

//   float dmb_up = 0;
//   float dmb_down = 0;
  
//   int b = h_red_up->GetNbinsX();
//   for(int i =1; i<=b; i++){
//     dmb_up = 0;
//     dmb_up = h_totdata_up->GetBinContent(i) - h_red_up->GetBinContent(i) - h_ttZ->GetBinContent(i) - h_WWZ->GetBinContent(i);
//     if(dmb_up>0.) h_data_up-> SetBinContent(i,dmb_up);
//     else h_data_up-> SetBinContent(i,0.);
//     dmb_down = 0;
//     dmb_down = h_totdata_down->GetBinContent(i) - h_red_down->GetBinContent(i) - h_ttZ->GetBinContent(i) - h_WWZ->GetBinContent(i);
//     if(dmb_down>0.) h_data_down-> SetBinContent(i,dmb_down);
//     else h_data_down-> SetBinContent(i,0.);

//     // std::cout << "bin " << i << " data up= "<< h_totdata_up->GetBinContent(i) << " +- " << h_totdata_up->GetBinError(i)<< " data down= "<< h_totdata_down->GetBinContent(i) << " +- " << h_totdata_down->GetBinError(i) << " red up = " << h_red_up->GetBinContent(i) << " +- " << h_red_up->GetBinError(i)  << "h_red down" << h_red_down->GetBinContent(i) << " +- " << h_red_down->GetBinError(i)<<   " ttZ = " << h_ttZ->GetBinContent(i) << " +- " << h_ttZ->GetBinError(i) <<  " ttWW = " << h_ttWW->GetBinContent(i) << " +- " << h_ttWW->GetBinError(i) <<   " WWZ = " << h_WWZ->GetBinContent(i)   << " +- " << h_WWZ->GetBinError(i) << " data-bkg up = " << dmb_up <<" +- " << h_data_up->GetBinError(i) << " " << h_data_ip->GetBinContent(i) << " data-bkg down = " << dmb_down <<" +- " << h_data_down->GetBinError(i) << " " << h_data_ip->GetBinContent(i)<< std::endl;
//   }
  
//   output->cd();   
//   h_data_up->Write(dataName_up.c_str(),TObject::kOverwrite);
//   h_data_down->Write(dataName_down.c_str(),TObject::kOverwrite);

//   h_totdata_up->Write(TotdataName_up.c_str(),TObject::kOverwrite); 
//   h_totdata_down->Write(TotdataName_down.c_str(),TObject::kOverwrite);

//   output->Close();
// }

//Plot data distributions
void DataToUnfold::Plot(string var,string finalstate, string path) 
{
  gROOT->Reset();  
  gROOT->SetStyle("Plain");   
  gStyle->SetOptStat(0);
  
  string xAxis; 
  string yAxis = "Events";
  
  if(finalstate == "4m") fs = "4#mu";
  else if(finalstate == "2e2m") fs = "2e2#mu";
  else fs = finalstate;
  
  file =  new TFile((var+"_test/DataToUnfold.root").c_str());
  file_syst = new TFile((var+"_test/DataToUnfold_syst.root").c_str());
 
  dataName = "DataminusBkg_"+var+"_ZZTo"+finalstate;
  TotdataName = "TotData_"+var+"_ZZTo"+finalstate;
  dataIrrpName = "DataminusBkg_irrp_"+var+"_ZZTo"+finalstate; 
  dataIrrmName = "DataminusBkg_irrm_"+var+"_ZZTo"+finalstate; 
  dataRedpName = "DataminusBkg_redp_"+var+"_ZZTo"+finalstate; 
  dataRedmName = "DataminusBkg_redm_"+var+"_ZZTo"+finalstate;
  
  h_data= (TH1*)file->Get(dataName.c_str());
  h_totdata= (TH1*)file->Get(TotdataName.c_str());
  h_data_irrp= (TH1*)file_syst->Get(dataIrrpName.c_str());
  h_data_irrm= (TH1*)file_syst->Get(dataIrrmName.c_str());
  h_data_redp= (TH1*)file_syst->Get(dataRedpName.c_str());
  h_data_redm= (TH1*)file_syst->Get(dataRedmName.c_str());
  
  if(var =="Mass"){
    xAxis = "reco m_{"+fs+"}";
    //max = matrix->GetBinContent(2,2)/2+3;
  }
  else if(var =="nJets"){
    xAxis = "reco Njets";
    // max = matrix->GetBinContent(1,1)/2;
  }
  else if(var =="nIncJets"){
    xAxis = "reco Njets";
    // max = matrix->GetBinContent(1,1)/2;
  }
  else if(var =="Mjj"){
    xAxis = "reco m_{jj}";
    //max = matrix->GetBinContent(2,2)*1.5;
  }
  else if(var =="Deta"){
    xAxis = "reco #Delta#eta_{jj}";
    //max = matrix->GetBinContent(2,2)*1.5;
  }
  else if(var =="nJets_Central"){
    xAxis = "reco Ncentraljets";
    //    max = matrix->GetBinContent(1,1)/3;
  }
  else if(var =="nIncJets_Central"){
    xAxis = "reco Ncentraljets";
    //    max = matrix->GetBinContent(1,1)/3;
  }
  else if(var =="Mjj_Central"){
    xAxis = "reco m_{jj}";
    //    max = matrix->GetBinContent(2,2)*1.5;
  }
  else if(var =="Deta_Central"){
    xAxis = "reco #Delta#eta_{jj}";
    //    max = matrix->GetBinContent(2,2)*1.5;
  }
  else if(var =="PtJet1"){
    xAxis = "reco p_{T}^{jet1}";
    //max = matrix->GetBinContent(1,1);
  }
  else if(var =="PtJet2"){
    xAxis = "reco p_{T}^{jet2}";
    //    max = matrix->GetBinContent(1,1);
  }
 else if(var =="EtaJet1"){
   xAxis = "reco |#eta^{jet1}|";
   //max = matrix->GetBinContent(2,2)/2;
 }
 else if(var =="EtaJet2"){
   xAxis = "reco |#eta^{jet2}|";
   //    max = matrix->GetBinContent(2,2)/2;
 }
   else if(var =="dRZZ"){
   xAxis = "reco #DeltaR(Z_1,Z_2)";
   //    max = matrix->GetBinContent(2,2)/2;
 }
   else if(var =="PtZZ"){
   xAxis = "reco p_{T}^{4#ell}";
   //    max = matrix->GetBinContent(2,2)/2;
 }
  TLegend *leg1 = new TLegend(0.65,0.65,0.6,0.85);
  leg1->SetFillColor(kWhite);
  leg1->SetBorderSize(0);
  leg1->SetTextSize(0.03); 
  
  TCanvas *c1 = new TCanvas ("c1","c1");
  h_totdata->Draw();
  h_totdata->SetLineColor(2);
  h_data->Draw("SAME");
  h_data->SetLineColor(1);
  h_totdata->SetMinimum(0.); 
  
  leg1->AddEntry(h_totdata,"full dataset","l"); 
  leg1->Draw(); 
  leg1->AddEntry(h_data,"data-background","l"); 
  leg1->Draw("SAME");
  
  h_data_redp->GetXaxis()->SetTitle(xAxis.c_str());
  
  TLegend *leg = new TLegend(0.65,0.65,0.6,0.85);
  leg->SetFillColor(kWhite);
  leg->SetBorderSize(0);
  leg->SetTextSize(0.03); 
  
  TCanvas *c = new TCanvas ("c","c");
  
  h_data_redp->Draw();
  h_data_redp->SetMinimum(0.);
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
  
  std::string SavePage = "~/www/PlotsVV/13TeV/";

  string png =SavePage+path+"/"+var+"/"+"Data_"+var+"_ZZTo" + finalstate + ".png";
  string pdf =SavePage+path+"/"+var+"/"+"Data_"+var+"_ZZTo" + finalstate + ".pdf";
  string png_syst =SavePage+path+"/"+var+"/"+"Data_syst_"+var+"_ZZTo" + finalstate + ".png";
  string pdf_syst =SavePage+path+"/"+var+"/"+"Data_syst_"+var+"_ZZTo" + finalstate + ".pdf";
  c->Print(png.c_str());
  c->Print(pdf.c_str());
  c1->Print(png_syst.c_str());
  c1->Print(pdf_syst.c_str());
}


