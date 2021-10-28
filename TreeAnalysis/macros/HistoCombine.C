///////////////////////////////////////////////////
// Executes the analysis on data samples         //
//                                               //
// Author: E. Racca (eleonora.racca@cern.ch)     //
///////////////////////////////////////////////////

#include "Riostream.h"
#include <TFile.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TString.h>

#include <string>
#include <vector>

using namespace std;

void Histocombine();
void PrintDatacard(TString outpath, TString suffix, TString histofile, vector<TString> namesdc, TString process);
void WZHisto(TString inpath, TString outpath, TString year, TString process);
void ZZHisto(TString inpath, TString outpath, TString year, TString process);


// NOTE: Use this inside macros folder
//       Create the folders ./Combine/datacards and ./Combine/histograms


void HistoCombine(){

  vector<TString> year = {"2016", "2017", "2018", "1618"};
  TString inpath = "\0";
  TString outpath = "Combine/";

  for(int y = 0; y < (int)year.size(); y++){
    cout << "\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << endl;
    cout << "Year: " << year[y] << endl;
    
    inpath = "../results/" + year[y] + "/";
    
    WZHisto(inpath, outpath, year[y], "WZto3l ");
    ZZHisto(inpath, outpath, year[y], "ZZto4l ");
  }

}



void WZHisto(TString inpath, TString outpath, TString year, TString process){

  // Variables declarations
  TFile* filetemp;
  TFile* filedc;
  vector<TFile*> fileout;
  vector<TFile*> samples;
  

  TH1F *histotemp;
  TH1F *histosum;
  vector<TH1F *> histograms;


  // names of the samples
  vector<TString> names = {"WLLJJ", "ZZJJTo4L", "WZTo3LNu", "ZZTo4l", "DYJetsToLL_M50", "TTTo2L2Nu", "TTZJets_M10_MLM", "TTWJetsToLNu", "gg", "data"};

  if(year == "2016")
    names[2] = "ZZTo4lamcatnlo";

  // names of the processes in the datacard
  vector<TString> namesdc = {"WZEW", "ZZEW", "WZQCD", "ZZQCD", "DY", "TT", "TTZ", "TTW", "ggF", "data_obs"};
  vector<TString> namesdcback = {"WZEW", "back", "data_obs"};

  // names of the graphs to combine
  vector<TString> namesgraphs = {"RecoJJ_massvsdeltaEtaabs_AC"};
  vector<int> divisors = {10}; // for rebinning

  // names of the datacards and histograms
  vector<TString> namesfiles = {"JJmasseta"};
  TString histofile = "\0";


  // Open results files
  for(int i = 0; i < (int)names.size(); i++){
    filetemp = TFile::Open(inpath + "WZAnalyzer_SR3L/" + names[i] + ".root");      
    samples.push_back(filetemp);
  }


  // Create and fill files in output
  for(int i = 0; i < (int)namesgraphs.size(); i++){
    cout << "\n-------------------------------------------------------------" << endl;
    cout << "Graph: " << namesgraphs[i] << endl;
  
    fileout.push_back(TFile::Open(outpath + "histograms/histoWZ" + namesfiles[i] + "-" + year + ".root", "RECREATE"));
    fileout[i] -> cd();

    cout << "Creating file" << endl;

    // 1D representation of 2D histograms
    if(namesgraphs[i] == "RecoJJ_massvsdeltaEtaabs_AC"){
      TH1F *mjjdeta = new TH1F("mjjdeta", "mjjdeta", 24, 0, 24);      
      TH2F *temp;
      
      for(int i = 0; i < (int)samples.size(); i++){
	temp = (TH2F*)samples[i]->Get(namesgraphs[i]);
	temp -> RebinX(2);
	
	for(int y = 0; y <= temp->GetNbinsY()-7; y++){
	  for(int x = 1; x <= temp -> GetNbinsX(); x++){
	    mjjdeta -> SetBinContent(x + y*temp->GetNbinsX(), temp->GetBinContent(x, y+7));
	  }
	}
	
	mjjdeta -> SetName(namesdc[i]);
	mjjdeta -> Write(namesdc[i], TObject::kWriteDelete);
      }

    }
    // 1D histograms
    else{
      for(int j = 0; j < (int)samples.size(); j++){
	histograms.push_back((TH1F*)samples[j]->Get(namesgraphs[i]));
	
	// Creating all the histograms
	if(histograms[j] == NULL){
	  cout << "Null histogram in sample: " << names[j] << endl;
	  continue;
	}
	
	histograms[j] -> Rebin(divisors[i]);
	histograms[j] -> SetName(namesdc[j]);
	histograms[j] -> Write(namesdc[j], TObject::kWriteDelete);
	
	
	// Creating sum histogram for datacards with only 1 background
	if(j == 0){
	  histosum = new TH1F(*histograms[j]);
	  histosum -> SetName("back");
	}
	else if(namesdc[j] != "WZEW" && namesdc[j] != "data_obs"){
	  histosum -> Add(histograms[j]);
	}
      }
      
      histosum -> Write("back", TObject::kWriteDelete);
      delete histosum;
      
      for(int j = 0; j < (int)histograms.size(); j++){
	delete histograms[j];
      }
      histograms.clear();
    }

    cout << "Creating datacards" << endl;
    
    histofile = "../VVjjAnalysis/histograms/histoWZ" + namesfiles[i] + "-" + year + ".root";
    
    PrintDatacard(outpath, namesfiles[i] + year + ".txt",  histofile, namesdc, process);
    PrintDatacard(outpath, namesfiles[i] + year + "_back.txt",  histofile, namesdcback, process);
  }
  

  // Closing all open files
  for(int i = 0; i < fileout.size(); i++){
    fileout[i]->Close();
  }
  
  for(int i = 0; i < samples.size(); i++){
    samples[i]->Close();
  }
  
  filetemp->Close();
  samples.clear();
  fileout.clear();
}



void ZZHisto(TString inpath, TString outpath, TString year, TString process){
  
  // Variables declarations
  TFile* filetemp;
  TFile* filedc;
  vector<TFile*> fileout;
  vector<TFile*> samples;
  

  TH1F *histotemp;
  TH1F *histosum;
  vector<TH1F *> histograms;

  
  vector<TString> names = {"WLLJJ", "ZZJJTo4L", "ZZTo4l", "DYTT", "TTZJets_M10_MLM", "TTWJetsToLNu", "gg", "data"};

  if(year == "2016")
    names[2] = "ZZTo4lamcatnlo";

  vector<TString> namesdc = {"WZEW", "ZZEW", "ZZQCD", "DYTT", "TTZ", "TTW", "ggF", "data_obs"};
  vector<TString> namesdcback = {"ZZEW", "back", "data_obs"};

  vector<TString> namesgraphs = {"MELA_discriminant"};
  vector<int> divisors = {10}; // for rebinning
  
  vector<TString> namesfiles = {"MELA"};
  TString histofile = "\0";


  // Open results files
  for(int i = 0; i < (int)names.size(); i++){
    if(names[i] != "DYTT")
      filetemp = TFile::Open(inpath + "ZZjjAnalyzer_SR/" + names[i] + ".root");
    else
      filetemp = TFile::Open(inpath + "ZZjjAnalyzer_CR/data.root");
      
    samples.push_back(filetemp);
  }


  // Create and fill files in output
  for(int i = 0; i < (int)namesgraphs.size(); i++){
    cout << "\n-------------------------------------------------------------" << endl;
    cout << "Graph: " << namesgraphs[i] << endl;
  
    fileout.push_back(TFile::Open(outpath + "histograms/histoZZ" + namesfiles[i] + "-" + year + ".root", "RECREATE"));
    fileout[i] -> cd();

    cout << "Creating file" << endl;

    for(int j = 0; j < (int)samples.size(); j++){
      histograms.push_back((TH1F*)samples[j]->Get(namesgraphs[i]));

      // Creating all the histograms
      if(histograms[j] == NULL){
	cout << "Null histogram in sample: " << names[j] << endl;
	continue;
      }

      histograms[j] -> Rebin(divisors[i]);
      histograms[j] -> SetName(namesdc[j]);
      histograms[j] -> Write(namesdc[j], TObject::kWriteDelete);


      // Creating sum histogram
      if(j == 0){
	histosum = new TH1F(*histograms[j]);
	histosum -> SetName("back");
      }
      else if(namesdc[j] != "WZEW" && namesdc[j] != "ZZEW" && namesdc[j] != "data_obs"){
	histosum -> Add(histograms[j]);
      }
    }

    histosum -> Write("back", TObject::kWriteDelete);
    delete histosum;
    
    for(int j = 0; j < (int)histograms.size(); j++){
      delete histograms[j];
    }
    histograms.clear();

    cout << "Creating datacards" << endl;
    
    histofile = "../VVjjAnalysis/histograms/histoZZ" + namesfiles[i] + "-" + year + ".root";
    
    PrintDatacard(outpath, namesfiles[i] + year + ".txt",  histofile, namesdc, process);
    PrintDatacard(outpath, namesfiles[i] + year + "_back.txt",  histofile, namesdcback, process);
  }
  

  // Closing all open files
  for(int i = 0; i < fileout.size(); i++){
    fileout[i]->Close();
  }
  
  for(int i = 0; i < samples.size(); i++){
    samples[i]->Close();
  }
  
  filetemp->Close();
  samples.clear();
  fileout.clear();
}



void PrintDatacard(TString outpath, TString suffix, TString histofile, vector<TString> namesdc, TString process){

  int nback = (int)namesdc.size()-3; // -3 for WZ, ZZ, data when WZ and ZZ are both signal
  
  ofstream datacardfile;

  datacardfile.open(outpath + "datacards/datacard" + suffix);

  //datacardfile << "imax 1 number of channels" << endl;
  //datacardfile << "jmax " << TString::Itoa(nback, 10) << " number of backgrounds" << endl;
  datacardfile << "imax * number of channels" << endl;
  datacardfile << "jmax * number of backgrounds" << endl;
  datacardfile << "kmax * number of nuisance parameters" << endl << endl;
  datacardfile << "-----------------" << endl << endl;
  datacardfile << "observation -1" << endl << endl;
  datacardfile << "-----------------" << endl << endl;
  datacardfile << "shapes * * " << histofile << " $PROCESS $PROCESS_$SYSTEMATIC" << endl << endl;
  datacardfile << "-----------------" << endl << endl;

  datacardfile << "bin ";
  for(int i = 0; i < (int)namesdc.size()-1; i++){
    datacardfile << process;
  }

  datacardfile << endl << "process ";
  for(int i = 0; i < (int)namesdc.size()-1; i++){
    datacardfile << namesdc[i] << " ";
  }

  datacardfile << endl << "process ";
  for(int i = -1; i < (int)namesdc.size()-2; i++){
    datacardfile << i << " ";
  }

  datacardfile << endl << "rate ";
  for(int i = 0; i < (int)namesdc.size()-1; i++){
    datacardfile << "-1 ";
  }
  
  datacardfile << endl << endl << "-----------------" << endl << endl;
  datacardfile << "#LIST OF SYSTEMATIC UNCERTAINTIES" << endl;

  if(process == "ZZto4l "){
    datacardfile << "QCD_scale               lnN 1.06  1.06  1.10  - -    -    1.09  " << endl;
    datacardfile << "PDF                     lnN 1.066 1.066 1.032 - -    -    1.05  " << endl;
    datacardfile << "CMS_lumi_13TeV          lnN 1.012 1.012 1.012 - -    -    1.012 " << endl;
    datacardfile << "CMS_lep_trig_eff_13TeV  lnN 1.025 1.025 1.025 - -    -    1.025 " << endl;
    datacardfile << "CMS_scale_j_13TeV       lnN 1.007 1.007 1.049 - -    -    1.036 " << endl;
    datacardfile << "CMS_res_j_13TeV         lnN 1.002 1.002 1.022 - -    -    1.01  " << endl;
    datacardfile << "CMS_PU_13TeV            lnN 1.003 1.003 1.002 - -    -    1.004 " << endl;
  
    datacardfile << "CMS_datadriven rateParam " + process + "DYTT 1 [0.67,1.33]" << endl;
    datacardfile << "MC_sample rateParam      " + process + "TTZ  1 [0.81,1.19]" << endl;
    datacardfile << "MC_sample rateParam      " + process + "TTW  1 [0.81,1.19]" << endl;
  }
  else if(process == "WZto3l "){
    datacardfile << "QCD_scale               lnN 1.08  1.015 1.147 1.015 - -  -  -    1.09  " << endl;
    datacardfile << "PDF                     lnN 1.078 1.010 1.044 1.010 - -  -  -    1.05  " << endl;
    datacardfile << "CMS_lumi_13TeV          lnN 1.012 1.012 1.012 1.012 - -  -  -    1.012 " << endl;
    datacardfile << "CMS_lep_trig_eff_13TeV  lnN 1.022 1.027 1.022 1.027 - -  -  -    1.025 " << endl;
    datacardfile << "CMS_scale_j_13TeV       lnN 1.03  1.02  1.047 1.02  - -  -  -    1.036 " << endl;
    datacardfile << "CMS_res_j_13TeV         lnN 1.01  1.041 1.01  1.041 - -  -  -    1.01  " << endl;
    datacardfile << "CMS_PU_13TeV            lnN 1.002 1.002 1.002 1.002 - -  -  -    1.002 " << endl;
  
    datacardfile << "MC_sample rateParam " + process + "DY  1 [0.80,1.20]" << endl;
    datacardfile << "MC_sample rateParam " + process + "TT  1 [0.80,1.20]" << endl;
    datacardfile << "MC_sample rateParam " + process + "TTZ 1 [0.80,1.20]" << endl;
    datacardfile << "MC_sample rateParam " + process + "TTW 1 [0.80,1.20]" << endl;
  }
  
  datacardfile.close();
}
