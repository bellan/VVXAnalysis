/**
 *  Macro to study VZZ cuts to background from results of WZZAnalyzer 
 *
 *	It may become necessary to change the path and/or to expand the samples list	
 *			
 *  $Date: 2019/10/12 09:29:47 
 *  $Revision: 0.2 $
 *
 *  \author C. Tarricone cristiano.tarrico@edu.unito.it
 */

#include <TCanvas.h>
#include <TFile.h>
#include <TH1F.h>
#include <TLegend.h>
#include <TMath.h>

#include <boost/foreach.hpp>
#define foreach BOOST_FOREACH

#define TEST_MODE

using namespace std;

void VZZCutEfficiency(){

  TString path =  "~/VVXAnalysis/TreeAnalysis/results/WZZAnalyzer_SR/";
  
  vector<TString> sampleNames = {"WZZ", "ZZZ","WZ","ggZZ2e2mu","ggZZ4e","ggZZ4mu","ZZTo4lamcatnlo"};
  vector<TString> regionNames = {"sign","bckg"};
  vector<TString> cutNumbers = {"0","1","2","3"};
  vector<TString> typeNames = {"Eta", "Phi", "Pt", "E", "Mass","Tot"};
  float compNum[]={0,0,0,0};
  float compDen[]={0,0,0,0};
  
  foreach(TString& sample, sampleNames){
    cout<<"Opening \""<<sample<<".root\"\n";

    TFile* result = TFile::Open(path + sample + ".root");

	int cutCounter=0;

	foreach(TString& cutNb, cutNumbers){
	  cout<<"\t Cut number: "<<cutNb<<"\n";
    
	  foreach(TString& type, typeNames){
	    cout<<"\t"<<type<<"\n";


	    if(sample=="WZZ" || sample=="ZZZ"){
      
	      foreach(TString& region, regionNames){
		cout<<"\t"<<region<<"\n";


	    
		TH1F* hNum = (TH1F*)result->Get("recoV"+type+"_"+region+cutNb);


		if(hNum == nullptr){
		  cout<<"Could not open recoV"<<type<<"_"<<region<<cutNb<<"\"\n";
		}

		TH1F* hDen = (TH1F*)result->Get("recoV"+type+"_"+region+"0");			
		if(hDen == nullptr){
		  cout<<"Could not open recoV"<<type<<"_"<<region<<"0""\"\n";
		}


		if(hNum != nullptr && hDen != nullptr){
		  TGraphAsymmErrors* hEff = new TGraphAsymmErrors(hNum, hDen, "cp");
		  hEff->SetTitle("afterCut"+cutNb+"_"+type+" ("+sample+" sample)");
		  hEff->GetYaxis()->SetRangeUser(0.,1.01);
		  TCanvas *cDrawing = new TCanvas("VSignal_vs_"+type+sample+region+cutNb,"VSignal_vs_"+type+sample+region+cutNb, 10,0,1280,1024);
		    cDrawing->cd();
		    hEff->Draw("AP");
		  
		  if(type=="Tot"){	    
		    if(region=="sign"){
		      compNum[cutCounter]+=hNum->GetBinContent(1);
		      cout<<compNum[cutCounter]<<endl;
		    }
		    else{
		      compDen[cutCounter]+=hNum->GetBinContent(1);
		      cout<<compDen[cutCounter]<<endl;

		    }

		  }
		}
	      }
	    }

	    
	    else{
	      TString region=regionNames.at(1);
	      
	    
	      TH1F* hNum = (TH1F*)result->Get("recoV"+type+"_"+region+cutNb);


	      if(hNum == nullptr){
		cout<<"Could not open recoV"<<type<<"_"<<region<<cutNb<<"\"\n";
	      }

	      TH1F* hDen = (TH1F*)result->Get("recoV"+type+"_"+region+"0");			

	      if(hDen == nullptr){
		cout<<"Could not open recoV"<<type<<"_"<<region<<"0""\"\n";
	      }

	    
	    
	    
	      if(hNum != nullptr && hDen != nullptr){
		TGraphAsymmErrors* hEff = new TGraphAsymmErrors(hNum, hDen, "cp");
		hEff->SetTitle("afterCut"+cutNb+"_"+type+" ("+sample+" sample)");
		hEff->GetYaxis()->SetRangeUser(0.,1.01);
		TCanvas *cDrawing = new TCanvas("VSignal_vs_"+type+sample+region+cutNb,"VSignal_vs_"+type+sample+region+cutNb/*, 10,0,1280,1024*/);
		  cDrawing->cd();
		  hEff->Draw("AP");
		  
		if(type=="Tot"){	    
		    compDen[cutCounter]+=hNum->GetBinContent(1);
		    cout<<compDen[cutCounter]<<endl;

		}
	      }
	    }
	  }
	  cutCounter++;

	}
  result->Close("R");

  }

  TH1F* hComp = new TH1F();  
  for(int i=0; i<cutNumbers.size(); i++){
    compDen[i]=TMath::Sqrt(compDen[i]);
    compNum[i]/=compDen[i];
    cout<<compNum[i]<<endl;
    hComp->SetBinContent(i+1, compNum[i]);
  }

    
}


