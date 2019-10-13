/**
 *  Macro to study VZZ cuts to background from results of WZZAnalyzer 
 *
 *	It may become necessary to change the path and/or to expand the samples list	
 *			
 *  $Date: 2019/10/12 09:29:47 
 *  $Revision: 0.1 $
 *
 *  \author C. Tarricone cristiano.tarrico@edu.unito.it
 */

#include <TCanvas.h>
#include <TFile.h>
#include <TH1F.h>
#include <TLegend.h>

#include <boost/foreach.hpp>
#define foreach BOOST_FOREACH

#define TEST_MODE

using namespace std;

void VZZCutEfficiency(){

  //TString otherPath = "~/VVXAnalysis/TreeAnalysis/results/WZZAnalyzer_MC/";
  TString path =  "~/VVXAnalysis/TreeAnalysis/results/WZZAnalyzer_SR/";
  
  vector<TString> sampleNames = {"WZZ", "ZZZ","WZ",/*"ZZTo2e2muJJ","ZZTo4eJJ","ZZTo4muJJ",*/"ggZZ2e2mu","ggZZ4e","ggZZ4mu","ZZTo4lamcatnlo"};
  vector<TString> regionNames = {"sign","bckg"};
  vector<TString> cutNumbers = {"1","2","3","4","5","6"};
  vector<TString> typeNames = {/*"Eta", "Phi", "Pt", "E",*/ "Mass","Tot"};

  
  foreach(TString& sample, sampleNames){
#ifdef TEST_MODE
    cout<<"Opening \""<<sample<<".root\"\n";
#endif

    TFile* result = TFile::Open(path + sample + ".root");

    if(sample=="WZZ" || sample=="ZZZ"){
      foreach(TString& region, regionNames){
#ifdef TEST_MODE
	cout<<"\t"<<region<<"\n";
#endif
	foreach(TString& cutNb, cutNumbers){
#ifdef TEST_MODE
	  cout<<"\t Cut number: "<<cutNb<<"\n";
#endif
    
	  foreach(TString& type, typeNames){
#ifdef TEST_MODE
	    cout<<"\t"<<type<<"\n";
#endif
	    TH1F* hNum = (TH1F*)result->Get("recoV"+type+"_"+region+cutNb);
	    if(hNum == nullptr){
#ifdef TEST_MODE
	      cout<<"Could not open recoV"<<type<<"_"<<region<<cutNb<<"\"\n";
#endif
	    }

	    TH1F* hDen = (TH1F*)result->Get("recoV"+type+"_"+region+"0");			
	    if(hDen == nullptr){
	      #
		cout<<"Could not open recoV"<<type<<"_"<<region<<"0""\"\n";
	    }
			
	    if(hNum != nullptr && hDen != nullptr){
	      TGraphAsymmErrors* hEff = new TGraphAsymmErrors(hNum, hDen, "cp");
	      hEff->SetTitle("afterCut"+cutNb+"_"+type+" ("+sample+" sample)");
	      hEff->GetYaxis()->SetRangeUser(0.,1.01);
	      TCanvas *cDrawing = new TCanvas("VSignal_vs_"+type+sample+region+cutNb,"VSignal_vs_"+type+sample+region+cutNb, 10,0,1280,1024);
	      cDrawing->cd();
	      hEff->Draw("AP");
	      //if(type=="Tot") hEff->Draw("TEXT SAME");

	    }
	  }
	}
      }
    
    }
    else{

      foreach(TString& cutNb, cutNumbers){
#ifdef TEST_MODE
	cout<<"\t Cut number: "<<cutNb<<"\n";
#endif
    
	foreach(TString& type, typeNames){
#ifdef TEST_MODE
	  cout<<"\t"<<type<<"\n";
#endif
	  TH1F* hNum = (TH1F*)result->Get("recoV"+type+"_num"+cutNb);
	  if(hNum == nullptr){
#ifdef TEST_MODE
	    cout<<"Could not open recoV"<<type<<"_num"<<cutNb<<"\"\n";
#endif
	  }

	  TH1F* hDen = (TH1F*)result->Get("recoV"+type+"_den");			
	  if(hDen == nullptr){
	    #
	      cout<<"Could not open recoV"<<type<<"_den""\"\n";
	  }
			
	  if(hNum != nullptr && hDen != nullptr){
	    TGraphAsymmErrors* hEff = new TGraphAsymmErrors(hNum, hDen, "cp");
	    hEff->SetTitle("afterCut"+cutNb+"_"+type+" ("+sample+" sample)");
	    hEff->GetYaxis()->SetRangeUser(0.,1.01);
	    TCanvas *cDrawing = new TCanvas("VSignal_vs_"+type+sample+cutNb,"VSignal_vs_"+type+sample+cutNb, 10,0,1280,1024);
	    cDrawing->cd();
	    hEff->Draw("AP");
	    //if(type=="Tot") hEff->Draw("TEXT SAME");

	  }
	}
      }
    }
    
    result->Close("R");
  
  }
}
