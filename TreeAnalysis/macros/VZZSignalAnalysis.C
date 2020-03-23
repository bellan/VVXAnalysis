/**
 *  Macro to study VZZ background in signal samples "WZZ" and "ZZZ" from results of WZZAnalyzer 
 *
 *
 *	Usage: 	root [-l] [-b] [-q] 'VZZSignalAnalysis.C("<sample name>")'		-l do not display banner 	-b run in background		-q close after finishing
 *	e.g.: 	root -l 'VZZSignalAnalysis.C("WZZ")'
 *	It may become necessary to change the path and/or to expand the samples list	
 *			
 *  $Date: 2019/10/12 08:45:26 
 *  $Revision: 0.3 $
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

void VZZSignalAnalysis(TString requestedSample){

  TString otherPath = "~/VVXAnalysis/TreeAnalysis/results/WZZAnalyzer_MC/";
  TString path =  "~/VVXAnalysis/TreeAnalysis/results/WZZAnalyzer_SR/";
  
  vector<TString> samples = {"WZZ", "ZZZ"/*,"ZZTo2e2muJJ","ZZTo4eJJ","ZZTo4muJJ"*/};
  vector<TString> parNames = {"W","Z"/*,"V"*/};
  vector<TString> typeNames = {/*"Eta", "Phi", "Pt", "E",*/ "Mass","Tot"};
  vector<TString> regionNames = {"sign", "bckg"};
  
  TString sampleName, name;
  if(requestedSample != "")
    for(int i = 0; i < samples.size(); i++)
      if(requestedSample == samples.at(i)){
	sampleName = requestedSample;
	name=parNames.at(i);
	cout<<name<<"\n";
      }

  if(sampleName == ""){
    cout<<"Unknown sample \""<<sampleName<<"\"\n";
    return;
  }
  
  cout<<"Opening \""<<sampleName<<".root\"\n";

  TFile* result = TFile::Open(path + sampleName + ".root");

  foreach(TString& region, regionNames){
#ifdef TEST_MODE
    cout<<"\t"<<region<<"\n";
#endif

  
    foreach(TString& type, typeNames){
#ifdef TEST_MODE
      cout<<"\t"<<type<<"\n";
#endif
      TH1F* hNum = (TH1F*)result->Get("genV"+type+"_"+sampleName+"_num_"+region);
      if(hNum == nullptr){
#ifdef TEST_MODE
	cout<<"Could not open genV"<<type<<"_"<<sampleName<<"_num_"<<region<<"\"\n";
#endif
      }

      TH1F* hDen = (TH1F*)result->Get("genV"+type+"_"+sampleName+"_den");			
      if(hDen == nullptr){
	#
	  cout<<"Could not open genV"<<type<<"_"<<sampleName<<"_den""\"\n";
      }
			
      if(hNum != nullptr && hDen != nullptr){
	TGraphAsymmErrors* hEff = new TGraphAsymmErrors(hNum, hDen, "cp");
	hEff->SetTitle(region+"_"+type+" ("+sampleName+" sample)");
	hEff->GetYaxis()->SetRangeUser(0.,1.01);
	TCanvas *cDrawing = new TCanvas(region+"_"+type,region+"_"+type, 10,0,1280,1024);
	cDrawing->cd();
	hEff->Draw("AP");
	//if(type=="Tot") hEff->Draw("TEXT SAME");

      }
    }
  }
  result->Close("R");
  
}
