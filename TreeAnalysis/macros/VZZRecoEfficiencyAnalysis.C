/**
 *  Macro to get efficiency of Bosons reconstruction from results of WZZAnalyzer 
 *
 *
 *	Usage: 	root [-l] [-b] [-q] 'VZZRecoEfficiencyAnalysis.C("<sample name>")'		-l do not display banner 	-b run in background		-q close after finishing
 *	e.g.: 	root -l 'VZZRecoEfficiencyAnalysis.C("WZZ")'
 *	It may become necessary to change the path and/or to expand the samples list	
 *			
 *  $Date: 2019/09/29 23:33:08 
 *  $Revision: 2.0 $
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

void VZZRecoEfficiencyAnalysis(TString requestedSample){

  //gStyle->SetOptStat(0);
  
  TString path = "~/VVXAnalysis/TreeAnalysis/results/WZZAnalyzer_MC/";
  TString otherPath =  "~/VVXAnalysis/TreeAnalysis/results/WZZAnalyzer_SR/";
  
  vector<TString> samples = {"WZZ", "ZZZ"};
  vector<TString> parNames = {"W","Z"};
  vector<TString> typeNames = {"Eta", "Phi", "Pt", "E", "Mass","Tot"};
  vector<TString> algNames = {"mW", "mZ", "maxVPt", "minTotPt", "mWZ"};

  
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

  /*foreach(TString& name, parNames){
    #ifdef TEST_MODE
    #endif
  */


  foreach(TString& alg, algNames){
#ifdef TEST_MODE
    cout<<"\t"<<alg<<"\n";
#endif
    foreach(TString& type, typeNames){
#ifdef TEST_MODE
      cout<<"\t"<<type<<"\n";
#endif
      TH1F* hNum = (TH1F*)result->Get("gen"+name+type+"_"+sampleName+"_"+alg+"_num");
      if(hNum == nullptr){
#ifdef TEST_MODE
	cout<<"Could not open gen"<<name<<type<<"_"<<sampleName<<"_"<<alg<<"_num""\"\n";
#endif
      }

      TH1F* hDen = (TH1F*)result->Get("gen"+name+type+"_"+sampleName+"_"+alg+"_den");			
      if(hDen == nullptr){
	#
	  cout<<"Could not open gen"<<name<<type<<"_"<<sampleName<<"_"<<alg<<"_den""\"\n";
      }
			
      if(hNum != nullptr && hDen != nullptr){
	TGraphAsymmErrors* hEff = new TGraphAsymmErrors(hNum, hDen, "cp");
	hEff->SetTitle(name+"Efficiency_vs_"+type+" ("+sampleName+" sample, algorithm "+alg+")");
	hEff->GetYaxis()->SetRangeUser(0.,1.01);
	TCanvas *cDrawing = new TCanvas(name+"ReconstructionEfficiency_vs_"+type+"_"+alg, name+"ReconstructionEfficiency_vs_"+type+"_"+alg, 10,0,1280,1024);
	cDrawing->cd();
	hEff->Draw("AP");
	//if(type=="Tot") hEff->Draw("TEXT SAME");

      }
    }
  }

  result->Close("R");
  
}
