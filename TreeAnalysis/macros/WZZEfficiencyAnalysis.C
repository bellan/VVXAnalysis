/**
 *  Macro to get efficiency from results of WZZAnalyzer 
 *		
 *  $Date: 2019/09/12 9:47:08 
 *  $Revision: 0.2 $
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

void WZZEfficiencyAnalysis(){

  vector<TString> parNames = {"Electrons", "Muons"/*, "Jets"*/};
  vector<TString> typeNames = {"Eta", "Phi", "Pt", "E"};

  TFile* result = TFile::Open("~/VVXAnalysis/TreeAnalysis/results/WZZAnalyzer_MC/WZZ.root");


  foreach(TString& name, parNames){
    #ifdef TEST_MODE
    cout<<name<<"\n";
    #endif
    foreach(TString& type, typeNames){
      #ifdef TEST_MODE
      cout<<"\t"<<type<<"\n";
      #endif
      TH1F* hNum = (TH1F*)result->Get("gen"+name+type+"_num");
      if(hNum == nullptr){
        #ifdef TEST_MODE
	cout<<"Could not open \gen"<<name<<type<<"_num""\"\n";
#endif
      }

      TH1F* hDen = (TH1F*)result->Get("gen"+name+type+"_den");
			
      if(hDen == nullptr){
	#
	cout<<"Could not open \gen"<<name<<type<<"_den""\"\n";
      }
			
      if(hNum != nullptr && hDen != nullptr){
	TGraphAsymmErrors* hEff = new TGraphAsymmErrors(hNum, hDen, "cp");
	hEff->SetTitle(name+"Efficiency_vs_"+type);
	hEff->GetYaxis()->SetRangeUser(0.,1.01);
	TCanvas *cDrawing = new TCanvas(name+"Efficiency_vs_"+type, name+"Efficiency_vs_"+type, 10,0,1280,1024);
	cDrawing->cd();
	hEff->Draw("AP");
      }
    }
  }

  result->Close("R");
}
