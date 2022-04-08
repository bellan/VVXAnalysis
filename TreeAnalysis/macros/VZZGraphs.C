/**
 *	Macro to study cuts and variable distributions from results of any Analyzer
 *	
 *	It may become necessary to change the path and/or signal and background	
 *	
 *	Usage: 	root [OPTIONS] 'VZZGraph.C("<category of graphs>", "<type of graphs>", "<requested graphs by name>")'
 *	e.g.: 	root -l 'VZZGraph.C("jetlep", "13", "all")'
 *
 *	It's possible to use it with standard values, by simply calling:	root -l WWosCutAnalysis.C
 *	
 *	<category of graphs> = {"cut", "jet", "lep"}
 *	<type of graphs> = {
 *		"0" sig and bkg same canvas (no stack)
 *		"1" sig and bkg same canvas (stacked)
 *		"2" sig, bkg and data same canvas (stacked) // Integral graps on same canvas (stacked)
 *		"3" sig/sqrt(bkg) (integrals)
 *	}
 *	<requested graphs by name> depends on the chosen category of graphs
 *	
 *  $Date: 2018/09/13 13:37:23 $
 *  $Revision: 1.0 $
 *  \author A. Mecca alberto.mecca@edu.unito.it
 */

#include "TCanvas.h"
#include "TFile.h"
#include "TKey.h"
#include "TROOT.h"
#include "TH1F.h"
#include "TLegend.h"
#include "THStack.h"
#include "TLine.h"
#include <iostream>
#include <bitset>
#include <string>
#include <algorithm>

#include <boost/foreach.hpp>
#define foreach BOOST_FOREACH

#define R_PATTERN "_r_" //If present in a graph's name, it means that a "reverse integral of sig/sqrt(bkg)" must be done
#define F_PATTERN "_f_" // ... forward integral of sig/sqrt(bkg)
//Obsolete
//#define NO_GRAPH "_N_" //If present in a graph's name, it means that it must be ignored
#define TEST_MODE

using std::cout;
using std::vector;
using std::bitset;
using std::string;
using std::pair;

TString nullTString = TString("");

struct GNames{
public:
  GNames(const TString& n, const TString& t, const vector<TString>& s,const vector<TString>& b) : graphName(n), graphTitle(t)
  {
    for(auto & si : s)
      sigGrNames.push_back( si.First('_') > 0 ? si(0, si.First('_')) : si );
    for(auto & bi : b)
      bkgGrNames.push_back( bi.First('_') > 0 ? bi(0, bi.First('_')) : bi );
  };
  TString graphName;   //Type of graph ("Cuts", "leadL pt", ...)
  TString graphTitle;  //Prettyfication
  vector <TString> sigGrNames;
  vector <TString> bkgGrNames;
  TString toString() const{
    TString res = "graphName: "+graphName+"\nSIG: ";
    foreach(const TString& s, sigGrNames) res += s+" \t";
    res += "\nBKG: ";
    foreach(const TString& b, bkgGrNames) res += b+" \t";
    return res;
  }
};

template <class TH=TH1F>
struct MyGraphs{
  //MyGraphs(vector<TH*> a, vector<TH*> b, TH* D): vhSigs(a), vhBkgs(b), data(D) {};
  MyGraphs(vector<TH*> a, vector<TH*> b, vector<TH*> c, vector<TH*> d, TH* D): vhSigs(a), vhBkgs(b), vhSigInt(c), vhBkgInt(d), data(D) {};
  vector<TH*> vhSigs;
  vector<TH*> vhBkgs;
  vector<TH*> vhSigInt;
  vector<TH*> vhBkgInt;
  TH* data;
  /*TString toString() const{
    return TString("vhSigs.size(): ")+vhSigs.size()+" - vhBkgs.size(): "+vhBkgs.size();
    }*/
};

TFile* openRootFile(const TString&/*path*/,const TString&/*name*/);
template <class TH1Z = TH1F>
TH1Z* openGraph(TFile*, const TString&, const TString& = nullTString);  // Points to the graph in the TFile
TH1F* integralGraph(const TH1F* const origin);  // Creates a copy which is modified
TH1F* reverseIntegrGraph(const TH1F* const origin);  // Integrates from the last bin
TH1F* sqrtGraph(const TH1F* const origin);  // Creates a copy which is modified
bool doThisGraphName(TString graphName, TString reqGrName);

template <size_t N, class TH=TH1F>
void printGraph(const MyGraphs<TH>& theGraphs, const GNames& theNames, const bitset<N>& typeToDo, const string&);

// Specific types of graph
template <class TH>
void compareCutGraphs(const MyGraphs<TH>&, const GNames&); //Uses helper functions
template <class TH>
void printGraphSame(const MyGraphs<TH>&, const GNames&, const string&);
template <class TH>
void printGraphStack(const MyGraphs<TH>&, const GNames&, const string&);
template <class TH>
void printGraphRatio(const MyGraphs<TH>&, const GNames&, const string& );
template <class TH>
void printGraphSqrt(const MyGraphs<TH>&, const GNames&, const string&);

//Struct building
template <class TH = TH1F>
MyGraphs<TH>* buildMyGraphs(const GNames& names, vector<TFile*>& signalFiles, vector<TFile*>& backgroundFiles, TFile* dataFile);

static const vector<Color_t> myColors =     {kRed,kYellow+8,kGreen-7,kOrange-3,kAzure+10,kBlue,  kViolet+1,kYellow  ,kYellow+2,kOrange+3,kGray+3, kMagenta+3, kGreen+2, kViolet+8, kOrange-4, kGreen+8, kBlack, kMagenta-5, kBlue+6, kRed-3, kAzure-4, kViolet-5};
static const vector<Color_t> myFillColors = {kRed,kYellow+8,kGreen-7,kOrange-3,kAzure+10,kBlue,  kViolet+1,kYellow  ,kYellow+2,kOrange+3,kGray+3, kMagenta+3, kGreen+2, kViolet+8, kOrange-4, kGreen+8, kBlack, kMagenta-5, kBlue+6, kRed-3, kAzure-4, kViolet-5};
//vector<Color_t> myFillColors(myColors.size(), 0);

const float  _FLIP = 1;
const std::vector<std::string> regions{"SR4P", "CR3P1F", "CR110", "CR101", "CR011", "CR100", "CR001", "CR010", "CR000"};

void VZZGraphs(string sReqCateg = string(""), string sReqType = "", string reqGrName = string("")){
  TH1::SetDefaultSumw2(true);
  
  for (auto & region : regions)   
    {
      cout << "\t----- Region: " << region << "-----\n";
      TString path = Form("../results/2016/VVGammaAnalyzer_%s/", region.c_str());
  
      vector<TString> signals = {"ZZGTo4LG"};
      vector<TString> backgrounds = {"ZZTo4l", "WZTo3LNu","ggTo4l", "ZGToLLG", "WGToLNuG", "triboson", "DYJetsToLL_M50", "WJetsToLNu", "TTWW", "TTZZ", "TTTo2L2Nu","tW","tZq","singleT","TTZJets_APV", "TTZToLLNuNu_M10" }; //add WZ2lnu and Wjets  //"ggTo2e2mu_Contin_MCFM701, "ZZZ", "WWZ", "WZZ", "ggTo4e_Contin_MCFM701","ggTo4mu_Contin_MCFM701
  

      std::transform(reqGrName.begin(), reqGrName.end(), reqGrName.begin(), ::tolower);
      std::transform(sReqCateg.begin(), sReqCateg.end(), sReqCateg.begin(), ::tolower);
  
      vector<std::pair<TString,TString>> variablePlots = {};

      bool fourlep_region  = false;
      bool threelep_region = false;

      if (region=="CR3P1F" || region=="CR2P2F" || region=="SR4P")
	fourlep_region = true;
      else if (region=="CR110" || region=="CR101" || region=="CR011" || region=="CR100" || region=="CR001" || region=="CR010" || region=="CR000" || region=="SR3P")
	threelep_region = true;

      if ( fourlep_region )
	{
	  variablePlots.push_back({"ZZ_mass_4e0m", "m_{4l} ;GeV/c^{2}"});
	  variablePlots.push_back({"ZZ_mass_2e2m", "m_{4l} ;GeV/c^{2}"});
	  variablePlots.push_back({"ZZ_mass_0e4m", "m_{4l} ;GeV/c^{2}"});
	  variablePlots.push_back({"ZZ_pT_4e0m",   "pT_{ZZ} ;GeV/c"});
	  variablePlots.push_back({"ZZ_pT_2e2m",   "pT_{ZZ} ;GeV/c"});
	  variablePlots.push_back({"ZZ_pT_0e4m",   "pT_{ZZ} ;GeV/c"});
	  // variablePlots.push_back({"2e2mu_mass", "m_{2e2mu} ;GeV/c^{2}"});
	  // variablePlots.push_back({"4e_mass"   , "m_{4e} ;GeV/c^{2}"});
	  // variablePlots.push_back({"4mu_mass"  , "m_{4mu} ;GeV/c^{2}"});
	  // variablePlots.push_back({"2e2mu_pT",   "pT_{2e2mu} ;GeV/c"});
	  // variablePlots.push_back({"4e_pT",   "pT_{4e} ;GeV/c"});
	  // variablePlots.push_back({"4mu_pT",   "pT_{4mu} ;GeV/c"});

	  variablePlots.push_back({"Z0_mass_4e0m", "m_{Z0} ;GeV/c^{2}"});
	  variablePlots.push_back({"Z0_mass_2e2m", "m_{Z0} ;GeV/c^{2}"});
	  variablePlots.push_back({"Z0_mass_0e4m", "m_{Z0} ;GeV/c^{2}"});
	  variablePlots.push_back({"Z1_mass_4e0m", "m_{Z1} ;GeV/c^{2}"});
	  variablePlots.push_back({"Z1_mass_2e2m", "m_{Z1} ;GeV/c^{2}"});
	  variablePlots.push_back({"Z1_mass_0e4m", "m_{Z1} ;GeV/c^{2}"});
	  variablePlots.push_back({"Z0_l0_pt_4e0m", "pT_{l0,Z0} ;GeV/c"});
	  variablePlots.push_back({"Z0_l0_pt_2e2m", "pT_{l0,Z0} ;GeV/c"});
	  variablePlots.push_back({"Z0_l0_pt_0e4m", "pT_{l0,Z0} ;GeV/c"});
	  variablePlots.push_back({"Z0_l1_pt_4e0m", "pT_{l1,Z0} ;GeV/c"});
	  variablePlots.push_back({"Z0_l1_pt_2e2m", "pT_{l1,Z0} ;GeV/c"});
	  variablePlots.push_back({"Z0_l1_pt_0e4m", "pT_{l1,Z0} ;GeV/c"});
	  variablePlots.push_back({"Z1_l0_pt_4e0m", "pT_{l0,Z1} ;GeV/c"});
	  variablePlots.push_back({"Z1_l0_pt_2e2m", "pT_{l0,Z1} ;GeV/c"});
	  variablePlots.push_back({"Z1_l0_pt_0e4m", "pT_{l0,Z1} ;GeV/c"});
	  variablePlots.push_back({"Z1_l1_pt_4e0m", "pT_{l1,Z1} ;GeV/c"});
	  variablePlots.push_back({"Z1_l1_pt_2e2m", "pT_{l1,Z1} ;GeV/c"});
	  variablePlots.push_back({"Z1_l1_pt_0e4m", "pT_{l1,Z1} ;GeV/c"});
	}


      else if( threelep_region )
	{
	  variablePlots.push_back({"ZW_massT_3e0m", "m_{T,3l}; GeV/c^{2}"});
	  variablePlots.push_back({"ZW_massT_2e1m", "m_{T,3l}; GeV/c^{2}"});
	  variablePlots.push_back({"ZW_massT_1e2m", "m_{T,3l}; GeV/c^{2}"});
	  variablePlots.push_back({"ZW_massT_0e3m", "m_{T,3l}; GeV/c^{2}"});
	  variablePlots.push_back({"Z_mass_3e0m",   "m_{Z}; GeV/c^{2}"   });
	  variablePlots.push_back({"Z_mass_2e1m",   "m_{Z}; GeV/c^{2}"   });
	  variablePlots.push_back({"Z_mass_1e2m",   "m_{Z}; GeV/c^{2}"   });
	  variablePlots.push_back({"Z_mass_0e3m",   "m_{Z}; GeV/c^{2}"   });
	  variablePlots.push_back({"W_massT_3e0m",  "m_{T,W}; GeV/c^{2}" });
	  variablePlots.push_back({"W_massT_2e1m",  "m_{T,W}; GeV/c^{2}" });
	  variablePlots.push_back({"W_massT_1e2m",  "m_{T,W}; GeV/c^{2}" });
	  variablePlots.push_back({"W_massT_0e3m",  "m_{T,W}; GeV/c^{2}" });
	  variablePlots.push_back({"Z_l0_pt_3e0m",  "pT_{l0,Z} ;GeV/c"   });
	  variablePlots.push_back({"Z_l0_pt_2e1m",  "pT_{l0,Z} ;GeV/c"   });
	  variablePlots.push_back({"Z_l0_pt_1e2m",  "pT_{l0,Z} ;GeV/c"   });
	  variablePlots.push_back({"Z_l0_pt_0e3m",  "pT_{l0,Z} ;GeV/c"   });
	  variablePlots.push_back({"Z_l1_pt_3e0m",  "pT_{l1,Z} ;GeV/c"   });
	  variablePlots.push_back({"Z_l1_pt_2e1m",  "pT_{l1,Z} ;GeV/c"   });
	  variablePlots.push_back({"Z_l1_pt_1e2m",  "pT_{l1,Z} ;GeV/c"   });
	  variablePlots.push_back({"Z_l1_pt_0e3m",  "pT_{l1,Z} ;GeV/c"   });
	  variablePlots.push_back({"W_l_pt_3e0m",   "pT_{l,W} ;GeV/c"    });
	  variablePlots.push_back({"W_l_pt_2e1m",   "pT_{l,W} ;GeV/c"    });
	  variablePlots.push_back({"W_l_pt_1e2m",   "pT_{l,W} ;GeV/c"    });
	  variablePlots.push_back({"W_l_pt_0e3m",   "pT_{l,W} ;GeV/c"    });
	  variablePlots.push_back({"W_MET_pt_3e0m", "MET ;GeV/c"         });
	  variablePlots.push_back({"W_MET_pt_2e1m", "MET ;GeV/c"         });
	  variablePlots.push_back({"W_MET_pt_1e2m", "MET ;GeV/c"         });
	  variablePlots.push_back({"W_MET_pt_0e3m", "MET ;GeV/c"         });
	}

      else if( region=="CRLFR" )
	{
	  variablePlots.push_back({"ZL_mass_3e0m", "m_{3l} ;GeV/c^{2}"});
	  variablePlots.push_back({"ZL_mass_2e1m", "m_{3l} ;GeV/c^{2}"});
	  variablePlots.push_back({"ZL_mass_1e2m", "m_{3l} ;GeV/c^{2}"});
	  variablePlots.push_back({"ZL_mass_0e3m", "m_{3l} ;GeV/c^{2}"});
	  variablePlots.push_back({"Z_mass_3e0m",  "m_{Z} ;GeV/c^{2}" });
	  variablePlots.push_back({"Z_mass_2e1m",  "m_{Z} ;GeV/c^{2}" });
	  variablePlots.push_back({"Z_mass_1e2m",  "m_{Z} ;GeV/c^{2}" });
	  variablePlots.push_back({"Z_mass_0e3m",  "m_{Z} ;GeV/c^{2}" });
	  variablePlots.push_back({"Z_l0_pt_3e0m", "pT_{l0,Z} ;GeV/c" });
	  variablePlots.push_back({"Z_l0_pt_2e1m", "pT_{l0,Z} ;GeV/c" });
	  variablePlots.push_back({"Z_l0_pt_1e2m", "pT_{l0,Z} ;GeV/c" });
	  variablePlots.push_back({"Z_l0_pt_0e3m", "pT_{l0,Z} ;GeV/c" });
	  variablePlots.push_back({"Z_l1_pt_3e0m", "pT_{l1,Z} ;GeV/c" });
	  variablePlots.push_back({"Z_l1_pt_2e1m", "pT_{l1,Z} ;GeV/c" });
	  variablePlots.push_back({"Z_l1_pt_1e2m", "pT_{l1,Z} ;GeV/c" });
	  variablePlots.push_back({"Z_l1_pt_0e3m", "pT_{l1,Z} ;GeV/c" });
	  variablePlots.push_back({"L_pt_3e0m",    "pT_{l} ;GeV/c"    });
	  variablePlots.push_back({"L_pt_2e1m",    "pT_{l} ;GeV/c"    });
	  variablePlots.push_back({"L_pt_1e2m",    "pT_{l} ;GeV/c"    });
	  variablePlots.push_back({"L_pt_0e3m",    "pT_{l} ;GeV/c"    });
	}
      //variablePlots.push_back({"AAA_cuts_u"," cutflow u"});
      variablePlots.push_back({"AAA_cuts"," cutflow weighted"});
  
      //#####end graph names#####
      TFile* dataFile = openRootFile(path, "data");
  
      vector<TFile*> signalFiles;
      for(size_t i = 0; i < signals.size(); i++){
	//cout<<"signals.at("<<i<<") = "<<signals.at(i)<<'\n';
	TFile* f = openRootFile(path, signals.at(i));
	//cout<<" \tnullptr? "<<(f == nullptr);
	if(f != nullptr)
	  signalFiles.push_back( f );
      }

      vector<TFile*> backgroundFiles;
      for(size_t i = 0; i < backgrounds.size(); i++){
	//cout<<"backgrounds.at("<<i<<") = "<<backgrounds.at(i)<<'\n';
	TFile* f = openRootFile(path, backgrounds.at(i));
	//cout<<" \tnullptr? "<<(backgroundFiles.at(i) == nullptr);
	if(f != nullptr)
	  backgroundFiles.push_back(f);
      }
  
#ifdef TEST_MODE
      cout<<"\nsignals: (size = "<<signals.size()<<')';
      foreach(TString& b, signals) cout<<"\n\t"<<b;
      cout<<"\nsignalFiles: (size = "<<signalFiles.size()<<')';
      foreach(TFile*& f, signalFiles) cout<<"\n\tfile found";
      cout<<'\n';
  
      cout<<"\nbackgrounds: (size = "<<backgrounds.size()<<')';
      foreach(TString& b, backgrounds) cout<<"\n\t"<<b;
      cout<<"\nbackgroundFiles: (size = "<<backgroundFiles.size()<<')';
      foreach(TFile*& f, backgroundFiles) cout<<"\n\tfile found";
      cout<<'\n';
#endif
  
      if(signalFiles.size() == 0){
	cout<<"No signal file found.\n";
	// exit(1);
      }
  
      cout<<'\n';
      /*
	TIter keyList(signalFiles.at(0)->GetListOfKeys());
	TKey *key;
	while ((key = (TKey*)keyList())) {
	TClass *cl = gROOT->GetClass(key->GetClassName());
	if (!cl->InheritsFrom("TH1")) continue;
	TString name ( key->ReadObj()->GetName() );
	if(name.Contains("weight_")) continue;
	variablePlots.push_back(name);//cout<<name<<'\n';
	}*/

      // ---------- Input interptretation ----------
      bitset<3> categToDo(string("10"));
      if(sReqCateg.find("all") != string::npos) categToDo.set();  //All bits to 1
      if(sReqCateg.find("cut") != string::npos) categToDo.set(0); //set to 1
      if(sReqCateg.find("var") != string::npos) categToDo.set(1);
  
      bitset<4> typeToDo;
      cout<<"sReqType = \""<<sReqType<<"\"\t";
      if(sReqType == string("")){ //Default
	cout<<"Setting default\n";
	//typeToDo.set(0);
	//typeToDo.set(1); //sig and bkg same canvas (stacked)
	typeToDo.set(2);   //stack + ratio
	//typeToDo.set(3); //sig/sqrt(bkg) (integrals)
      }
      else{
	cout<<'\n';
	if(sReqType.find('0') != string::npos) typeToDo.set(0); //sig and bkg no stack
	if(sReqType.find('1') != string::npos) typeToDo.set(1); //Stack
	if(sReqType.find('2') != string::npos) typeToDo.set(2); //Stack + data + ratio
	if(sReqType.find('3') != string::npos) typeToDo.set(3); //sig/sqrt(bkg) (integrals)
      }
  
      cout<<"categToDo: "<<categToDo<<" \ttypeToDo: "<<typeToDo<<'\n';

      TFile* outFile = nullptr; //Dummy file, created to be the currente dir of the histograms, so that I can close the source files
      outFile = TFile::Open(path + "graphs.root", "RECREATE");

      cout<<'\n';

      //################################EVENT PLOTS ################################
      // if(categToDo.test(0)){

      // 	//#####Struct creation #####
      // 	GNames theNames("AAA cuts w", "Cuts", signals, backgrounds);
      // 	MyGraphs<TH1F> theGraphs = *buildMyGraphs(theNames, signalFiles, backgroundFiles, dataFile);
      // 	//#####--------------- #####

      // 	compareCutGraphs(theGraphs, theNames);//Integration makes no sense here
      // }


      //################################VARIABLE PLOTS ################################
      if(categToDo.test(1)){
	for(std::pair<TString, TString> nameTitle : variablePlots){
#ifdef TEST_MODE
	  cout<<"Varible plots: looping over \""<<nameTitle.first<<"\"\n";
#endif
	  if(! doThisGraphName(nameTitle.first, reqGrName) ) continue;

	  //#####Struct creation #####
	  GNames theNames(nameTitle.first, nameTitle.second, signals, backgrounds);
	  MyGraphs<TH1F> theGraphs = *buildMyGraphs(theNames, signalFiles, backgroundFiles, dataFile);
	  //#####--------------- #####

	  printGraph(theGraphs, theNames, typeToDo, region);
	  /*
	    if(graphName == "mll"){
	    int finalCut = (int)(420.)/10;
	    double sigError = 0.;
	    double sigIntegral = theGraphs.hSig->IntegralAndError(finalCut,-1,sigError);
	    cout<<"\n\nfinalCut :"<<finalCut*10<<'\n'<<signals<<": \t"<< sigIntegral;
	    cout<<" \t+- "<<sigError<<'\n';

	    double integral, error;
	    double totBkg = 0.;
	    double varBkg = 0.;
	    for(size_t k = 0; k < theGraphs.vhBkgs.size(); ++k){
	    integral = theGraphs.vhBkgs.at(k)->IntegralAndError(finalCut,-1,error);
	    cout<<backgrounds.at(k)<<": "<< integral;
	    cout<<" \t+- "<<error<<'\n';
	    totBkg += integral;
	    varBkg += error*error;
	    }
	    cout<<"\nTotal Background: "<<totBkg<<" +- "<<sqrt(varBkg);
	    }*/
	}
	cout<<'\n';
      }

      cout<<'\n';

      //outFile->Close();
      for(size_t i = 0; i < signalFiles.size(); i++)
	signalFiles.at(i)->Close();
      for(size_t i = 0; i < backgroundFiles.size(); i++)
	backgroundFiles.at(i)->Close();
    }


}//end loop regions


// ################################################################################################
TFile* openRootFile(const TString& path, const TString& name){
#ifdef TEST_MODE
  cout<<"\nOpening \""<<name<<".root\" \t";
#endif
  TFile* file = nullptr;
  file = TFile::Open(path + name + ".root", "READ");
  if(file == nullptr){
    cout<<"Error: file is nullptr\n";
    return nullptr;
  }
  if(file->IsZombie()){
    cout<<"\""<<name<<".root\" does not exist\n";
    file->Close();
    delete file;
  }
  return file;
}

template <class TH1Z = TH1F>
TH1Z* openGraph(TFile* file, const TString& graphName, const TString& fileName){
  if(file == nullptr) return nullptr;
  TH1Z* hSig = (TH1Z*)file->Get(graphName); 
#ifdef TEST_MODE
  if(hSig == nullptr) 
    cout<<"Could not open \""<<graphName<<"\" from \""<<fileName<<".root"<<"\"\n";
  else{
    cout<<"Retrieved \""<<graphName<<"\" from \""<<fileName<<".root\""<<"\tEntries: "<<hSig->GetEntries()<<'\n';
  }
#endif
  return hSig;
}


// ########################## Utils ####################################

TH1F* integralGraph(const TH1F* const origin){
  if(origin == nullptr){
    cout<<" \tError: nullptr in integralGraph \t";
    return nullptr;
  }
  TString name = origin->GetName();
  TH1F* newGraph = (TH1F*)origin->Clone((name+" I").Data());
  int nbins = origin->GetNbinsX();
  double error = 0.;
  double total = origin->IntegralAndError(-1,-1,error);
  for(int i = 1; i <= nbins; i++)
    newGraph->SetBinContent(i, origin->Integral(i,-1));
  
  //cout<<"Integral: "<<total<<" +- "<<error<<'\n';
  return newGraph;
}

TH1F* reverseIntegrGraph(const TH1F* const origin){
  if(origin == nullptr){
    cout<<" \tError: nullptr in reverseIntegrGraph \t";
    return nullptr;
  }
  TString name = origin->GetName();
  TH1F* newGraph = (TH1F*)origin->Clone((name+" rI").Data());
  int nbins = origin->GetNbinsX();
  for(int i = nbins; i >= 1; i--)
    newGraph->SetBinContent(i, origin->Integral(1, i));
  
  return newGraph;
}

TH1F* sqrtGraph(const TH1F* const origin){
  if(origin == nullptr){
    cout<<" \tError: nullptr in sqrtGraph \t";
    return nullptr;
  }
  TString name = origin->GetName();
  TH1F* newGraph = (TH1F*)origin->Clone((name+" Sqrt").Data());
  int nbins = origin->GetNbinsX();
  for(int i = 1; i <= nbins; i++){
    float s = sqrt(newGraph->GetBinContent(i));
    newGraph->SetBinContent(i, s);
  }
  return newGraph;
}

bool doThisGraphName(TString graphName, TString reqGrName){
  //if(SgraphName.find(NO_GRAPH) != string::npos) return false;
  graphName.ToLower();
  reqGrName.ToLower();
  std::string SgraphName = std::string(graphName.Data());
  std::string SreqGrName = std::string(reqGrName.Data());
  if(SgraphName.find(F_PATTERN) == string::npos && SgraphName.find(R_PATTERN) ==string::npos){
    //return false;
  }
  
  /*cout<<"\nSgraphName = \""<<SgraphName<<"\" \tSreqGrName = \""<<SreqGrName<<'\"';
    cout<<"\nSreqGrName != string(\"\") "<<(SreqGrName != string(""));
    cout<<"\nSgraphName.find(SreqGrName) == string::npos "<<(SgraphName.find(SreqGrName) == string::npos);
    cout<<"\nSreqGrName.find(\"all\") == string::npos "<<(SreqGrName.find("all") == string::npos);*/
  if(SreqGrName != string("")){
    if(
       SgraphName.find(SreqGrName) == string::npos && 
       SreqGrName.find("all") == string::npos
       ){
      //cout<<"\nReturning false\n";
      return false;
    }
  } //else cout<<"\nReturning true\n";
  cout<<Form("\n----- doing \"%s\" -----\n", graphName.Data());
  return true;
}

// ########################## Graphs ####################################
// ----- Service stuff  -----
template <class TH = TH1F, typename f = float>
void cutsNotScaled(const MyGraphs<TH>& theGraphs, const GNames& names){
  TCanvas* cCutSame = new TCanvas("Cuts same", "Cuts same", 10,0,1280,1024); //names.graphName + " same"
  
  float yMax = (theGraphs.vhSigs.at(0))->GetBinContent(1);
  for(size_t i = 1; i < theGraphs.vhSigs.size(); i++){
    float newY = theGraphs.vhSigs.at(i)->GetBinContent(1);
    yMax = (newY > yMax ? newY : yMax);
  }
  for(size_t i = 0; i < theGraphs.vhBkgs.size(); i++){
    float newY = theGraphs.vhBkgs.at(i)->GetBinContent(1);
    yMax = (newY > yMax ? newY : yMax);
  }
  
  TLegend* legendCut = new TLegend(0.78, 0.94, 0.98, 0.74);
  //legendCut->AddEntry(hSigC2, names.sigGrName);
  
  for(size_t i = 0; i < theGraphs.vhBkgs.size(); i++){
    TH* cg = (TH*)(theGraphs.vhBkgs.at(i)->Clone(names.bkgGrNames.at(i)+" C"));
    cg->SetTitle(names.graphName+" (not scaled) "/*+names.bkgGrNames.at(i)*/);
    cg->SetLineColor(myColors.at(theGraphs.vhSigs.size() + i));
    //cg->SetFillColor(myFillColors.at(theGraphs.vhSigs.size() + i));

    cg->GetYaxis()->SetRangeUser(1.,yMax*1.05);
    cCutSame->cd();
    cg->Draw("SAME HIST");
    legendCut->AddEntry(cg, names.bkgGrNames.at(i));
  }
  
  for(size_t i = 0; i < theGraphs.vhSigs.size(); i++){
    TH* cg = (TH*)(theGraphs.vhSigs.at(i)->Clone(names.sigGrNames.at(i)+" C"));
    cg->SetTitle("Cuts (not scaled)");
    cg->SetLineColor(myColors.at(i));
    //cg->SetFillColor(myFillColors.at(i));

    cg->GetYaxis()->SetRangeUser(1.,yMax*1.05);
    cCutSame->cd();
    cg->Draw("SAME HIST");
    legendCut->AddEntry(cg, names.sigGrNames.at(i));
  }
  
  cCutSame->cd();
  legendCut->Draw("SAME HIST");
}

template <class TH = TH1F>
void cutsScaled(const MyGraphs<TH>& theGraphs, const GNames& names, UInt_t nbins, float xm, float xM){
  TCanvas* cCutScaled = new TCanvas(names.graphName+" same norm", names.graphName+" same norm", 10,0,1280,1024);
  TLegend* legendCutNorm = new TLegend(0.78, 0.94, 0.98, 0.74);
  
  for(size_t i = 0; i < theGraphs.vhBkgs.size(); i++){
    if(theGraphs.vhBkgs.at(i) == nullptr){
      cout<<"\tError: missing \""<<names.graphName<<"\" from \""<<names.bkgGrNames.at(0)<<'\n';
      continue;//return;
    }
    TH* cgNorm = (TH*)(theGraphs.vhBkgs.at(i)->Clone(names.bkgGrNames.at(i)+" norm"));
    cgNorm->Scale(1./cgNorm->GetBinContent(1));
    //cgNorm->SetMinimum(0.);
    cgNorm->SetLineColor(myColors.at(theGraphs.vhSigs.size() + i));
    //cgNorm->SetFillColor(myFillColors.at(theGraphs.vhSigs.size() + i));
    //cgNorm->SetFillColorAlpha(myColors.at(theGraphs.vhSigs.size() + i),35);
    cCutScaled->cd();
    cgNorm->Draw("SAME HIST");
    legendCutNorm->AddEntry(cgNorm, names.bkgGrNames.at(i));
  }
  
  for(size_t i = 0; i < theGraphs.vhSigs.size(); i++){
    if(theGraphs.vhSigs.at(i) == nullptr){
      cout<<"\tError: missing \""<<names.graphName<<"\" from \""<<names.sigGrNames.at(0)<<'\n';
      continue;//return;
    }
    TH* cgNorm = (TH*)(theGraphs.vhSigs.at(i)->Clone(names.sigGrNames.at(i)+" norm"));
    cgNorm->Scale(1./cgNorm->GetBinContent(1));
    cgNorm->SetLineColor(myColors.at(i));
    cCutScaled->cd();
    cgNorm->Draw("SAME HIST");
    legendCutNorm->AddEntry(cgNorm, names.sigGrNames.at(i));
  }
  
  legendCutNorm->Draw("SAME HIST");
}

template <class TH = TH1F>
void cutSigSqrtBkg(const MyGraphs<TH>& theGraphs, const GNames& names, UInt_t nbins,float xm,float xM){
  if(theGraphs.vhSigs.at(0) == nullptr){
    cout<<"\tError: missing \""<<names.graphName<<"\" from \""<<names.sigGrNames.at(0)<< " --> can't create CutGraph (sig/sqrt(bkg))";
    return;
  }
  TCanvas* cCut = new TCanvas("Cuts sig/#sqrt{bkg}", "Cuts sig/#sqrt{bkg}", 10,0,1280,1024);
  TAxis* sig0XAxis = theGraphs.vhSigs.at(0)->GetXaxis();
  
  //Preparing Graphs
  TH* hSigSum = new TH("Cuts sig","Cuts sig", nbins, xm, xM);
  for(UInt_t i = 1; i<= nbins; i++)
    hSigSum->GetXaxis()->SetBinLabel(i, sig0XAxis->GetBinLabel(i));
  
  TH* hBkgSum = new TH("Cuts bkg","Cuts bkg", nbins, xm, xM);
  for(UInt_t i = 1; i<= nbins; i++)
    hBkgSum->GetXaxis()->SetBinLabel(i, sig0XAxis->GetBinLabel(i));
  
  //Summing graphs
  for(size_t i = 0; i < theGraphs.vhSigs.size(); i++){
    TH* cg = (TH*)(theGraphs.vhSigs.at(i)->Clone(names.sigGrNames.at(i)+" C"));
    cg->SetTitle("Cuts sqrt "+names.sigGrNames.at(i));
    cg->SetLineColor(myColors.at(i));
    hSigSum->Add(cg);
    //cout<<"\nAdded: \""<<names.sigGrNames.at(i)<<'\"';
  }
  
  for(size_t i = 0; i < theGraphs.vhBkgs.size(); i++){
    TH* cg = (TH*)(theGraphs.vhBkgs.at(i)->Clone(names.bkgGrNames.at(i)+" C"));
    cg->SetTitle("Cuts sqrt "+names.bkgGrNames.at(i));
    cg->SetLineColor(myColors.at(theGraphs.vhSigs.size() + i));
    hBkgSum->Add(cg);
    //cout<<"\nAdded: \""<<names.bkgGrNames.at(i)<<'\"';
  }
  
  cCut->cd();
  hSigSum->Divide(sqrtGraph(hBkgSum));
  hSigSum->SetTitle("Cuts sig/#sqrt{bkg};;# Events");
  hSigSum->SetMinimum(0.);
  hSigSum->Draw("HIST TEXT");
}
// ----- End service stuff---

template<class TH = TH1F>
void compareCutGraphs(const MyGraphs<TH>& theGraphs, const GNames& names){
#ifdef TEST_MODE
  cout<<"\n----compareCutGraphs: "<<names.graphName;
#endif
  if(theGraphs.vhSigs.at(0) == nullptr){
    cout<<"\tError: missing \""<<names.graphName<<"\" from \""<<names.sigGrNames.at(0)<< " --> can't create CutGraph";
    return;
  }
  //Check range
  TAxis* sig0XAxis = theGraphs.vhSigs.at(0)->GetXaxis();
  Double_t xm = sig0XAxis->GetBinLowEdge(1);
  Int_t nbins = sig0XAxis->GetNbins();
  Double_t xM = sig0XAxis->GetBinLowEdge(nbins+1);
  
  cutsNotScaled(theGraphs, names);
  cutsScaled(theGraphs, names, nbins, xm, xM);
  cutSigSqrtBkg(theGraphs, names, nbins, xm, xM);
}

// ##################################################################
template <class TH>
void printGraphSame(const MyGraphs<TH>& graphs, const GNames& names, const string& region){
#ifdef TEST_MODE
  cout<<"-----printGraphSame: "<<names.graphName;
#endif
  
  TCanvas *c0 = new TCanvas(names.graphName+" nostack", names.graphName+" nostack", 10,0,1280,1024);
  TLegend* legend0 = new TLegend(0.75,0.68,0.98,0.95);
  

  for(size_t i = 0; i < graphs.vhBkgs.size(); i++){
    if(graphs.vhBkgs.at(i) == nullptr){
      cout<<"\tError: missing \""<<names.graphName<<"\" from \""<<names.bkgGrNames.at(i)<<'\n';
      continue;
    }
    TH* cg = (TH*)(graphs.vhBkgs.at(i)->Clone(names.bkgGrNames.at(i)));
    TString newTitle = names.graphTitle + ";# Events";
    newTitle.ReplaceAll(F_PATTERN, "").ReplaceAll(R_PATTERN, "");
    cg->SetTitle(newTitle);
    cg->SetLineColor(myColors.at(i + graphs.vhSigs.size()));
    legend0->AddEntry(cg, names.bkgGrNames.at(i));
    cg->Draw("SAME HISTO");
    //cout<<"\nAdded: \""<<names.bkgGrNames.at(i)<<'\"';
  }
  
  for(size_t i = 0; i < graphs.vhSigs.size(); i++){
    if(graphs.vhSigs.at(i) == nullptr){
      cout<<"\tError: missing \""<<names.graphName<<"\" from \""<<names.sigGrNames.at(i)<<'\n';
      continue;
    }
    TH* cg = (TH*)(graphs.vhSigs.at(i)->Clone(names.sigGrNames.at(i)));
    TString newTitle = names.graphTitle + ";# Events";
    newTitle.ReplaceAll(F_PATTERN, "").ReplaceAll(R_PATTERN, "");
    cg->SetTitle(newTitle);
    cg->SetLineColor(myColors.at(i));
    legend0->AddEntry(cg, names.sigGrNames.at(i));
    cg->Draw("SAME HISTO");
    //cout<<"\nAdded: \""<<names.bkgGrNames.at(i)<<'\"';
  }
  
  legend0->Draw("SAME HISTO");
  //c0->BuildLegend(0.7,0.8,0.95,0.95);//It's a kind of magic
  TString name = names.graphName;
  name.ReplaceAll(' ', '_');
  c0->SaveAs(Form("plots/%s_%s_nostack.png", region.c_str(),name.Data()));

}

template <class TH>
void printGraphStack(const MyGraphs<TH>& graphs, const GNames& names,const string& region){
#ifdef TEST_MODE
  cout<<"\n-----printGraphStack: "<<names.graphName<<'\n';
#endif
  TCanvas *c1 = new TCanvas(names.graphName+" stack", names.graphName+" stack", 10,0,1280,1024);
  THStack* stack1 = new THStack(names.graphName+" stack", names.graphTitle+";Events");
  
  for(size_t i = 0; i < graphs.vhSigs.size(); i++){
    if(graphs.vhSigs.at(i) == nullptr){
      cout<<"\tError: missing \""<<names.graphName<<"\" from \""<<names.sigGrNames.at(i)<<'\n';
      continue;
    }
    TH* cg = (TH*)(graphs.vhSigs.at(i)->Clone(names.sigGrNames.at(i)));
    //cg->Scale(_FLIP);
    cg->SetTitle(names.sigGrNames.at(i));
    cg->SetLineColor(myColors.at(i));
    cg->SetFillColor(myFillColors.at(i));
    stack1->Add(cg);
  }
  
  for(size_t i = 0; i < graphs.vhBkgs.size(); i++){
    if(graphs.vhBkgs.at(i) == nullptr){
      cout<<"\tError: missing \""<<names.graphName<<"\" from \""<<names.bkgGrNames.at(i)<<'\n';
      continue;
    }
    TH* cg = (TH*)(graphs.vhBkgs.at(i)->Clone(names.bkgGrNames.at(i)));
    cg->Scale(_FLIP);
    std::string name_gg(names.bkgGrNames.at(i));
    if (name_gg.find("ggTo") != std::string::npos) cg->Scale(0.001);
    cg->SetTitle(names.bkgGrNames.at(i));
    cg->SetLineColor(myColors.at(i + graphs.vhSigs.size()));
    cg->SetFillColor(myFillColors.at(i + graphs.vhSigs.size()));
    stack1->Add(cg);
  }
  
  stack1->SetMaximum(stack1->GetMaximum() *1.2);
  stack1->Draw("HIST");
  
  c1->BuildLegend(0.65,0.70,0.88,0.88);//It's a kind of magic

  TString name = names.graphName;
  name.ReplaceAll(' ', '_');
  c1->SaveAs(Form("plots/%s_%s_stack.png", region.c_str(),name.Data()));
}

template <class TH>
void printGraphRatio(const MyGraphs<TH>& graphs, const GNames& names,const string& region){
#ifdef TEST_MODE
  cout<<"-----printGraphRatio: "<<names.graphName<<'\n';
#endif
  TCanvas *c1 = new TCanvas(names.graphName+" ratio", names.graphName+" ratio", 10,0,1280,1024);
  THStack* stack1 = new THStack(names.graphName+" stack", names.graphTitle+";Events");
  
  for(size_t i = 0; i < graphs.vhSigs.size(); i++){
    if(graphs.vhSigs.at(i) == nullptr){
      cout<<"\tError: missing \""<<names.graphName<<"\" from \""<<names.sigGrNames.at(i)<<'\n';
      continue;
    }
    TH* cg = (TH*)(graphs.vhSigs.at(i)->Clone(names.sigGrNames.at(i)));
    //cg->Scale(_FLIP);
    cg->SetTitle(names.sigGrNames.at(i));
    cg->SetLineColor(myColors.at(i));
    cg->SetFillColor(myFillColors.at(i));
    stack1->Add(cg);
  }
  
  for(size_t i = 0; i < graphs.vhBkgs.size(); i++){
    if(graphs.vhBkgs.at(i) == nullptr){
      cout<<"\tError: missing \""<<names.graphName<<"\" from \""<<names.bkgGrNames.at(i)<<'\n';
      continue;
    }
    TH* cg = (TH*)(graphs.vhBkgs.at(i)->Clone(names.bkgGrNames.at(i)));
    cg->Scale(_FLIP);
    std::string name_gg(names.bkgGrNames.at(i));
    if (name_gg.find("ggTo") != std::string::npos) cg->Scale(0.001);
    cg->SetTitle(names.bkgGrNames.at(i));
    cg->SetLineColor(myColors.at(i + graphs.vhSigs.size()));
    cg->SetFillColor(myFillColors.at(i + graphs.vhSigs.size()));
    stack1->Add(cg);
  }
  
  //graphs.data->Scale(_FLIP);
  if(graphs.data == nullptr){
    cout<<"\tError: missing \""<<names.graphName<<"\" from \""<<"DATA"<<"\"\n";
    return;
  }
  graphs.data->SetTitle("data");
  graphs.data->SetMarkerStyle(kFullCircle);
  graphs.data->SetMarkerColor(kBlack);
  graphs.data->SetLineColor(kBlack);
  
  stack1->SetMaximum(max(stack1->GetMaximum(), graphs.data->GetMaximum()) *1.2);
  
  //cout<<"sum name "<<sum->GetName()<<"   "<<sum->GetTitle()<<'\n';
  TH1F* sum = (TH1F*)stack1->GetStack()->Last()->Clone("sum");
  TRatioPlot* ratio = new TRatioPlot(graphs.data, sum); //"divsym"
  ratio->Draw();
  ratio->GetLowerRefGraph()->SetMarkerStyle(kFullCircle);
  ratio->GetLowerRefGraph()->SetMarkerColor(kBlack);
  ratio->GetLowerRefGraph()->SetLineColor(kBlack);
  //ratio->GetLowerRefGraph()->GetYaxis()->SetNdivisions(6, 12, 60, true);
  //ratio->GetUpperPad()->SetLogy();
  ratio->GetLowerRefGraph()->SetMaximum(2.5);
  ratio->GetLowerRefGraph()->SetMinimum(0.);
  
  ratio -> GetUpperPad()->Clear();
  ratio -> GetUpperPad()->cd();
  stack1->Draw("HIST A");
  graphs.data->Draw("SAME");
  ratio -> GetUpperPad()->BuildLegend(0.65,0.60,0.88,0.88);//It's a kind of magic
  //ratio -> GetUpperPad()->Update();

  TString name = names.graphName;
  name.ReplaceAll(' ', '_');
  c1->SaveAs(Form("plots/%s_%s_ratio.png", region.c_str(),name.Data()));
}

template <class TH>
void printGraphSqrt(const MyGraphs<TH>& graphs, const GNames& names,const string& region){
#ifdef TEST_MODE
  cout<<"\n----printGraphSqrt: "<<names.graphName;
#endif
  TString thisName = names.graphName+" I(S)/#sqrt{I(B)}";
  TCanvas *cIS = new TCanvas(thisName.Data(), thisName.Data(), 10,0,1280,1024);
  
  //Takes the first signal as a reference for the range
  if(graphs.vhSigInt.at(0) == nullptr){
    cout<<"\tError: missing \""<<names.graphName<<"\" from \""<<names.sigGrNames.at(0)<<" --> can't create the sig/sqrt(bkg) graph";
    return;
  }
  TH* hSig0 = (TH*)graphs.vhSigInt.at(0);
  //Check range
  TAxis* sig0XAxis = hSig0->GetXaxis();
  Double_t xm = sig0XAxis->GetBinLowEdge(1);
  Int_t nbins = sig0XAxis->GetNbins();
  Double_t xM = sig0XAxis->GetBinLowEdge(nbins+1);
  
  
  TH* hSigSqrt = new TH(names.graphName+" I(S)/#sqrt{I(B)}", names.graphName+" I(S)/#sqrt{I(B)}", nbins, xm, xM); // Sum of all signals. Will be divided by sqrt(sum(bkg))
  if(strcmp(hSig0->GetXaxis()->GetBinLabel(1), "") != 0)
    for(int b = 1; b <= hSig0->GetNbinsX(); ++b)
      hSigSqrt->GetXaxis()->SetBinLabel(b, hSig0->GetXaxis()->GetBinLabel(b));
  
  for(size_t i = 0; i < graphs.vhSigInt.size(); i++){
    if(graphs.vhSigInt.at(i) == nullptr){
      cout<<"\tError: missing \""<<names.graphName<<"\" from \""<<names.sigGrNames.at(i)<<'\n';
      continue;
    }
    TH* cg = (TH*)( graphs.vhSigInt.at(i)->Clone(names.sigGrNames.at(i)+" Integral") ); //Copy integral of signal
    cg->SetTitle(names.sigGrNames.at(i)+" Integral");
    //cg->SetLineColor(i);
    hSigSqrt->Add(cg);
    //cout<<"\nAdded: \""<<names.sigGrNames.at(i)<<'\"';
  }
  
  
  TH* hBkgSqrt = new TH(names.graphName+" I(allBkgs)",names.graphName+" I(allBkgs)", nbins, xm, xM);
  if(strcmp(hSig0->GetXaxis()->GetBinLabel(1), "") != 0)
    for(int b = 1; b <= hSig0->GetNbinsX(); ++b)
      hBkgSqrt->GetXaxis()->SetBinLabel(b, hSig0->GetXaxis()->GetBinLabel(b));
  
  for(size_t i = 0; i < graphs.vhBkgInt.size(); i++){
    if(graphs.vhBkgInt.at(i) == nullptr){
      cout<<"\tError: missing \""<<names.graphName<<"\" from \""<<names.bkgGrNames.at(i)<<'\n';
      continue;
    }
    TH* cg = (TH*)( graphs.vhBkgInt.at(i)->Clone(names.bkgGrNames.at(i)+" Integral") ); //Copy integral of background
    cg->SetTitle(names.bkgGrNames.at(i)+" Integral");
    //cg->SetLineColor(graphs.vhSigs.size() + i);
    hBkgSqrt->Add(cg);
    //cout<<"\nAdded: \""<<names.bkgGrNames.at(i)<<'\"';
  }
  
  hSigSqrt->Divide(sqrtGraph(hBkgSqrt)); //Division and sqrt of background combined
  
  cIS->cd();
  if(names.graphName.Contains(R_PATTERN)){
    hSigSqrt->SetLineColor(kRed);
    hSigSqrt->SetFillColor(kRed-7);
  } else {
    hSigSqrt->SetFillColor(kBlue-7);
  }
  
  TString newTitle = hSigSqrt->GetTitle();
  newTitle.ReplaceAll(F_PATTERN, "").ReplaceAll(R_PATTERN, "");
  hSigSqrt->SetTitle(newTitle);
  hSigSqrt->Draw("HIST");
  //cIS->BuildLegend(0.7,0.8,0.95,0.95);

  TString name = names.graphName;
  name.ReplaceAll(' ', '_');
  cIS->SaveAs(Form("plots/%s_%s_sqrt.png", region.c_str(),name.Data()));
}

// ########################## All the graphs in one place ####################################
template <size_t N, class TH=TH1F>
void printGraph(const MyGraphs<TH>& theGraphs, const GNames& theNames, const bitset<N>& typeToDo, const string& region){
#ifdef TEST_MODE
  cout<<"-----printGraph for \""<<theNames.graphName<<"\"\n";
#endif
  if(typeToDo.test(0))
    printGraphSame(theGraphs, theNames, region);
  
  if(typeToDo.test(1))
    printGraphStack(theGraphs, theNames, region);
  
  if(typeToDo.test(2))
    printGraphRatio(theGraphs,theNames, region);
  
  if(typeToDo.test(3))
    printGraphSqrt(theGraphs,theNames, region);
#ifdef TEST_MODE
  cout<<"-----returning from printing \""<<theNames.graphName<<"\"\n";
#endif
  return;
}

template <class TH = TH1F>
MyGraphs<TH>* buildMyGraphs(const GNames& names, vector<TFile*>& signalFiles, vector<TFile*>& backgroundFiles, TFile* dataFile){
#ifdef TEST_MODE
  cout<<"-----buildMyGraphs for \""<<names.graphName<<"\"\n";
#endif
  vector<TH*> vhSignals;
  vector<TH*> vhSigIntegr;
  vector<TH*> vhBackgrounds;
  vector<TH*> vhBackgrIntegr;
  TH* data;
  
  bool isReverse = names.graphName.Contains(R_PATTERN);
  TString strippedGrName(names.graphName);//graphName will be overridden at the end of the func
  strippedGrName.ReplaceAll(F_PATTERN, "").ReplaceAll(R_PATTERN, "");  // Replacing PATTERNS
  
  for(size_t i = 0; i < signalFiles.size(); i++){
    TH* pippo = openGraph(signalFiles.at(i), names.graphName, names.sigGrNames.at(i));
    //cout<<"\tSignal \""<<names.sigGrNames.at(i)<<"\" --> "<<(pippo!=nullptr)<<"\n";
    if(pippo == nullptr){
      cout<<"\tError: missing \""<<names.graphName<<"\" from \""<<names.sigGrNames.at(i)<<"\"\n";
      continue;//return
    }
    TH* pippoC = (TH*) ((TH1*)pippo->Clone(strippedGrName))->Rebin(1);
    vhSignals.push_back( pippoC );
    if(isReverse)
      vhSigIntegr.push_back( reverseIntegrGraph(pippoC) );
    else
      vhSigIntegr.push_back( integralGraph(pippoC) );
  }
  
  for(size_t i = 0; i < backgroundFiles.size(); i++){
    TH* pippo = openGraph(backgroundFiles.at(i), names.graphName, names.bkgGrNames.at(i));
    //cout<<"\tBackground \""<<names.bkgGrNames.at(i)<<"\" --> "<<(pippo!=nullptr)<<"\n";
    if(pippo == nullptr){
      cout<<"\tError: missing \""<<names.graphName<<"\" from \""<<names.bkgGrNames.at(i)<<"\"\n";
      continue;//return;
    }
    TH* pippoC = (TH*) ((TH1*) pippo->Clone(strippedGrName))->Rebin(1);
    vhBackgrounds.push_back( pippoC );
    if(isReverse)
      vhBackgrIntegr.push_back( reverseIntegrGraph(pippoC) );
    else
      vhBackgrIntegr.push_back( integralGraph(pippoC) );
  }
  
  data = (TH*) openGraph(dataFile, names.graphName, "data");
  if(data) data->Rebin(1);
  
#ifdef TEST_MODE
  cout<<"-----Returning from \""<<names.graphName<<"\"\n";
#endif
  //names.graphName = strippedGrName;
  return new MyGraphs<TH>(vhSignals, vhBackgrounds, vhSigIntegr, vhBackgrIntegr, data);
}
