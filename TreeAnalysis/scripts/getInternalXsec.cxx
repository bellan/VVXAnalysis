#include <iostream>
#include <utility>
#include <numeric>

#include "TTree.h"
#include "TBranch.h"
#include "TFile.h"
#include "TH1D.h"
#include "TCanvas.h"

using std::pair;
using std::make_pair;

//pair<pair<double, double>, pair<double, double>> getInternalXsec(const char* fname){
void getInternalXsec(const char* fname){
  TFile fin(fname);
  if(fin.IsZombie() || !fin.IsOpen())
    return; // make_pair(make_pair(0, 0), make_pair(0, 0));
  
  TTree* tree = (TTree*) fin.Get("treePlanter/HollyTree");
  if(! tree)
    return; // make_pair(make_pair(0, 0), make_pair(0, 0));
  
  double internalCrossSection_;
  double externalCrossSection_;
  
  TBranch *b_internalCrossSection    = 0;
  TBranch *b_externalCrossSection    = 0;
  
  tree->SetBranchAddress("internalCrossSection"   , &internalCrossSection_   , &b_internalCrossSection   );
  tree->SetBranchAddress("externalCrossSection"   , &externalCrossSection_   , &b_externalCrossSection   );

  Long64_t nentries = tree->GetEntries();
  
  // TH1D* h_internal = new TH1D("h_internal", "h_internal", 1,0,1);  // Dummy container to have GetMean() and GetStdDev()
  // TH1D* h_external = new TH1D("h_external", "h_external", 1,0,1);
  
  std::vector<double> v_internal; v_internal.reserve(nentries);
  std::vector<double> v_external; v_external.reserve(nentries);
  
  for(Long64_t e = 0; e < nentries; e++){
    //Long64_t tentry = t->LoadTree(i); tree->GetEntry(tentry);
    tree->GetEntry(e);
    v_internal.push_back(internalCrossSection_);  // h_internal->Fill(internalCrossSection_);
    v_external.push_back(externalCrossSection_);  // h_external->Fill(externalCrossSection_);
  }
  
  // std::cout<<h_internal->GetMean() <<';'<< h_internal->GetStdDev() <<';'<< h_external->GetMean() <<';'<< h_external->GetStdDev() <<'\n';
  double mean_int = 0., std_int = 0., mean_ext = 0., std_ext = 0.;
  mean_int = std::accumulate(v_internal.begin(), v_internal.end(), 0.);
  mean_ext = std::accumulate(v_external.begin(), v_external.end(), 0.);
  mean_int /= ( nentries > 0 ? nentries : 1 );
  mean_ext /= ( nentries > 0 ? nentries : 1 );
  //std_int=std::accumulate(v_internal.begin(),v_internal.end(),0.,[&mean_int](double tot,double d){return tot+(d-mean_int)*(d-mean_int);});
  //std_ext=std::accumulate(v_external.begin(),v_external.end(),0.,[&mean_ext](double tot,double d){return tot+(d-mean_ext)*(d-mean_ext);});
  for(const double& d : v_internal) std_int += (d - mean_int)*(d - mean_int);
  for(const double& d : v_external) std_ext += (d - mean_ext)*(d - mean_ext);
  std_int /= ( nentries > 1 ? nentries - 1 : 1 );
  std_ext /= ( nentries > 1 ? nentries - 1 : 1 );
  std::cout<< mean_int <<';'<< std_int <<';'<< mean_ext <<';'<< std_ext <<'\n';
    
  // auto result = make_pair( make_pair(h_internal->GetMean(), h_internal->GetStdDev()),
  // 			   make_pair(h_external->GetMean(), h_external->GetStdDev()) );
  //delete h_internal; delete h_external;
  // return result;
}
