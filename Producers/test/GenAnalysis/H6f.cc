#include "H6f.h"
#include <DataFormats/HepMCCandidate/interface/GenParticle.h>
#include <FWCore/ServiceRegistry/interface/Service.h>
#include <CommonTools/UtilAlgos/interface/TFileService.h>
#include "TFile.h"
#include <vector>
#include <algorithm>


using namespace std;
using namespace reco;


H6f:: H6f(TString name_) : name(name_) {
  // Book histograms
  edm::Service<TFileService> fileService;
  h6fMass  = fileService->make<TH1F>(name+"_h6fMass", name+"_h6fMass", 300, 0., 3000.);
  h2l0Mass = fileService->make<TH1F>(name+"_h2l0Mass", name+"_h2l0Mass", 300, 0., 150.);
  h2l1Mass = fileService->make<TH1F>(name+"_h2l1Mass", name+"_h2l1Mass", 300, 0., 150.);
  hjjMass  = fileService->make<TH1F>(name+"_hjjMass", name+"_hjjMass", 300, 0., 150.);
  h4lMass  = fileService->make<TH1F>(name+"_h4lMass", name+"_h4lMass", 300, 0., 3000.);
  h2l0Pt   = fileService->make<TH1F>(name+"_h2l0Pt", name+"_h2l0Pt", 300, 0., 300.);
  h2l1Pt   = fileService->make<TH1F>(name+"_h2l1Pt", name+"_h2l1Pt", 300, 0., 300.);
  hjjPt    = fileService->make<TH1F>(name+"_hjjPt", name+"_hjjPt", 300, 0., 300.);
  h4lPt_1  = fileService->make<TH1F>(name+"_h4lPt_1", name+"_h4lPt_1", 300, 0., 300.);
  h4lPt_2  = fileService->make<TH1F>(name+"_h4lPt_2", name+"_h4lPt_2", 300, 0., 300.);
  h4lPt_3  = fileService->make<TH1F>(name+"_h4lPt_3", name+"_h4lPt_3", 300, 0., 300.);
  h4lPt_4  = fileService->make<TH1F>(name+"_h4lPt_4", name+"_h4lPt_4", 300, 0., 300.);
  hjPt_1   = fileService->make<TH1F>(name+"_hjPt_1", name+"_hjPt_1", 300, 0., 300. );
  hjPt_2   = fileService->make<TH1F>(name+"_hjPt_2", name+"_hjPt_2", 300, 0., 300. );


}



H6f::H6f(TString name_, TFile* file) : name(name_) {
  h6fMass  = (TH1F*) file->Get("GenAnalyzer/"+ name+"_h6fMass");
  h2l0Mass = (TH1F*) file->Get("GenAnalyzer/"+ name+"_h2l0Mass");
  h2l1Mass = (TH1F*) file->Get("GenAnalyzer/"+ name+"_h2l1Mass");
  hjjMass  = (TH1F*) file->Get("GenAnalyzer/"+ name+"_hjjMass");
  h4lMass  = (TH1F*) file->Get("GenAnalyzer/"+ name+"_h4lMass");
  h2l0Pt   = (TH1F*) file->Get("GenAnalyzer/"+ name+"_h2l0Pt");
  h2l1Pt   = (TH1F*) file->Get("GenAnalyzer/"+ name+"_h2l1Pt");
  hjjPt    = (TH1F*) file->Get("GenAnalyzer/"+ name+"_hjjPt");
  h4lPt_1  = (TH1F*) file->Get("GenAnalyzer/"+ name+"_h4lPt_1");
  h4lPt_2  = (TH1F*) file->Get("GenAnalyzer/"+ name+"_h4lPt_2");
  h4lPt_3  = (TH1F*) file->Get("GenAnalyzer/"+ name+"_h4lPt_3");  
  h4lPt_4  = (TH1F*) file->Get("GenAnalyzer/"+ name+"_h4lPt_4");
  hjPt_1   = (TH1F*) file->Get("GenAnalyzer/"+ name+"_hjPt_1");
  hjPt_2   = (TH1F*) file->Get("GenAnalyzer/"+ name+"_hjPt_2");
    
}


void H6f::Fill(Boson *V0, Boson *V1, Boson *V2) {

  typedef Candidate::LorentzVector LorentzVector;
  
  Candidate::LorentzVector p_2lV0 = V0->p4daughter1() + V0->p4daughter2();
  Candidate::LorentzVector p_2lV1 = V1->p4daughter1() + V1->p4daughter2();
  Candidate::LorentzVector p_jj   = V2->p4daughter1() + V2->p4daughter2();
  Candidate::LorentzVector p_4l   = p_2lV0 + p_2lV1;
  Candidate::LorentzVector p_6f   = p_2lV0 + p_2lV1 + p_jj;

 
  
  h6fMass->Fill(p_6f.mass());

  h2l0Mass->Fill(p_2lV0.mass());
  h2l1Mass->Fill(p_2lV1.mass());
  hjjMass->Fill(p_jj.mass());

  h4lMass->Fill(p_4l.mass());
    
  h2l0Pt->Fill(p_2lV0.pt());
  h2l1Pt->Fill(p_2lV1.pt());
  hjjPt->Fill(p_jj.pt());

  vector<float> pts_lep; 
  pts_lep.push_back(V0->p4daughter1().pt());
  pts_lep.push_back(V0->p4daughter2().pt());
  pts_lep.push_back(V1->p4daughter1().pt());
  pts_lep.push_back(V1->p4daughter2().pt());
  
  sort(pts_lep.begin(),pts_lep.end());
  h4lPt_1->Fill(pts_lep[3]);
  h4lPt_2->Fill(pts_lep[2]);
  h4lPt_3->Fill(pts_lep[1]);
  h4lPt_4->Fill(pts_lep[0]);
    
  vector<float> pts_j;
  pts_j.push_back(V2->p4daughter1().pt());
  pts_j.push_back(V2->p4daughter2().pt());
  
  sort(pts_j.begin(),pts_j.end());
  hjPt_1->Fill(pts_j[1]);
  hjPt_2->Fill(pts_j[0]);

}


void H6f::Scale(float w) {

  h2l0Mass->Scale(w);
  h2l1Mass->Scale(w);
  hjjMass->Scale(w);
  
  h6fMass->Scale(w);

  h4lMass->Scale(w);
  
  h2l0Pt->Scale(w);
  h2l1Pt->Scale(w);
  hjjPt->Scale(w);

  h4lPt_1->Scale(w);
  h4lPt_2->Scale(w);
  h4lPt_3->Scale(w);
  h4lPt_4->Scale(w);
  hjPt_1->Scale(w);
  hjPt_2->Scale(w);

}

void H6f::SetLineColor(Color_t c) {

  h2l0Mass->SetLineColor(c);
  h2l1Mass->SetLineColor(c);
  hjjMass->SetLineColor(c);
  
  h6fMass->SetLineColor(c);
  
  h4lMass->SetLineColor(c);  

  h2l0Pt->SetLineColor(c);
  h2l1Pt->SetLineColor(c);
  hjjPt->SetLineColor(c);
  
  h4lPt_1->SetLineColor(c);
  h4lPt_2->SetLineColor(c);
  h4lPt_3->SetLineColor(c);
  h4lPt_4->SetLineColor(c);
  hjPt_1->SetLineColor(c);
  hjPt_2->SetLineColor(c);

}
