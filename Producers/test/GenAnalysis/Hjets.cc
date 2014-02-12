#include "Hjets.h"
#include <DataFormats/HepMCCandidate/interface/GenParticle.h>
#include <FWCore/ServiceRegistry/interface/Service.h>
#include <CommonTools/UtilAlgos/interface/TFileService.h>
#include <DataFormats/Math/interface/deltaPhi.h>
#include <DataFormats/Math/interface/deltaR.h>
#include "TFile.h"
#include <vector>
#include <algorithm>

using namespace std;
using namespace reco;


Hjets:: Hjets(TString name_) : name(name_) {
  // Book histograms
  edm::Service<TFileService> fileService;

  hjjDeta = fileService->make<TH1F>(name+"_hjjDeta", name+"_hjjDeta", 100, 0., 20.);
  hjjDtheta = fileService->make<TH1F>(name+"_hjjDtheta", name+"_hjjDtheta", 100, 0., 20.);
  hjjDphi = fileService->make<TH1F>(name+"_hjjDphi", name+"_hjjDphi", 100, 0., 20.);  
  hjjDR = fileService->make<TH1F>(name+"_hjjDR", name+"_hjjDR", 100, 0., 20.);

}

Hjets::Hjets(TString name_, TFile* file) : name(name_) {

  hjjDeta = (TH1F*) file->Get("GenAnalyzer/"+ name+"_hjjDeta");
  hjjDtheta = (TH1F*) file->Get("GenAnalyzer/"+ name+"_hjjDtheta");
  hjjDphi = (TH1F*) file->Get("GenAnalyzer/"+ name+"_hjjDphi");
  hjjDR = (TH1F*) file->Get("GenAnalyzer/"+ name+"_hjjDR");

}

void Hjets::Filljet(std::vector<const reco::Candidate *> theGenq) {

  const Candidate* j0 = theGenq[0];
  const Candidate* j1 = theGenq[1];


  float Deta = abs((j0->p4().eta())-(j1->p4().eta()));
  float Dtheta = abs((j0->p4().theta())-(j1->p4().theta()));
  double Dphi = reco::deltaPhi(j0->p4().phi(), j1->p4().phi());
  double DR = reco::deltaR(j0->p4().eta(), j0->p4().phi(), j1->p4().eta(), j1->p4().phi());
  
  hjjDeta->Fill(Deta);
  hjjDtheta->Fill(Dtheta);
  hjjDphi->Fill(Dphi);
  hjjDR->Fill(DR);

}

void Hjets::Scale(float w) {

  hjjDeta->Scale(w); 
  hjjDtheta->Scale(w);

  hjjDphi->Scale(w);
  hjjDR->Scale(w);
  
}

void Hjets::SetLineColor(Color_t c) {

  hjjDeta->SetLineColor(c);
  hjjDtheta->SetLineColor(c);
  
  hjjDphi->SetLineColor(c);
  hjjDR->SetLineColor(c);
  
}
