#include "Hbos.h"
#include <DataFormats/Math/interface/deltaR.h>
#include "TFile.h"
#include <vector>
#include <algorithm>

#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "FWCore/ServiceRegistry/interface/Service.h"

using namespace std;
using namespace reco;

Hbos::Hbos(TString name_) : name(name_) {
  edm::Service<TFileService> fileService;

  hZ0Mass   = fileService->make<TH1F>(name+"_hZ0Mass", name+"_hZ0Mass", 300, 0., 140.);
  hZ1Mass   = fileService->make<TH1F>(name+"_hZ1Mass", name+"_hZ1Mass", 300, 0., 140.);
  hVMass    = fileService->make<TH1F>(name+"_hVMass", name+"_hVMass", 300, 0., 140.);

  hZPt_1    = fileService->make<TH1F>(name+"_hZPt_1", name+"_hZPt_1", 300, 0., 300.);
  hZPt_2    = fileService->make<TH1F>(name+"_hZPt_2", name+"_hZPt_2", 300, 0., 300.);
  hVPt      = fileService->make<TH1F>(name+"_hVPt", name+"_hVPt", 300, 0., 300.);
  hZZPt     = fileService->make<TH1F>(name+"_hZZPt", name+"_hZZPt", 300, 0., 300.);

  hZVDR     = fileService->make<TH1F>(name+"_hZVDR", name+"_hZVDR", 90, 0., 30. );
  hZZDR     = fileService->make<TH1F>(name+"_hZZDR", name+"_hZZDR", 90, 0., 30. );
  hZZ_VDR   = fileService->make<TH1F>(name+"_hZZ_VDR", name+"_hZZ_VDR", 90, 0., 30. );

  hZZDeta   = fileService->make<TH1F>(name+"_hZZDeta", name+"_hZZDeta", 90, 0., 30. );
  hZZ_VDeta = fileService->make<TH1F>(name+"_hZZ_VDeta", name+"_hZZ_VDeta", 90, 0., 30. );

}

Hbos::Hbos(TString name_, TFile* file) : name(name_) {

  hZ0Mass   = (TH1F*) file->Get("MyAnalyzer/"+ name+"_hZ0Mass"); 
  hZ1Mass   = (TH1F*) file->Get("MyAnalyzer/"+ name+"_hZ1Mass");
  hVMass    = (TH1F*) file->Get("MyAnalyzer/"+ name+"_hVMass");

  hZPt_1    = (TH1F*) file->Get("MyAnalyzer/"+ name+"_hZPt_1");
  hZPt_2    = (TH1F*) file->Get("MyAnalyzer/"+ name+"_hZPt_2");
  hVPt      = (TH1F*) file->Get("MyAnalyzer/"+ name+"_hVPt");
  hZZPt     = (TH1F*) file->Get("MyAnalyzer/"+ name+"_hZZPt");

  hZVDR     = (TH1F*) file->Get("MyAnalyzer/"+ name+"_hZVDR");
  hZZDR     = (TH1F*) file->Get("MyAnalyzer/"+ name+"_hZZDR");
  hZZ_VDR   = (TH1F*) file->Get("MyAnalyzer/"+ name+"_hZZ_VDR");

  hZZDeta   = (TH1F*) file->Get("MyAnalyzer/"+ name+"_hZZDeta");
  hZZ_VDeta = (TH1F*) file->Get("MyAnalyzer/"+ name+"_hZZ_VDeta");

}


void Hbos::FillBos(Boson *Z0, Boson *Z1, Boson *V){
  
  LorentzVector ZZp4 = Z0->p4() + Z1->p4();

  vector<float> pts_Z;
  for (int i=0; i<2; i++ ) {
    pts_Z.push_back((Z0->p4()).pt());
    pts_Z.push_back((Z1->p4()).pt());
  }
  sort(pts_Z.begin(),pts_Z.end());
  
  double DR1      = reco::deltaR(V->p4().eta(), V->p4().phi(), Z0->p4().eta(), Z0->p4().phi());
  double DR2      = reco::deltaR(V->p4().eta(), V->p4().phi(), Z1->p4().eta(), Z1->p4().phi());
  double DR_ZV    = std::min(DR1,DR2);
  double DR_ZZ    = reco::deltaR(Z0->p4().eta(), Z0->p4().phi(), Z1->p4().eta(), Z1->p4().phi());
  double DR_ZZ_V  = reco::deltaR(ZZp4.eta(), ZZp4.phi(), V->p4().eta(), V->p4().phi());
  float Deta_ZZ   = abs((Z0->p4().eta()) - (Z1->p4().eta()));
  float Deta_ZZ_V = abs(ZZp4.eta() - (V->p4().eta()));
  
  hZ0Mass->Fill(Z0->p4().mass());
  hZ1Mass->Fill(Z1->p4().mass());
  hVMass->Fill(V->p4().mass()); 
  
  hZPt_1->Fill(pts_Z[1]);
  hZPt_2->Fill(pts_Z[0]); 
  hVPt->Fill(V->p4().pt());
  hZZPt->Fill(ZZp4.pt()); 
 
  hZVDR->Fill(DR_ZV);
  hZZDR->Fill(DR_ZZ);
  hZZ_VDR->Fill(DR_ZZ_V);

  hZZDeta->Fill(Deta_ZZ);
  hZZ_VDeta->Fill(Deta_ZZ_V);

}

void Hbos::Scale(float w) {

  hZ0Mass->Scale(w);
  hZ1Mass->Scale(w);
  hVMass->Scale(w);
  hZPt_1->Scale(w);
  hZPt_2->Scale(w);
  hVPt->Scale(w);
  hZZPt->Scale(w);
  hZVDR->Scale(w);
  hZZDR->Scale(w);
  hZZ_VDR->Scale(w);
  hZZDeta->Scale(w);
  hZZ_VDeta->Scale(w);

}


void Hbos::SetLineColor(Color_t c) {

  hZ0Mass->SetLineColor(c);
  hZ1Mass->SetLineColor(c);
  hVMass->SetLineColor(c);
  hZPt_1->SetLineColor(c);
  hZPt_2->SetLineColor(c);
  hVPt->SetLineColor(c);
  hZZPt->SetLineColor(c);
  hZVDR->SetLineColor(c);
  hZZDR->SetLineColor(c); 
  hZZ_VDR->SetLineColor(c); 
  hZZDeta->SetLineColor(c); 
  hZZ_VDeta->SetLineColor(c);

}
