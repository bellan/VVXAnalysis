// -*- C++ -*-
//
// Package:    VVXAnalyzer/kfactorProducer
// Class:      kfactorProducer
// 
/**\class kfactorProducer kfactorProducer.cc 

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Gian Luca Pinna Angioni
//         Created:  Tue, 12 Apr 2016 15:57:17 GMT
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/StreamID.h"

#include <DataFormats/PatCandidates/interface/PFParticle.h>
#include "TLorentzVector.h"
#include "TSpline.h"
#include <CommonTools/UtilAlgos/interface/TFileService.h>
#include <SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h>
#include "ZZAnalysis/AnalysisStep/interface/EwkCorrections.h"
#include "ZZAnalysis/AnalysisStep/src/kFactors.C"

//
// class declaration
//

class kfactorProducer : public edm::stream::EDProducer<> {
   public:
      explicit kfactorProducer(const edm::ParameterSet&);
      ~kfactorProducer();

  //  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
  
private:
  //virtual void beginStream(edm::StreamID) override;
  virtual void produce(edm::Event&, const edm::EventSetup&) override;
  //  virtual void endStream() override;
    
  // ----------member data ---------------------------
  edm::EDGetTokenT<edm::View<reco::Candidate> > theGenVBCollectionToken;
  edm::EDGetTokenT<edm::View<reco::Candidate> > genToken_;
  edm::EDGetTokenT<GenEventInfoProduct> genInfoToken;
  
  bool isMC_ = 0;
  
  Float_t KFactorQCDqqZZ_dPhi = 1.;
  Float_t KFactorQCDqqZZ_M    = 1.;
  Float_t KFactorQCDqqZZ_Pt   = 1.;
  Float_t KFactorggZZ         = 1.;
  Float_t KFactorEWKqqZZ      = 1.;
  
  
  std::string fipPath;
  TSpline3* spkfactor;
  std::vector<std::vector<float>> ewkTable;
};


kfactorProducer::kfactorProducer(const edm::ParameterSet& config):
  isMC_ (config.getUntrackedParameter<bool>("isMC",false))
{
  genToken_ = consumes<edm::View<reco::Candidate> >(config.getParameter<edm::InputTag>("src"));
  genInfoToken = consumes<GenEventInfoProduct>(edm::InputTag("generator"));
  theGenVBCollectionToken = consumes<edm::View<reco::Candidate> >(config.getUntrackedParameter<edm::InputTag>("GenVBCollection", edm::InputTag("genCategory","vectorBosons")));

  edm::FileInPath ggzzFIP("ZZAnalysis/AnalysisStep/data/kfactors/Kfactor_Collected_ggHZZ_2l2l_NNLO_NNPDF_NarrowWidth_13TeV.root");
  fipPath=ggzzFIP.fullPath();
  TFile* ggZZKFactorFile = TFile::Open(fipPath.data());
  spkfactor = (TSpline3*)ggZZKFactorFile->Get("sp_kfactor_Nominal")->Clone();
  ggZZKFactorFile->Close();


  // Read EWK K-factor table from file
  edm::FileInPath ewkFIP("ZZAnalysis/AnalysisStep/data/kfactors/ZZ_EwkCorrections.dat");
  fipPath=ewkFIP.fullPath();
  ewkTable = EwkCorrections::readFile_and_loadEwkTable(fipPath.data());

  produces<float>("ggZZ"    );
  produces<float>("qqZZM"  );
  produces<float>("qqZZPt" );
  produces<float>("qqZZdPhi");
  produces<float>("EWKqqZZ" );
  
}


kfactorProducer::~kfactorProducer(){
 
}


//
// member functions
//

// ------------ method called to produce the data  ------------
void
kfactorProducer::produce(edm::Event& event, const edm::EventSetup& iSetup)
{
   using namespace edm;

   if(isMC_){
     edm::Handle<edm::View<reco::Candidate> > genParticles;
     event.getByToken(genToken_, genParticles);
     
     edm::Handle<edm::View<reco::Candidate> > genVBParticles;
     event.getByToken(theGenVBCollectionToken,  genVBParticles);
     
     Float_t m4l = (genVBParticles->at(0).p4()+genVBParticles->at(3).p4()).M();
     Float_t pt4l = (genVBParticles->at(0).p4()+genVBParticles->at(3).p4()).Pt();
     
     // ggZZ 
     
     KFactorggZZ = (float)spkfactor->Eval(m4l);
     
     // EWK qqZZ
 
     edm::Handle<GenEventInfoProduct> genInfo;
     event.getByToken(genInfoToken, genInfo);    
     GenEventInfoProduct genInfoP = *(genInfo.product());
     
     TLorentzVector GENZ1Vec,GENZ2Vec;
     
     if(genVBParticles->at(0).pdgId()!=23 || genVBParticles->at(3).pdgId()!= 23){
       std::cout<<"Error: The GenParticles used for kfactor have not pdg 23 23 but "<<genVBParticles->at(0).pdgId()<<" "<<genVBParticles->at(3).pdgId()<<std::endl;
       abort();
     }

     GENZ1Vec.SetPtEtaPhiM(genVBParticles->at(0).pt(),genVBParticles->at(0).eta(),genVBParticles->at(0).phi(),genVBParticles->at(0).mass());
     GENZ2Vec.SetPtEtaPhiM(genVBParticles->at(3).pt(),genVBParticles->at(3).eta(),genVBParticles->at(3).phi(),genVBParticles->at(3).mass());
     KFactorEWKqqZZ = EwkCorrections::getEwkCorrections(genParticles, ewkTable, genInfoP,GENZ1Vec,GENZ2Vec);
     
     // QCD qqZZ      

     bool sameflavor=(genVBParticles->at(0).daughter(0)->pdgId()*genVBParticles->at(0).daughter(1)->pdgId() == genVBParticles->at(3).daughter(0)->pdgId()*genVBParticles->at(3).daughter(1)->pdgId());
     
     // last argument is the order. Check it.
     KFactorQCDqqZZ_dPhi = kfactor_qqZZ_qcd_dPhi( fabs(genVBParticles->at(0).phi() - genVBParticles->at(3).phi()), (sameflavor)?1:2);  
     KFactorQCDqqZZ_M    = kfactor_qqZZ_qcd_M   ( m4l, (sameflavor)?1:2 ,2);
     KFactorQCDqqZZ_Pt   = kfactor_qqZZ_qcd_Pt  ( pt4l, (sameflavor)?1:2 );
     
   }
   
   std::auto_ptr<Float_t> result_ggZZ        (new float(KFactorggZZ));
   std::auto_ptr<Float_t> result_QCDqqZZ_dPhi(new float(KFactorQCDqqZZ_dPhi));
   std::auto_ptr<Float_t> result_QCDqqZZ_M   (new float(KFactorQCDqqZZ_M  ));
   std::auto_ptr<Float_t> result_QCDqqZZ_Pt  (new float(KFactorQCDqqZZ_Pt ));
   std::auto_ptr<Float_t> result_EWKqqZZ     (new float(KFactorEWKqqZZ));
   
   event.put(std::move(result_ggZZ        ),"ggZZ");
   event.put(std::move(result_QCDqqZZ_dPhi),"qqZZM"  );
   event.put(std::move(result_QCDqqZZ_M   ),"qqZZPt" );
   event.put(std::move(result_QCDqqZZ_Pt  ),"qqZZdPhi");
   event.put(std::move(result_EWKqqZZ     ),"EWKqqZZ" );
   
}

//define this as a plug-in
DEFINE_FWK_MODULE(kfactorProducer);
