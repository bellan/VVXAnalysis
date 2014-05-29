#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/ServiceRegistry/interface/Service.h"

#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/Math/interface/deltaPhi.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonSelectors.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/CompositeCandidate.h"
#include "CommonTools/Utils/interface/StringCutObjectSelector.h"

#include "AnalysisDataFormats/CMGTools/interface/PFJet.h"
#include "DataFormats/Math/interface/Vector3D.h"
#include "VVXAnalysis/Commons/interface/Utilities.h"
#include "VVXAnalysis/Commons/interface/Colours.h"

#include <boost/foreach.hpp>
#define foreach BOOST_FOREACH

using namespace colour;

class  JetsWithLeptonsRemover: public edm::EDProducer {
public:
  explicit JetsWithLeptonsRemover(const edm::ParameterSet & iConfig);
  virtual ~JetsWithLeptonsRemover() { }

  virtual void produce(edm::Event & iEvent, const edm::EventSetup & iSetup);
  bool isGood(const cmg::PFJet jet, int year) const;
  
private:
  /// Labels for input collections
  edm::InputTag jetSrc_;
  edm::InputTag muonSrc_;
  edm::InputTag electronSrc_;
  edm::InputTag diBosonSrc_;
  bool tagOnly_;
  double enFractionAllowed_;
  bool cleaningFromDiboson_;
 
  /// Preselection cut
  StringCutObjectSelector<cmg::PFJet> preselectionJ_;
  StringCutObjectSelector<pat::CompositeCandidate> preselectionVV_;

  // Some istograms for monitoring
  bool  activateDebugPrintOuts_;
  bool  doDebugPlots_;
  TH1F *hNLeptonJets;
  TH1F *hDeltaPt_jet_lepton;
  TH1F *hDeltaPt_jetcomp_lepton;
  TH1F *hDeltaPt_jet_fsr;
  TH1F *hDeltaPt_jetcomp_fsr;
  TH1F *hDeltaPhi_jet_lepton;
  TH1F *hDeltaPhi_jet_fsr;
  TH1F *hDeltaEta_jet_lepton;
  TH1F *hDeltaEta_jetcomp_lepton;
  TH1F *hDeltaEta_jet_fsr;
  TH1F *hDeltaEta_jetcomp_fsr;
};


JetsWithLeptonsRemover::JetsWithLeptonsRemover(const edm::ParameterSet & iConfig)
  : jetSrc_           (iConfig.getParameter<edm::InputTag>("Jets"))
  , muonSrc_          (iConfig.getParameter<edm::InputTag>("Muons"))
  , electronSrc_      (iConfig.getParameter<edm::InputTag>("Electrons"))
  , diBosonSrc_       (iConfig.getParameter<edm::InputTag>("Diboson"))  
  , tagOnly_          (iConfig.getParameter<bool>("TagOnly")) // Tag muons and electrons - matching jets, BUT still discard the one that matches diboson grand-daughters
  , enFractionAllowed_(iConfig.getParameter<double>("EnergyFractionAllowed"))
  , preselectionJ_    (iConfig.getParameter<std::string>("JetPreselection"))
  , preselectionVV_   (iConfig.getParameter<std::string>("DiBosonPreselection"))

  , activateDebugPrintOuts_ (iConfig.getUntrackedParameter<bool>("DebugPrintOuts",false))   
  , doDebugPlots_           (iConfig.getUntrackedParameter<bool>("DebugPlots",false)) 
{
  produces<std::vector<cmg::PFJet> >(); 

  if(diBosonSrc_.label() != "") cleaningFromDiboson_ = true;
  else cleaningFromDiboson_ = false;

  if(doDebugPlots_){
    edm::Service<TFileService> fs;
    hNLeptonJets              = fs->make<TH1F>("hNLeptonJets"            , "Number of lepton-jets found",  10,   0, 10);
    hDeltaPt_jet_lepton       = fs->make<TH1F>("hDeltaPt_jet_lepton"     , "#Delta p_T (jet, l)"        , 100, -50, 50);
    hDeltaPt_jetcomp_lepton   = fs->make<TH1F>("hDeltaPt_jetcomp_lepton" , "#Delta p_T (jetcomp, l)"    , 100, -50, 50);
    hDeltaPt_jet_fsr          = fs->make<TH1F>("hDeltaPt_jet_fsr"        , "#Delta p_T (jet, fsr)"      , 100, -50, 50);
    hDeltaPt_jetcomp_fsr      = fs->make<TH1F>("hDeltaPt_jetcomp_fsr"    , "#Delta p_T (jetcomp, fsr)"  , 100, -50, 50);
    hDeltaPhi_jet_lepton      = fs->make<TH1F>("hDeltaPhi_jet_lepton"    , "#Delta #phi (jet, l)"       , 100,  -4,  4);
    hDeltaPhi_jet_fsr	      = fs->make<TH1F>("hDeltaPhi_jet_fsr"       , "#Delta #phi (jet, fsr)"     , 100,  -4,  4);
    hDeltaEta_jet_lepton      = fs->make<TH1F>("hDeltaEta_jet_lepton"    , "#Delta #eta (jet, l)"       , 100,  -5,  5);
    hDeltaEta_jetcomp_lepton  = fs->make<TH1F>("hDeltaEta_jetcomp_lepton", "#Delta #eta (jetcomp, l)"   , 100,  -5,  5);
    hDeltaEta_jet_fsr	      = fs->make<TH1F>("hDeltaEta_jet_fsr"       , "#Delta #eta (jet, fsr)"     , 100,  -5,  5);
    hDeltaEta_jetcomp_fsr     = fs->make<TH1F>("hDeltaEta_jetcomp_fsr"   , "#Delta #eta (jetcomp, fsr)" , 100,  -5,  5);
  }
}

void JetsWithLeptonsRemover::produce(edm::Event & event, const edm::EventSetup & iSetup) {
  using namespace edm;
  using namespace std;
  
  Handle<edm::View<cmg::PFJet> > jets;
  event.getByLabel(jetSrc_, jets);


  //std::cout<<"----------- Muon -----------"<<std::endl;
  //foreach(const pat::Muon& muon, *muons)
  //  std::cout<<"pt: " << muon.pt() << " eta: " << muon.eta() << " phi: " << muon.phi() << " p: " << muon.p() <<std::endl;
  

  if(activateDebugPrintOuts_) std::cout << "\n\n----------- NEW EVENT ----------- number of jets: " << jets->size() << std::endl;
  int passPresel = 0;
  int numLepJets = 0;
  auto_ptr<vector<cmg::PFJet> > out(new vector<cmg::PFJet>());
  foreach(const cmg::PFJet& jet, *jets){
    
    if(!preselectionJ_(jet) && isGood(jet,2012)) continue;
    ++passPresel;



    if(activateDebugPrintOuts_) std::cout<<"\n+++++ Jet +++++ pt: " << jet.pt() << " eta: " << jet.eta() << " phi: " << jet.phi() << std::endl;

    bool leptonjet = false;

    const cmg::PFJetComponent mucomp     = jet.component(reco::PFCandidate::ParticleType::mu);
    const cmg::PFJetComponent ecomp      = jet.component(reco::PFCandidate::ParticleType::e);
   

    //if((mucomp.number() == 0 && ecomp.number() == 0) || (mucomp.fraction() < enFractionAllowed_ && ecomp.fraction() < enFractionAllowed_)){
    //   out->push_back(jet);
    //   continue;
    // }

    math::XYZVectorD vm(mucomp.pt(), 0, sqrt(mucomp.energy()*mucomp.energy() - mucomp.pt()*mucomp.pt()));
    double mucomp_abseta = vm.eta();
    double mucomp_pt      = mucomp.pt();
    
    math::XYZVectorD ve(ecomp.pt(), 0, sqrt(ecomp.energy()*ecomp.energy() - ecomp.pt()*ecomp.pt()));
    double ecomp_abseta = ve.eta();
    double ecomp_pt     = ecomp.pt();
    
    if(activateDebugPrintOuts_){
      std::cout << "== mu comp == num: " << mucomp.number() << " fraction: " << mucomp.fraction() << " pt: " << mucomp_pt << " eta: " << mucomp_abseta << std::endl; 
      std::cout << "== e comp == num: "  << ecomp.number()  << " fraction: " << ecomp.fraction()  << " pt: " << ecomp_pt << " eta: " << ecomp_abseta << std::endl; 
    }

    if(cleaningFromDiboson_){
      edm::Handle<edm::View<pat::CompositeCandidate> > VV   ; event.getByLabel(diBosonSrc_, VV);
      pat::CompositeCandidate bestVV; bool found = false;
      foreach(const pat::CompositeCandidate &vv, *VV)
	if (preselectionVV_(vv)){
	  bestVV = vv; 
	  found = true;
	  break;
	}
      
      if(found) { 
	for(int i=0; i<2; ++i){
	  
	  const pat::CompositeCandidate *v =  dynamic_cast<const pat::CompositeCandidate*>(bestVV.daughter(i)->masterClone().get());
	  

	  for(int j=0; j<2; ++j){
	      double lepcomp_pt     = -999999999;
	      double lepcomp_abseta = -999999999;
	      if(abs(v->daughter(j)->pdgId()) == 13){
		lepcomp_pt     = mucomp_pt;
		lepcomp_abseta = mucomp_abseta;
	      }
	      else if(abs(v->daughter(j)->pdgId()) == 11){
		lepcomp_pt     = ecomp_pt;
		lepcomp_abseta = ecomp_abseta;
	      }
	      else std::cout << "Do not know what to do ... do you know what are you doing?" << std::endl;

	      if(activateDebugPrintOuts_)
		std::cout << "ID lepton: " << v->daughter(j)->pdgId() << " pt: " << v->daughter(j)->pt()   << " eta: " << v->daughter(j)->eta() << " phi: " << v->daughter(j)->phi() << " p: " << v->daughter(j)->p() << std::endl;
	      if(physmath::isAlmostEqual(v->daughter(j)->pt(), lepcomp_pt, 0.1) && fabs(fabs(v->daughter(j)->eta()) - lepcomp_abseta) < 0.01){
		leptonjet = true;
		if(activateDebugPrintOuts_) std::cout << Green("\t\t !!! Found a matching lepton-jet !!!")<<std::endl;
		if(doDebugPlots_){
		  hDeltaPt_jet_lepton     ->Fill(v->daughter(j)->pt()  - jet.pt());
		  hDeltaPt_jetcomp_lepton ->Fill(v->daughter(j)->pt()  - lepcomp_pt);    
		  hDeltaPhi_jet_lepton    ->Fill(v->daughter(j)->phi() - jet.phi());
		  hDeltaEta_jet_lepton    ->Fill(v->daughter(j)->eta() - jet.eta());  
		  hDeltaEta_jetcomp_lepton->Fill(fabs(v->daughter(j)->eta()) - lepcomp_abseta);    
		}
	      }
	  }
	  
	  // Check if the jet matches FSR photons
	  if(v->hasUserFloat("dauWithFSR") && v->userFloat("dauWithFSR") >= 0){
	    
	    const cmg::PFJetComponent photoncomp = jet.component(reco::PFCandidate::ParticleType::gamma);
	    
	    double photon_en_frac = v->daughter(2)->energy()/jet.energy();
	    if(activateDebugPrintOuts_)
	      std::cout << "Sister of " << v->userFloat("dauWithFSR") << " (" <<  v->daughter(v->userFloat("dauWithFSR"))->pdgId() << "), pt: "
			<< v->daughter(2)->pt()   << " eta: " << v->daughter(2)->eta() << " phi: " << v->daughter(2)->phi() << " p: " << v->daughter(2)->p()
			<< " Photon energy fraction in the jet: " <<  photon_en_frac 
			<< std::endl;
	    if(photoncomp.number() > 0 && photon_en_frac > 0.5 && reco::deltaR(*v->daughter(2), jet) < 0.05){
	      leptonjet = true; 
	      if(activateDebugPrintOuts_) std::cout << Blue("\t\t !!! Found a matching FSR lepton-jet !!!")<<std::endl;	  
	      if(doDebugPlots_){
		math::XYZVectorD vp(photoncomp.pt(), 0, sqrt(photoncomp.energy()*photoncomp.energy() - photoncomp.pt()*photoncomp.pt()));
		hDeltaPt_jet_fsr     ->Fill(v->daughter(2)->pt()  - jet.pt());
		hDeltaPt_jetcomp_fsr ->Fill(v->daughter(2)->pt()  - photoncomp.pt());    
		hDeltaPhi_jet_fsr    ->Fill(v->daughter(2)->phi() - jet.phi());
		hDeltaEta_jet_fsr    ->Fill(v->daughter(2)->eta() - jet.eta());  
		hDeltaEta_jetcomp_fsr->Fill(fabs(v->daughter(2)->eta()) - vp.eta()); 
	      }
	    }
	  }
	}
      }
      if(leptonjet) ++numLepJets;
    }

    


    // ----------------------------------------------------------------------------------------------------------------------------------------------------------------
    // If it is a jet matching any of the leptons from the di-boson, or from a FSR photon, then it will not be loaded inside the jet cleaned container (no matter what)
    if(leptonjet) continue;
    // ----------------------------------------------------------------------------------------------------------------------------------------------------------------




    // Now, check for jets originated from extra leptons in the event. Here, the jets can be tagged only, for a later study of them.
 
    // Check for muon-originated jets
    if(!leptonjet && mucomp.number() > 0 && mucomp.fraction() >= enFractionAllowed_){
      edm::Handle<pat::MuonCollection>       muons       ; event.getByLabel(muonSrc_    ,     muons);
      
      //std::cout<<"-- comp -- pt: " << mucomp_pt << " eta: " << mucomp_abseta                           << " p: " << mucomp.energy() << std::endl;

      foreach(const pat::Muon& muon, *muons){
	// Perform the matching only if the lepton pass certain requirements

	//std::cout<< (muon.pt()-mucomp_pt)/muon.pt() << " " << (abs(muon.eta())-mucomp_abseta)/abs(muon.eta()) << std::endl;
	  
	if(physmath::isAlmostEqual(muon.pt(), mucomp_pt, 0.1) && fabs(fabs(muon.eta()) - mucomp_abseta) < 0.01){
	  leptonjet = true;
	  //std::cout<<"\t\t !!! Found a muon-jet matching !!!"<<std::endl;
	  //std::cout<<"-- muon -- pt: " << muon.pt()   << " eta: " << muon.eta()    << " phi: " << muon.phi() << " p: " << muon.p()        << std::endl;  
	}
      }
    }

    // Check for electron-originated jets
    if(!leptonjet && ecomp.number() > 0 && ecomp.fraction() >= enFractionAllowed_){
      edm::Handle<pat::ElectronCollection>   electrons   ; event.getByLabel(electronSrc_, electrons); 

      //std::cout<<"-- comp -- pt: "     << ecomp_pt      << " eta: " << ecomp_abseta                                    << " p: " << ecomp.energy() << std::endl;      

      foreach(const pat::Electron& electron, *electrons){

	// Perform the matching only if the lepton pass certain requirements

	if(physmath::isAlmostEqual(electron.pt(), ecomp_pt, 0.1) && fabs(fabs(electron.eta() - ecomp_abseta)) < 0.01){
	  leptonjet = true;
      	  //std::cout<<"\t\t !!! Found a electron-jet matching !!!"<<std::endl;
	  //std::cout<<"-- electron -- pt: " << electron.pt()   << " eta: " << electron.eta()    << " phi: " << electron.phi() << " p: " << electron.p()   << std::endl;
	}
      }
    }

    
    // FIXME asso map!

    if(tagOnly_ || !leptonjet) out->push_back(jet);
  }

  if(doDebugPlots_) hNLeptonJets->Fill(numLepJets);
  
  if(activateDebugPrintOuts_){
    std::cout<<"Pass Presel: "<<passPresel<<" pass cleaning: "<<out->size() << std::endl;
    std::cout<<"-------------------------------------------------------------------------"<< std::endl;
  }

  event.put(out);
}



bool JetsWithLeptonsRemover::isGood(const cmg::PFJet jet, int year) const {
  
  float jeta=fabs(jet.eta());
  
  //		bool looseJetID = (jet.getSelection("cuts_looseJetId") > 0) ;
  bool looseJetID = (jet.component(5).fraction() < 0.99 && 
		     jet.component(4).fraction() < 0.99 && 
		     jet.nConstituents() > 1 && 
		     ( jet.component(1).fraction() > 0 || jeta > 2.4 )  &&
		     ( ( jet.component(1).number() + jet.component(2).number() + jet.component(3).number() ) > 0 || jeta > 2.4 ) &&
		     ( jet.component(2).fraction() < 0.99 || jeta > 2.4 ) );
  
  if(!looseJetID) return false;

  bool passPU = true; //jet.passPuJetId("full", PileupJetIdentifier::kLoose);
  
  //HARD CODED implementation of JetMET V00-03-04 WPs - for synch only
  //#4 Eta Categories  0-2.5 2.5-2.75 2.75-3.0 3.0-5.0
  //Pt010_Loose    = cms.vdouble(-0.95,-0.96,-0.94,-0.95),
  //Pt1020_Loose   = cms.vdouble(-0.95,-0.96,-0.94,-0.95),
  //Pt2030_Loose   = cms.vdouble(-0.63,-0.60,-0.55,-0.45),
  //Pt3050_Loose   = cms.vdouble(-0.63,-0.60,-0.55,-0.45), 
  float jpt=jet.pt();
  float jpumva=0.;
  if(year==2012)jpumva=jet.puMva("full53x");
  else jpumva=jet.puMva("full");
  if(jpt>20){
    if(jeta>3.){
      if(jpumva<=-0.45)passPU=false;
    }else if(jeta>2.75){
      if(jpumva<=-0.55)passPU=false;
    }else if(jeta>2.5){
      if(jpumva<=-0.6)passPU=false;
    }else if(jpumva<=-0.63)passPU=false;
  }else{
    if(jeta>3.){
      if(jpumva<=-0.95)passPU=false;
    }else if(jeta>2.75){
      if(jpumva<=-0.94)passPU=false;
    }else if(jeta>2.5){
      if(jpumva<=-0.96)passPU=false;
    }else if(jpumva<=-0.95)passPU=false;
  }

  return looseJetID && passPU;
}

#include "FWCore/Framework/interface/MakerMacros.h"

DEFINE_FWK_MODULE(JetsWithLeptonsRemover);
