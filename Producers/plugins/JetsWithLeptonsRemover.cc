#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"

#include "DataFormats/Math/interface/deltaR.h"

#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonSelectors.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "CommonTools/Utils/interface/StringCutObjectSelector.h"

#include "AnalysisDataFormats/CMGTools/interface/PFJet.h"
#include "DataFormats/Math/interface/Vector3D.h"
#include "VVXAnalysis/Commons/interface/Utilities.h"

#include <boost/foreach.hpp>
#define foreach BOOST_FOREACH

class  JetsWithLeptonsRemover: public edm::EDProducer {
public:
  explicit JetsWithLeptonsRemover(const edm::ParameterSet & iConfig);
  virtual ~JetsWithLeptonsRemover() { }

  virtual void produce(edm::Event & iEvent, const edm::EventSetup & iSetup);
  
  
private:
  /// Labels for input collections
  edm::InputTag jetSrc_;
  edm::InputTag muonSrc_;
  edm::InputTag electronSrc_;
  double enFractionAllowed_;

  /// Preselection cut
  StringCutObjectSelector<cmg::PFJet> preselection_;
};


JetsWithLeptonsRemover::JetsWithLeptonsRemover(const edm::ParameterSet & iConfig)
  : jetSrc_           (iConfig.getParameter<edm::InputTag>("Jets"))
  , muonSrc_          (iConfig.getParameter<edm::InputTag>("Muons"))
  , electronSrc_      (iConfig.getParameter<edm::InputTag>("Electrons"))
  , enFractionAllowed_(iConfig.getParameter<double>("EnergyFractionAllowed"))
  , preselection_     (iConfig.getParameter<std::string>("Preselection"))
{
  produces<std::vector<cmg::PFJet> >(); 
}

void JetsWithLeptonsRemover::produce(edm::Event & event, const edm::EventSetup & iSetup) {
  using namespace edm;
  using namespace std;
  
  Handle<edm::View<cmg::PFJet> > jets;
  event.getByLabel(jetSrc_, jets);


  edm::Handle<pat::MuonCollection>       muons       ; event.getByLabel(muonSrc_    ,     muons);
  edm::Handle<pat::ElectronCollection>   electrons   ; event.getByLabel(electronSrc_, electrons);

  //std::cout<<"----------- Muon -----------"<<std::endl;
  //foreach(const pat::Muon& muon, *muons)
  //  std::cout<<"pt: " << muon.pt() << " eta: " << muon.eta() << " phi: " << muon.phi() << " p: " << muon.p() <<std::endl;
  

  //std::cout<<"----------- Jets -----------"<<std::endl;
  int passPresel = 0;
  auto_ptr<vector<cmg::PFJet> > out(new vector<cmg::PFJet>());
  foreach(const cmg::PFJet& jet, *jets){
    
    if(!preselection_(jet)) continue;
    ++passPresel;

    //std::cout<<"\n+++++ Jet +++++ pt: " << jet.pt() << " eta: " << jet.eta() << " phi: " << jet.phi() << std::endl;

    bool leptonjet = false;

    const cmg::PFJetComponent mucomp = jet.component(reco::PFCandidate::ParticleType::mu);
    const cmg::PFJetComponent ecomp = jet.component(reco::PFCandidate::ParticleType::e);
    //std::cout<<"== mu comp == num: " << mucomp.number() << " fraction: " << mucomp.fraction() << std::endl; 
    //std::cout<<"== e comp == num: "  << ecomp.number()  << " fraction: " << ecomp.fraction()  << std::endl; 

    if((mucomp.number() == 0 && ecomp.number() == 0) || (mucomp.fraction() < enFractionAllowed_ && ecomp.fraction() < enFractionAllowed_)){
       out->push_back(jet);
       continue;
    }

    if(mucomp.number() > 0 && mucomp.fraction() >= enFractionAllowed_){
      math::XYZVectorD v(mucomp.pt(), 0, sqrt(mucomp.energy()*mucomp.energy() - mucomp.pt()*mucomp.pt()));
      double mucomp_abseta = v.eta();
      double mucomp_pt      = mucomp.pt();
      //std::cout<<"-- comp -- pt: " << mucomp_pt << " eta: " << mucomp_abseta                           << " p: " << mucomp.energy() << std::endl;

      foreach(const pat::Muon& muon, *muons){
	//std::cout<< (muon.pt()-mucomp_pt)/muon.pt() << " " << (abs(muon.eta())-mucomp_abseta)/abs(muon.eta()) << std::endl;
	  
	if(physmath::isAlmostEqual(muon.pt(), mucomp_pt, 0.1) && fabs(fabs(muon.eta()) - mucomp_abseta) < 0.01){
	  leptonjet = true;
	  //std::cout<<"\t\t !!! Found a muon-jet matching !!!"<<std::endl;
	  //std::cout<<"-- muon -- pt: " << muon.pt()   << " eta: " << muon.eta()    << " phi: " << muon.phi() << " p: " << muon.p()        << std::endl;
	}
      }
    }

    if(ecomp.number() > 0 && ecomp.fraction() >= enFractionAllowed_){
      math::XYZVectorD v(ecomp.pt(), 0, sqrt(ecomp.energy()*ecomp.energy() - ecomp.pt()*ecomp.pt()));
      double ecomp_abseta = v.eta();
      double ecomp_pt     = ecomp.pt();
      //std::cout<<"-- comp -- pt: "     << ecomp_pt      << " eta: " << ecomp_abseta                                    << " p: " << ecomp.energy() << std::endl;      

      foreach(const pat::Electron& electron, *electrons)
	if(physmath::isAlmostEqual(electron.pt(), ecomp_pt, 0.1) && fabs(fabs(electron.eta() - ecomp_abseta)) < 0.01){
	  leptonjet = true;
      	  //std::cout<<"\t\t !!! Found a electron-jet matching !!!"<<std::endl;
	  //std::cout<<"-- electron -- pt: " << electron.pt()   << " eta: " << electron.eta()    << " phi: " << electron.phi() << " p: " << electron.p()   << std::endl;
	}
    }

    if(!leptonjet) out->push_back(jet);
  }
  //std::cout<<"Pass Presel: "<<passPresel<<" pass cleaning: "<<out->size() << std::endl;
  event.put(out);
}

 

#include "FWCore/Framework/interface/MakerMacros.h"

DEFINE_FWK_MODULE(JetsWithLeptonsRemover);
