/** \class dumpData
 *
 *  Dump all userFloat values attached to relevant collections of candidates.
 *
 *  $Date: 2013/06/06 15:40:38 $
 *  $Revision: 1.11 $
 *  \authors N. Amapane - Torino, A. Mecca - Torino
 */


#include <FWCore/Framework/interface/Frameworkfwd.h>
#include <FWCore/Framework/interface/EDAnalyzer.h>
#include <FWCore/Framework/interface/ESHandle.h>
#include <FWCore/Framework/interface/Event.h>
#include <FWCore/ParameterSet/interface/ParameterSet.h>
#include <DataFormats/Candidate/interface/Candidate.h>
#include <DataFormats/PatCandidates/interface/CompositeCandidate.h>
#include <DataFormats/PatCandidates/interface/Muon.h>
#include <DataFormats/PatCandidates/interface/Electron.h>
#include <DataFormats/PatCandidates/interface/Jet.h>
#include "DataFormats/Common/interface/ValueMap.h"
#include <ZZAnalysis/AnalysisStep/interface/PhotonFwd.h>
#include <DataFormats/Common/interface/TriggerResults.h>
#include <FWCore/Common/interface/TriggerNames.h>
#include <ZZAnalysis/AnalysisStep/interface/LeptonIsoHelper.h>

#include <iostream>
#include <iterator>
#include <string>

using namespace std;
using namespace edm;
using namespace reco;


class dumpData: public edm::EDAnalyzer {
public:
  dumpData(const ParameterSet& pset);

  void analyze(const edm::Event & event, const edm::EventSetup& eventSetup);
 
  virtual void beginJob() {};
  virtual void endJob() {};  

  void dumpCandidates(const View<pat::CompositeCandidate>& cands);
  template<typename T> void dumpUserVal(const T& cand);

  bool dumpJets;
  bool dumpPhotons;
  // bool dumpSieie, dumpChIso, dumpNeIso, dumpPhIso;
  edm::EDGetTokenT< vector<reco::Vertex> > vtxToken;
  std::map< std::string, edm::EDGetTokenT<pat::JetCollection> > jetNamedTokens;
  
  edm::EDGetTokenT< pat::PhotonCollection> photonToken;
  // edm::EDGetTokenT<edm::Handle<edm::View<pat::Photon>> > photonToken;
  // edm::EDGetTokenT< edm::ValueMap<float> > sIeIeToken_;
  // edm::EDGetTokenT< edm::ValueMap<float> > phoChIsoToken_;
  // edm::EDGetTokenT< edm::ValueMap<float> > phoNeIsoToken_;
  // edm::EDGetTokenT< edm::ValueMap<float> > phoPhIsoToken_;
  vector<string> collNames;
  vector<string> muCollNames;
  vector<string> eleCollNames;
  vector<edm::EDGetTokenT<pat::MuonCollection> > muCandidateSrcTokens;
  vector<edm::EDGetTokenT<pat::ElectronCollection> > eleCandidateSrcTokens;
  vector<edm::EDGetTokenT<edm::View<pat::CompositeCandidate> > > candidateSrcTokens;
  bool listTriggers;
  edm::EDGetTokenT<edm::TriggerResults> triggerResultToken;

};


dumpData::dumpData(const ParameterSet& pset):
  // dumpJets(   pset.existsAs<InputTag>("jetSrc"   )),
  dumpJets(pset.existsAs<ParameterSet>("jetSrc")),
  dumpPhotons(pset.existsAs<InputTag>("photonSrc")),
  // dumpSieie(  pset.existsAs<InputTag>("full5x5SigmaIEtaIEtaMap"  )),
  // dumpChIso(  pset.existsAs<InputTag>("phoChargedIsolation"      )),
  // dumpNeIso(  pset.existsAs<InputTag>("phoNeutralHadronIsolation")),
  // dumpPhIso(  pset.existsAs<InputTag>("phoPhotonIsolation"       )),
  vtxToken(consumes<vector<reco::Vertex>>(edm::InputTag("goodPrimaryVertices"))),
  // jetToken(    dumpJets    ? consumes<pat::JetCollection>(   pset.getParameter<InputTag>("jetSrc"   )) : edm::EDGetTokenT<pat::JetCollection>()    ),
  // photonToken( dumpPhotons ? consumes<edm::Handle<edm::View<pat::Photon> >>(pset.getParameter<InputTag>("photonSrc")) : edm::EDGetTokenT<edm::Handle<edm::View<pat::Photon> >>() ),
  photonToken( dumpPhotons ? consumes<pat::PhotonCollection>(pset.getParameter<InputTag>("photonSrc")) : edm::EDGetTokenT<pat::PhotonCollection>() ),
  // sIeIeToken_(   dumpSieie ? consumes<edm::ValueMap<float>>(pset.getParameter<InputTag>("full5x5SigmaIEtaIEtaMap"  )) : edm::EDGetTokenT<edm::ValueMap<float>>()),
  // phoChIsoToken_(dumpChIso ? consumes<edm::ValueMap<float>>(pset.getParameter<InputTag>("phoChargedIsolation"      )) : edm::EDGetTokenT<edm::ValueMap<float>>()),
  // phoNeIsoToken_(dumpNeIso ? consumes<edm::ValueMap<float>>(pset.getParameter<InputTag>("phoNeutralHadronIsolation")) : edm::EDGetTokenT<edm::ValueMap<float>>()),
  // phoPhIsoToken_(dumpPhIso ? consumes<edm::ValueMap<float>>(pset.getParameter<InputTag>("phoPhotonIsolation"       )) : edm::EDGetTokenT<edm::ValueMap<float>>()),
  listTriggers(pset.getUntrackedParameter<bool>("dumpTrigger",false))
{
  ParameterSet muCollps = pset.getParameter<ParameterSet>("muonSrcs");
  ParameterSet eleCollps = pset.getParameter<ParameterSet>("electronSrcs");
  ParameterSet collps = pset.getParameter<ParameterSet>("candidateSrcs");

  muCollNames = muCollps.getParameterNamesForType<InputTag>();
  for( unsigned i=0; i<muCollNames.size(); ++i) {
    muCandidateSrcTokens.push_back(consumes<pat::MuonCollection>(muCollps.getParameter<InputTag>(muCollNames[i])));
  }

  eleCollNames = eleCollps.getParameterNamesForType<InputTag>();
  for( unsigned i=0; i<eleCollNames.size(); ++i) {
    eleCandidateSrcTokens.push_back(consumes<pat::ElectronCollection>(eleCollps.getParameter<InputTag>(eleCollNames[i])));
  }

  collNames = collps.getParameterNamesForType<InputTag>();
  for( unsigned i=0; i<collNames.size(); ++i) {
    candidateSrcTokens.push_back(consumes<edm::View<pat::CompositeCandidate> >(collps.getParameter<InputTag>(collNames[i])));
  }

  triggerResultToken = consumes<edm::TriggerResults>(edm::InputTag("TriggerResults","","HLT"));
  
  ParameterSet jetCollections = pset.getParameter<ParameterSet>("jetSrc");
  std::vector<std::string> jetNames = jetCollections.getParameterNamesForType<InputTag>();
  for(const string& name : jetNames)
    jetNamedTokens[name] = consumes<pat::JetCollection>(jetCollections.getParameter<InputTag>(name));
  
  cout << "[dumpData] Succesfully constructed\n";
  cout << "\tMuon collections     (" <<  muCollNames.size() << "):";
  for(auto n:  muCollNames){ cout << n << ' ' ;} cout << '\n';
  cout << "\tElectron collections (" << eleCollNames.size() << "):";
  for(auto n: eleCollNames){ cout << n << ' ' ;} cout << '\n';
  cout << "\tCompositeCandidates collections (" << collNames.size() << "):";
  for(auto n: collNames){ cout << n << ' ' ;} cout << '\n';
  cout << "\tPhoton collection (" << dumpPhotons << "):";
  if(dumpPhotons){ cout << pset.getParameter<InputTag>("photonSrc") ;} cout << '\n';
  // cout << "\tJet collection    (" << dumpJets   << "):";
  // if(dumpJets)   { cout << pset.getParameter<InputTag>("jetSrc")    ;} cout << '\n';
  // cout << "\tPhoton value maps: sigmaIeIe=" << dumpSieie << " chIso=" << dumpChIso << " neIso=" << dumpNeIso << " phIso=" << dumpPhIso << '\n';
}


//Explicit instantiation of template function
// template void dumpData::dumpUserVal<const pat::Muon>(const pat::Muon& cand);
// template void dumpData::dumpUserVal<const pat::Electron>(const pat::Electron& cand);
// template void dumpData::dumpUserVal<const pat::CompositeCandidate>(const pat::CompositeCandidate& cand);


void dumpData::analyze(const Event & event, const EventSetup& eventSetup){

  int irun=event.id().run();
  long long int ievt=event.id().event(); 
  int ils =event.luminosityBlock();
  cout << "Dump for event " << irun << ":" << ils << ":" << ievt << endl; 

  bool dumpVertices=false;
  if (dumpVertices) {  
    edm::Handle<vector<reco::Vertex> > vtxs;
    event.getByToken(vtxToken, vtxs);

    for( vector<reco::Vertex>::const_iterator vtx =vtxs->begin(); vtx != vtxs->end(); ++vtx ) {
      float rho = vtx->position().rho();
      float z = vtx->z();
      float isFake =vtx->isFake();
      float ndof = vtx->ndof();
    
      cout << "VTX: " << rho << " " << z << " " << isFake << " " << ndof<< " " 
	   << (!isFake && ndof > 4 && abs(z) <= 24 && rho <= 2) << endl;
    }
  }
  
  

  unsigned int nColls = muCollNames.size();
  for(unsigned i=0; i<muCollNames.size(); ++i) {
    Handle<pat::MuonCollection> muons;
    int j = nColls-i-1;
    event.getByToken(muCandidateSrcTokens[j],muons);
    
    cout << muCollNames[j] << ": " << muons->size() << endl;

    for( pat::MuonCollection::const_iterator lep =muons->begin(); lep != muons->end(); ++lep ) {
      int i = distance(muons->begin(),lep);

      int genID=0;
      float genPT=0.;
      const reco::GenParticle * gp =lep->genLepton();
      if (gp) {
	genID=gp->pdgId();
	genPT=gp->pt();
      }

//   float PFChargedHadIso   = lep->pfIsolationR03().sumChargedHadronPt;
//   float PFNeutralHadIso   = lep->pfIsolationR03().sumNeutralHadronEt;
//   float PFPhotonIso       = lep->pfIsolationR03().sumPhotonEt;
//   float PFPUChargedHadIso = lep->pfIsolationR03().sumPUPt;

   float combRelIsoPF = LeptonIsoHelper::combRelIsoPF(2018, 2018, 0, *lep);

      //--- SIP, dxy, dz
      float IP      = std::abs(lep->dB(pat::Muon::PV3D));
      float IPError = lep->edB(pat::Muon::PV3D);
      float SIP     = IP/IPError;

      float dxy = 999.;
      float dz  = 999.;
      const Vertex* vertex = 0;
      edm::Handle<vector<reco::Vertex> > vtxs;
      event.getByToken(vtxToken, vtxs);
      if (vtxs->size()>0) {
         vertex = &(vtxs->front());
         dxy = fabs(lep->muonBestTrack()->dxy(vertex->position()));
         dz  = fabs(lep->muonBestTrack()->dz(vertex->position()));
      }
      cout << "#" << i << " mu"  << ((lep->charge()>0)?"+ ":"- ") << " pt= " << lep->pt() << " eta= " << lep->eta() << " phi= " << lep->phi() << " GLB= " << lep->isGlobalMuon() << " TK= " << lep->isTrackerMuon() << " matches= " << lep->numberOfMatches() << " BTT= " << lep->muonBestTrackType() << " t0_nDof: " << lep->time().nDof << " t0(ns): " << lep->time().timeAtIpInOut << " genID= " << genID <<  " genPT= " << genPT << " combRelIsoPF=" << combRelIsoPF << " SIP=" << SIP << " dxy=" << dxy << " dz=" << dz << " isPFMuon= " << lep->isPFMuon() << " muonBestTrackType= " << lep->muonBestTrackType();

//	 << " BTPT: " <<  lep->muonBestTrack()->pt() << " " << lep->innerTrack()->pt() << " " <<  lep->innerTrack()->eta() << " " << lep->innerTrack()->phi();
      // dumpUserVal(*lep);
      if (lep->hasUserData("FSRCandidates")){
	const PhotonPtrVector* fsrEle = lep->userData<PhotonPtrVector>("FSRCandidates");
	if (fsrEle->size()) {
	  cout << " Photons: pT=";	
	  for (PhotonPtrVector::const_iterator g = fsrEle->begin(); g!=fsrEle->end(); ++g) {
	    cout << " " << (*g)->pt();
	  }
	}
      }
      cout << endl;
    }
  }

  nColls = eleCollNames.size();
  for(unsigned i=0; i<eleCollNames.size(); ++i) {
    Handle<pat::ElectronCollection> electrons;
    int j = nColls-i-1;
    event.getByToken(eleCandidateSrcTokens[j],electrons);
    cout << eleCollNames[j] << ": " << electrons->size() << endl;
    for( pat::ElectronCollection::const_iterator lep = electrons->begin(); lep != electrons->end(); ++lep ) {
      int i = distance(electrons->begin(),lep);

      int genID=0;
      float genPT=0.;
      const reco::GenParticle * gp =lep->genLepton();
      if (gp) {
	genID=gp->pdgId();
	genPT=gp->pt();
      }

      cout << "#" << i << " e"  << ((lep->charge()>0)?"+  ":"-  ") << " pt= " << lep->pt() << " eta= " << lep->eta() << " phi= " << lep->phi() << " genID= " << genID <<  " genPT= " << genPT;

      // dumpUserVal(*lep);
      if (lep->hasUserData("FSRCandidates")){
	const PhotonPtrVector* fsrEle = lep->userData<PhotonPtrVector>("FSRCandidates");
	if (fsrEle->size()) {
	  cout << " Photon pTs:"; // fsrEle->size() << endl;
	  for (PhotonPtrVector::const_iterator g = fsrEle->begin(); g!=fsrEle->end(); ++g) {
	    cout << " (pt=" << (*g)->pt() ;//<< " isFromMu=" << (*g)->isFromMuon() << ")";
	  }
	}
      }
      cout << endl;
    }
  }
  
  nColls = collNames.size();
  for(unsigned i=0; i<collNames.size(); ++i) {
    Handle<View<pat::CompositeCandidate> > coll;
    int j = nColls-i-1;
    event.getByToken(candidateSrcTokens[j],coll);
    if(coll.failedToGet()) { // protection for filtered collections like TLE and RSE
      continue;
    }
    cout << collNames[j] << ": " << coll->size() << endl;
    dumpCandidates(*coll);
  }


  // if (dumpJets) {
  //   Handle<pat::JetCollection> jets;
  //   event.getByToken(jetToken, jets);

  //   cout << "Jets (only for pT>30):" << endl;
  //   for( pat::JetCollection::const_iterator jet = jets->begin(); jet != jets->end(); ++jet ) {
  //     if(jet->pt()>30){
  // 	int i = distance(jets->begin(),jet);
  // 	cout << "#" << i << " pt=" << jet->pt() << " eta=" << jet->eta() << " phi=" << jet->phi() << " combinedInclusiveSecondaryVertexV2BJetTags=" << jet->bDiscriminator("combinedInclusiveSecondaryVertexV2BJetTags");
  // 	dumpUserVal(*jet);
  // 	cout << endl;
  //     }
  //   }
  Handle<pat::JetCollection> jets;
  for(auto it = jetNamedTokens.begin() ; it != jetNamedTokens.end() ; ++it){
    event.getByToken(it->second, jets);
    
    cout << "##### Jets: " << it->first << " #####\n";
    for( pat::JetCollection::const_iterator jet = jets->begin(); jet != jets->end(); ++jet ) {
      if(jet->pt()>30){
  	int i = distance(jets->begin(),jet);
  	cout << "#" << i << " pt=" << jet->pt() << " eta=" << jet->eta() << " phi=" << jet->phi() << " combinedInclusiveSecondaryVertexV2BJetTags=" << jet->bDiscriminator("combinedInclusiveSecondaryVertexV2BJetTags") << std::endl;
  	// dumpUserVal(*jet);
	cout << "bDiscriminators: {\n";
	const std::vector<std::pair<std::string, float> >& pairDiscri = jet->getPairDiscri();
	for(const std::pair<std::string, float>& discriminator : pairDiscri)
	  cout << '\t' << discriminator.first << ": " << discriminator.second << '\n';
	cout << "}\n";
      }
    }
    cout << "########################################" << std::endl;
  }
  
  if(dumpPhotons) {
    // edm::Handle<edm::View<pat::Photon> > photons;
    Handle<pat::PhotonCollection> photons;
    event.getByToken(photonToken, photons);
    
    //   edm::Handle<edm::ValueMap<float> > sieieMap;
    //   edm::Handle<edm::ValueMap<float> > chIsoMap;
    //   edm::Handle<edm::ValueMap<float> > neIsoMap;
    //   edm::Handle<edm::ValueMap<float> > phIsoMap;
    //   if(dumpSieie) event.getByToken(sIeIeToken_, sieieMap);
    //   if(dumpChIso) event.getByToken(phoChIsoToken_, chIsoMap);
    //   if(dumpNeIso) event.getByToken(phoNeIsoToken_, neIsoMap);
    //   if(dumpPhIso) event.getByToken(phoPhIsoToken_, phIsoMap);
    
    cout<< "##### Photons #####\n";
    // for( edm::View<pat::Photon>::const_iterator photon = photons->begin(); photon != photons->end(); ++photon ){
    for( pat::PhotonCollection::const_iterator photon = photons->begin(); photon != photons->end(); ++photon ){
      int i = std::distance(photons->begin(),photon);
      // for(size_t i  = 0; i < photons->size() ; ++i){
      // 	auto photon = photons->ptrAt(i);
      cout << '#' << i << " pt=" << photon->pt() << " eta=" << photon->eta() << " phi=" << photon->phi() << '\n';
      dumpUserVal(*photon);
      // 	cout << " valueMaps: {"
      // 	     << " sieie=" << (dumpSieie ? (*sieieMap)[ photon ] : -999.)
      // 	     << " chIso=" << (dumpChIso ? (*chIsoMap)[ photon ] : -999.)
      // 	     << " neIso=" << (dumpNeIso ? (*chIsoMap)[ photon ] : -999.)
      // 	     << " phIso=" << (dumpPhIso ? (*chIsoMap)[ photon ] : -999.)
      // 	     << '}' << endl;
    }
    cout << std::endl;
  }


  // Print passing triggers
//  if (listTriggers) {
//    Handle<TriggerResults> triggerResults;
//    if (event.getByToken(triggerResultToken, triggerResults)) {
//      edm::TriggerNames const* trigNames_;  
//      trigNames_ = &event.triggerNames(*triggerResults);
//      cout << "Trigger bits:" << endl;
//      for (unsigned int i=0; i<triggerResults->size(); i++) {
//   if (triggerResults->accept(i)) cout << "   " <<
//     trigNames_->triggerName(i) << endl;
//      }
//    }
//  }
}


void dumpData::dumpCandidates(const View<pat::CompositeCandidate>& cands) {
  for( View<pat::CompositeCandidate>::const_iterator cand = cands.begin(); cand != cands.end(); ++ cand ) {
    int i = distance(cands.begin(),cand);
    cout << "#" << i << " mass: " << cand->mass() << " m0=" << cand->daughter(0)->mass() << " m1=" << cand->daughter(1)->mass()
	 << " pt0: " <<  cand->daughter(0)->pt() << " pt1: " <<  cand->daughter(1)->pt()
	 << " id1: " << cand->daughter(0)->pdgId() << " id2: " << cand->daughter(1)->pdgId();
    dumpUserVal(*cand);
    cout << endl;
  }
}

template<typename T> 
void dumpData::dumpUserVal(const T& cand) {
  const std::vector<std::string> & userLabels = cand.userFloatNames();
  //  copy(userLabels.begin(), userLabels.end(), ostream_iterator<string>(cout, " "));
  cout << "userFloats: {\n";
  for (std::vector<std::string>::const_iterator name = userLabels.begin(); name!= userLabels.end(); ++name){
    cout << '\t' << *name << "=" << cand.userFloat(*name) << '\n';
  }
  cout << "}\n";

  cout << "userInts: {\n";
  const std::vector<std::string> & userILabels = cand.userIntNames();
  //  copy(userLabels.begin(), userLabels.end(), ostream_iterator<string>(cout, " "));
  for (std::vector<std::string>::const_iterator name = userILabels.begin(); name!= userILabels.end(); ++name){
    cout << '\t' << *name << "=" << cand.userInt(*name) << '\n';
  }
  cout << "}\n";

}


#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(dumpData);

