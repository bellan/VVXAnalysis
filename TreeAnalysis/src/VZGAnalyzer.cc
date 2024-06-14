
#include "VVXAnalysis/TreeAnalysis/interface/VZGAnalyzer.h"
#include "VVXAnalysis/Commons/interface/SignalDefinitions.h"
#include "VVXAnalysis/Commons/interface/Utilities.h"
#include "VVXAnalysis/Commons/interface/Comparators.h"

#include "VVXAnalysis/Commons/interface/GenVBHelper.h"


#include <boost/foreach.hpp>
#define foreach BOOST_FOREACH

// #include <TString.h>

#include <boost/assign/std/vector.hpp>
using namespace boost::assign;

using std::cout;
using std::endl;

using namespace phys;

bool verbose = false;
bool signalSample = true;
double etacut=4.7;
double ptcut=30;
double LumiSF=1.;//7.035;


double cosOmega(TLorentzVector a, TLorentzVector b){
  return a.CosTheta()*b.CosTheta()+TMath::Sin(a.Theta())*TMath::Sin(b.Theta())*TMath::Cos(a.Phi()-b.Phi());
}

// ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- 
double Legendre(int l, double x){
  double P=1;
  switch(l)
    {
    case 0:
      P=1;
      break;
      
    case 1:
      P=x;
      break;
  
    case 2:
      P=0.5*(3*x*x-1);
      break;

    case 3:
      P=0.5*x*(5*x*x-3);      
      break;
  
    case 4:
      P=(1./8.)*(35*x*x*x*x-30*x*x+3);
      break;
      
    case 5:
      P=(1./8.)*x*(63*x*x*x*x-70*x*x+15);
      break;
    
    case 6:
      P=(1./16.)*(231*pow(x,6)-315*pow(x,4)+105*x*x-5);
      break;
    
    default:
      P=1;
      break;
    }
  return P;
}

// ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- 
double WeightsFWM(char wi, TLorentzVector a, TLorentzVector b){
  //double W=1;
  double den, num;
  double ya=a.Rapidity();
  double yb=b.Rapidity();
  double ymean=(ya+yb)/2.;

  switch(wi)
    {
    case 's':
      num = a.Mag()*b.Mag();
      den = (a+b)*(a+b);
      break;

    case 'p':
      num = a.Mag()*b.Mag();
      den = (a+b).Mag2();
      break;
  
    case 't':
      num = a.Pt()*b.Pt();
      den = (a.Pt()+b.Pt())*(a.Pt()+b.Pt());
      break;

    case 'z':
      num = a.Pz()*b.Pz();
      den = ( (a+b).Pz() )*( (a+b).Pz() );
      break;
  
    case 'y':
      ya=1./fabs(ya-ymean);
      yb=1./fabs(yb-ymean);
      num =ya*yb;
      den =(ya+yb)*(ya+yb);
      break;
    
    default:
      num=1.;
      den=1.;
      break;
    }
  return num/den;
}

/*
// ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- 
double FWM(int l, char wi, TLorentzVector* objects){
  double H=0; 
  int N =objects->size();
  for(int i=0; i<N; i++)
    for(int j=0; j<N; j++)
      H+=WeightsFWM(wi,objects[i],objects[j])*Legendre(l, cosOmega(objects[i],objects[j]));
  //return H;
  return H/2.;
  
}
*/

// ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- 
double FWM(int l, char wi, TLorentzVector a, TLorentzVector b){
  return WeightsFWM(wi,a,b)*Legendre( l, cosOmega(a,b) );
}

// ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- 
double SumFWM(int l, char wi, std::vector <TLorentzVector> objects){//, int N){
  double H=0; 
  int N =objects.size();
  for(int i=0; i<N-1; i++)
    for(int j=i+1; j<N; j++)
      H+=FWM(l,wi,objects[i],objects[j]);

  return H;
}


bool KinematicsOK(phys::Particle p, float pt,float eta)
{
  if (fabs(p.eta()) < eta && fabs(p.pt()) > pt) return true;
  return false;
}


bool VZGAnalyzer::LeptonicSignalConstraint()
{
  //int LeptonicZdecays=genVBHelper_.ZtoChLep().size();//+genVBHelper_.ZtoNeutrinos().size();
  if (genVBHelper_.ZtoChLep().size()==1
      && KinematicsOK(genVBHelper_.ZtoChLep()[0].daughter(0), 5.,2.5)
      && KinematicsOK(genVBHelper_.ZtoChLep()[0].daughter(1), 5.,2.5)
      && genVBHelper_.ZtoChLep()[0].mass()>60 && genVBHelper_.ZtoChLep()[0].mass()<120)
    return true;
  return false;
}
bool VZGAnalyzer::HadronicSignalConstraint()
{
  if ((genVBHelper_.ZtoQ().size()==1
       && KinematicsOK(genVBHelper_.ZtoQ()[0].daughter(0), ptcut,etacut)
       && KinematicsOK(genVBHelper_.ZtoQ()[0].daughter(1), ptcut,etacut)
       && genVBHelper_.ZtoQ()[0].mass()>50 && genVBHelper_.ZtoQ()[0].mass()<120)
      ||(genVBHelper_.WtoQ().size()==1
	 && KinematicsOK(genVBHelper_.WtoQ()[0].daughter(0), ptcut,etacut)
	 && KinematicsOK(genVBHelper_.WtoQ()[0].daughter(1), ptcut,etacut)
	 && genVBHelper_.WtoQ()[0].mass()>50
	 && genVBHelper_.WtoQ()[0].mass()<120))
    return true;
  return false;
}

bool VZGAnalyzer::PhotonSignalConstraint()
{
       std::vector<phys::Particle> selectedphotons;
       for (auto p : *genParticles)
              if (p.id() == 22 && KinematicsOK(p, 20, 2.4) && p.genStatusFlags().test(phys::isPrompt) &&  p.genStatusFlags().test(phys::fromHardProcess))
                  selectedphotons.push_back(p);
       if(verbose==true) std::cout<< "Number of selected gen photons = "<<selectedphotons.size()<<std::endl;
       if (selectedphotons.size()>=1)
	 return true;
       return false;
}

bool VZGAnalyzer::IN_GENsignalDef()
{
  if (LeptonicSignalConstraint() && HadronicSignalConstraint() && PhotonSignalConstraint())
    return true;
  return false;
}


bool VZGAnalyzer::baselineRequirements()
{
  //----------------------------------------Building jj pairs ----------------------------------------//

  std::vector<phys::Jet> selectedRECOjets;
  std::vector<phys::Boson<phys::Jet>> DiJets;

  foreach (const phys::Jet &jet, *jets)
    if (KinematicsOK(jet, ptcut, etacut)) // KinematicsOK(jet)
      selectedRECOjets.push_back(jet);

  for (size_t i = 0; i < selectedRECOjets.size(); i++)
    {
      phys::Jet jetA = selectedRECOjets[i];

      for (size_t j = i + 1; j < selectedRECOjets.size(); j++)
	{
	  phys::Jet jetB = selectedRECOjets[j];
	  float mjj = (jetA.p4() + jetB.p4()).M();

	  if (verbose == true) std::cout << "mjj= " << mjj << std::endl;

	  if (mjj > 50 && mjj < 120) 
	    DiJets.push_back(phys::Boson<phys::Jet>(jetA, jetB));
                     
	}
    }

  if(verbose == true) std::cout << "DiJets size: " << DiJets.size() << std::endl;
  //-------------------------------------Requirements on Photons----------------------------------------//

  std::vector<phys::Photon> selectedphotons;
  std::vector<phys::Photon> selectedKinPhotons;
  std::vector<phys::Photon> selectedVLPhotons;
  std::vector<phys::Photon> selectedLoosePhotons;

  for (auto p : *photons)
    {
      if (p.id() == 22 && KinematicsOK(p, 20, 2.4) && !p.hasPixelSeed() && p.passElectronVeto())
	{
	  selectedKinPhotons.push_back(p);
	  if(p.cutBasedID(Photon::IdWp::VeryLoose))
	    {
	      selectedVLPhotons.push_back(p);
	      if (p.cutBasedIDLoose())
		{
		  selectedLoosePhotons.push_back(p);
		  selectedphotons.push_back(p);
		}
	    }
	}
    }
  theHistograms->fill("Atleast1KinPhoton", "Atleast1KinPhoton", 2, 0, 2, selectedKinPhotons.size() >0 , theWeight*LumiSF);
  theHistograms->fill("Atleast1VLPhoton", "AtleastVLPhoton", 2, 0, 2, selectedVLPhotons.size() >0 , theWeight*LumiSF);
  theHistograms->fill("Atleast1LoosePhoton", "Atleast1LoosePhoton", 2, 0, 2, selectedLoosePhotons.size() >0 , theWeight*LumiSF);

    
  if(verbose == true) std::cout<< "Number of selected RECO photons = "<<selectedphotons.size()<<std::endl;

  bool goodZ= (Z->mass() > 60 && Z->mass() < 120 && KinematicsOK(Z->daughter(0), 5, 2.5) && KinematicsOK(Z->daughter(1), 5, 2.5)); // KinematicsOK(jet)

       
  if (DiJets.size() > 0 && selectedphotons.size() >0 && goodZ )
    return true;
  return false;
}


Bool_t VZGAnalyzer::cut(Int_t n, phys::Boson<phys::Jet> recoV, std::vector<phys::Photon> selectedPhotons)
{ // returns false if the event has to be cut

    phys::Photon mostEnergeticPhoton;
    if (selectedPhotons.size() > 0)
      {
	std::stable_sort(selectedPhotons.begin(), selectedPhotons.end(), phys::EComparator());
	mostEnergeticPhoton = selectedPhotons[0];
      }


    std::pair<phys::Photon, phys::Jet> nearestRECOjetstoPhoton;

    if(abs(physmath::deltaR(recoV.daughter(0), mostEnergeticPhoton))<abs(physmath::deltaR(recoV.daughter(1), mostEnergeticPhoton)) )
      nearestRECOjetstoPhoton={mostEnergeticPhoton, recoV.daughter(0)};
    else if (abs(physmath::deltaR(recoV.daughter(0), mostEnergeticPhoton))>abs(physmath::deltaR(recoV.daughter(1), mostEnergeticPhoton)) )
      nearestRECOjetstoPhoton={mostEnergeticPhoton, recoV.daughter(1)};
  
  std::vector<TLorentzVector> jjG;
  
  std::vector<TLorentzVector> lljjG;

  if(selectedPhotons.size()>0)
    {
      jjG.push_back(recoV.daughter(0).p4());
      jjG.push_back(recoV.daughter(1).p4());
      jjG.push_back(selectedPhotons.at(0).p4());

      lljjG.push_back(Z->daughter(0).p4());
      lljjG.push_back(Z->daughter(1).p4());
      lljjG.push_back(jjG.at(0));
      lljjG.push_back(jjG.at(1));
      lljjG.push_back(jjG.at(2));
    }
  bool baseline = (jets->size() > 1
	&& recoV.mass() > 50 && recoV.mass() < 120
	&& KinematicsOK(recoV.daughter(0), ptcut,etacut)
	&& KinematicsOK(recoV.daughter(1), ptcut,etacut)
	&& Z->mass() > 60 && Z->mass() < 120
	&& KinematicsOK(Z->daughter(0), 5.,etacut)
	&& KinematicsOK(Z->daughter(1), 5.,etacut)
        && selectedPhotons.size()>0);

  switch (n)
  {
  case 1://baseline
    if (baseline)
      return true;
    break;

  case 2://baseline+mjj=mV+-15
    if (baseline
	&& recoV.mass() > 65 && recoV.mass() < 115)
      return true;
    break;

  case 3://baseline+mjj=mV+-15 & dRjj 2.4
    if (baseline
	&& abs(physmath::deltaR(recoV.daughter(0),recoV.daughter(1)))<2.4
	&& recoV.mass() > 65 && recoV.mass() < 115)
      return true;
    break;

  case 4://baseline+mjj=mV+-15 & dRjj 2.4 & 2D deltas 
    if (baseline
	//&& abs(physmath::deltaR(nearestRECOjetstoPhoton.first, nearestRECOjetstoPhoton.second))> 0.6 + (nearestRECOjetstoPhoton.first.pt()-nearestRECOjetstoPhoton.second.pt())/100
	&& abs(physmath::deltaR(nearestRECOjetstoPhoton.first, nearestRECOjetstoPhoton.second))> 0.1
	&& abs(physmath::deltaR(recoV.daughter(0),recoV.daughter(1)))<2.4
	&& recoV.mass() > 65 && recoV.mass() < 115)
      return true;
    break;

  case 5://baseline+mjj=mV+-15 & dRjj 2.4 & FWM T0 
    if (baseline
	&& abs(physmath::deltaR(nearestRECOjetstoPhoton.first, nearestRECOjetstoPhoton.second))> 0.1
	//	&& 1.6 < SumFWM(0, 't', lljjG) < 2
	&& abs(physmath::deltaR(recoV.daughter(0),recoV.daughter(1)))<2.4
	&& KinematicsOK(recoV.daughter(0), 40, etacut)
	&& recoV.mass() > 65 && recoV.mass() < 115)
      return true;
    break;


  case 6://baseline+mjj=mV+-15+ptj0>40 & dRjj 2.4 
    if (baseline
	&& KinematicsOK(recoV.daughter(0), 40,etacut)
	//&& KinematicsOK(recoV.daughter(1), 40,etacut)
	&& abs(physmath::deltaR(recoV.daughter(0),recoV.daughter(1)))<2.4
	&& recoV.mass() > 70 && recoV.mass() < 105)
      return true;
    break;
    

  case 7://baseline+mjj=mV+-20+ptj0>40 & dRjj 2.4
    if (baseline
	&& KinematicsOK(recoV.daughter(0), 40,etacut)
	&& KinematicsOK(recoV.daughter(1), 40,etacut)
	&& abs(physmath::deltaR(recoV.daughter(0),recoV.daughter(1)))<2.4
	&& recoV.mass() > 70 && recoV.mass() < 105
	&& abs(physmath::deltaR(nearestRECOjetstoPhoton.first, nearestRECOjetstoPhoton.second))> 1.67 + (nearestRECOjetstoPhoton.first.pt()-nearestRECOjetstoPhoton.second.pt())*3.67/15
	&& abs(physmath::deltaR(nearestRECOjetstoPhoton.first, nearestRECOjetstoPhoton.second))> 1.85 - (nearestRECOjetstoPhoton.first.pt()-nearestRECOjetstoPhoton.second.pt())*0.65/20)
      
      return true;
    break;

    
    /*    
  case 6://baseline+mjj=mV+-15+ptj>40 & dRjj 2.4 & > 2jets
    if (baseline
	&& jets->size() > 2
	&& KinematicsOK(recoV.daughter(0), 40,etacut)
	&& KinematicsOK(recoV.daughter(1), 40,etacut)
	&& abs(physmath::deltaR(recoV.daughter(0),recoV.daughter(1)))<2.4
	&& recoV.mass() > 65 && recoV.mass() < 115)
      return true;
    break;
  case 7://baseline+mjj=mV+-15+ptj>40 & dRjj 2.4 & > 2jets
    if (baseline
	&& jets->size() > 3
	&& KinematicsOK(recoV.daughter(0), 40,etacut)
	&& KinematicsOK(recoV.daughter(1), 40,etacut)
	&& abs(physmath::deltaR(recoV.daughter(0),recoV.daughter(1)))<2.4
	&& recoV.mass() > 65 && recoV.mass() < 115)
      return true;
    break;
    */
  default:
    return true;
  }
  return false;
}

void VZGAnalyzer::analyze()
{ // It's the only member function running each event.

  if(verbose==true)
    {
      cout << "----------------------------------------------------------------" << endl;
      cout << "Run: " << run << " event: " << event << endl;
    }

  genAnalyze();

  if (IN_GENsignalDef())    QuarksToJets();
  
  phys::Boson<phys::Jet> recoV;
  Reconstruct(&recoV);
  std::vector<phys::Photon> selectedphotons;
  PhotonSelection(&selectedphotons);

  // std::vector<phys::Jet> selectedRECOjets;
  // foreach (const phys::Jet &jet, *jets)
  // {
  //   if (KinematicsOK(jet, ptcut, etacut)) // KinematicsOK(jet)
  //   {
  //     selectedRECOjets.push_back(jet);
  //   }
  // }

  // phys::Photon mostEnergeticPhoton;
  //  if (selectedphotons.size() > 0)
  // {
  //   std::stable_sort(selectedphotons.begin(), selectedphotons.end(), phys::EComparator());
  //   mostEnergeticPhoton = selectedphotons[0];
  // }

  // theHistograms->fill("pt_mostenergeticphoton", "pt_mostenergeticphoton", 30, 0, 300, mostEnergeticPhoton.pt(), (theWeight*LumiSF));
  // std::vector<std::pair<phys::Photon, phys::Jet>> nearestRECOjetstoPhoton;

  // if (selectedRECOjets.size() > 0)
  // {
  //   std::stable_sort(selectedRECOjets.begin(), selectedRECOjets.end(), phys::DeltaRComparator(mostEnergeticPhoton));
  //   phys::Jet nearestRECOjet = selectedRECOjets.at(0);
  //   nearestRECOjetstoPhoton.push_back({mostEnergeticPhoton, nearestRECOjet});
  // }
  // for (auto pair : nearestRECOjetstoPhoton)
  // {
  //   ResolutionPlots(pair.first, pair.second, "Photon_vs_jet_", theWeight*LumiSF);
  //   theHistograms->fill("DeltaR_mostEnergeticPhoton_vs_BestMatchedRECOJet", "DeltaR_mostEnergeticPhoton_vs_BestMatchedRECOJet; #DeltaR", 50, 0, 5, abs(physmath::deltaR(pair.first, pair.second)), theWeight*LumiSF);
  //   theHistograms->fill("DeltaR_vs_Deltapt", "DeltaR_vs_Deltapt;#Delta pt [GeV/c] ; #DeltaR", 20, -100, 100, 20, 0, 2, pair.first.pt()-pair.second.pt(),abs(physmath::deltaR(pair.first, pair.second)), theWeight*LumiSF);
  // }

  /*
  std::vector<phys::Jet> kinRECOjets;
  foreach (const phys::Jet &jet, *jets)
    if (KinematicsOK(jet, ptcut, etacut)) // KinematicsOK(jet)
      kinRECOjets.push_back(jet);

  
  theHistograms->fill("#kinRECOjets", "#kinRECOjets", 8, 0, 8, kinRECOjets.size(), theWeight*LumiSF);
  theHistograms->fill("AtLeast_2kinRecoJets", "AtLeast_2kinRecoJets", 2, 0, 2, kinRECOjets.size()>1, theWeight*LumiSF);


  phys::Photon mostEnergeticPhoton;
  if (selectedphotons.size() > 0)
  {
    std::stable_sort(selectedphotons.begin(), selectedphotons.end(), phys::EComparator());
    mostEnergeticPhoton = selectedphotons[0];
  }

  theHistograms->fill("pt_mostenergeticphoton", "pt_mostenergeticphoton", 30, 0, 300, mostEnergeticPhoton.pt(), (theWeight*LumiSF));
  std::vector<std::pair<phys::Photon, phys::Jet>> nearestRECOjetstoPhoton;

  if (kinRECOjets.size() > 0)
  {
    std::stable_sort(kinRECOjets.begin(), kinRECOjets.end(), phys::DeltaRComparator(mostEnergeticPhoton));
    phys::Jet nearestRECOjet = kinRECOjets.at(0);
    nearestRECOjetstoPhoton.push_back({mostEnergeticPhoton, nearestRECOjet});
  }
  for (auto pair : nearestRECOjetstoPhoton)
  {
    //ResolutionPlots(pair.first, pair.second, "Photon_vs_jet_", theWeight*LumiSF);
    theHistograms->fill("DeltaR_mostEnergeticPhoton_vs_BestMatchedRECOJet", "DeltaR_mostEnergeticPhoton_vs_BestMatchedRECOJet; #DeltaR", 50, 0, 5, abs(physmath::deltaR(pair.first, pair.second)), theWeight*LumiSF);
    theHistograms->fill("DeltaR_vs_Deltapt", "DeltaR_vs_Deltapt;#Delta pt [GeV/c] ; #DeltaR", 20, -100, 100, 20, 0, 2, pair.first.pt()-pair.second.pt(),abs(physmath::deltaR(pair.first, pair.second)), theWeight*LumiSF);
  }
  */


  //---------------ALL-THE-EVENTS---------------//

 
  //  printHistos(0, "all", recoV,selectedphotons);

  //---------------SIGNAL-EVENTS---------------//

  if (IN_GENsignalDef())
  {
    // std::vector phys::Photon selectedphotons;
    // for (auto p : *photons)
    // {
    //   if (p.id() == 22 && KinematicsOK(p, 20, 2.4) && !p.hasPixelSeed() && p.passElectronVeto())
    //   {
    //     selectedphotons->push_back(p);
    //   }
    // }



    theHistograms->fill("Signal_V_vs_Z_pt", "Signal_V_vs_Z_pt;V p_{t} [GeV/c]; #Z p_{t} [GeV/c]", 30, 0, 300, 30, 0, 300, recoV.pt(), Z->pt(), theWeight*LumiSF);
    printHistos(0, "sign", recoV, selectedphotons);
  } // not signalConstraint anymore

  //---------------BACKGROUND-EVENTS---------------//

  else    printHistos(0, "bckg", recoV,selectedphotons);

}


void VZGAnalyzer::genPhotonsAnalyzer()
{}


void VZGAnalyzer::genVBAnalyzer()

{
  theHistograms->fill("nZtoChLep"    ,"Number of Z->ll per event" , 7,0,7, genVBHelper_.ZtoChLep().size());
  theHistograms->fill("nZtoNeutrinos","Number of Z->nn per event" , 7,0,7, genVBHelper_.ZtoNeutrinos().size());
  theHistograms->fill("nWtoLep"      ,"Number of W->lnu per event", 7,0,7, genVBHelper_.WtoLep().size());
  theHistograms->fill("nZtoQ"        ,"Number of Z->qq per event" , 7,0,7, genVBHelper_.ZtoQ().size());
  theHistograms->fill("nWtoQ"        ,"Number of W->qq' per ev""ent", 7,0,7, genVBHelper_.WtoQ().size());

  int nVBs = genVBHelper_.ZtoChLep().size() + genVBHelper_.ZtoNeutrinos().size() + genVBHelper_.WtoLep().size() + genVBHelper_.ZtoQ().size() + genVBHelper_.WtoQ().size();
  theHistograms->fill("nVBs", "Number of VB per event", 7,0,7, nVBs);
  for (auto VB : genVBHelper_.ZtoQ())
  {
    theHistograms->fill("ZtoQ mass", "ZtoQ mass", 10, 50, 120, VB.mass());
    theHistograms->fill("VB mass", "VB mass", 10, 50, 120, VB.mass());
  }
  for (auto VB : genVBHelper_.WtoQ())
  {
    theHistograms->fill("WtoQ mass", "WtoQ mass", 10, 50, 120, VB.mass());
    theHistograms->fill("VB mass", "VB mass", 10, 50, 120, VB.mass());
  }
  for (auto VB : genVBHelper_.ZtoChLep())
  {
    theHistograms->fill("ZtoChLep mass", "ZtoChLep mass", 10, 50, 120, VB.mass());
    theHistograms->fill("VB mass", "VB mass", 10, 50, 120, VB.mass());
  }
  for (auto VB : genVBHelper_.WtoLep())
  {
    theHistograms->fill("WtoLep mass", "WtoLep mass", 10, 50, 120, VB.mass());
    theHistograms->fill("VB mass", "VB mass", 10, 50, 120, VB.mass());
  }
  for (auto VB : genVBHelper_.ZtoNeutrinos())
  {
    theHistograms->fill("ZtoNeutrinos mass", "ZtoNeutrinos mass", 10, 50, 120, VB.mass());
    theHistograms->fill("VB mass", "VB mass", 10, 50, 120, VB.mass());
  }


  theHistograms->fill("nZtoChLep Alternative"    ,"Number of Z->ll per event Alternative" , 7,0,7, genZlepCandidates_->size());
  theHistograms->fill("nWtoLep Alternative"      ,"Number of W->lnu per event Alternative", 7,0,7, genWlepCandidates_->size());
  theHistograms->fill("nZtoQ Alternative"        ,"Number of Z->qq per event Alternative" , 7,0,7, genZhadCandidates_->size());
  theHistograms->fill("nWtoQ Alternative"        ,"Number of W->qq' per event Alternative", 7,0,7, genWlepCandidates_->size());

  int nVBs_alternative = genZlepCandidates_->size() + genWlepCandidates_->size() + genZhadCandidates_->size() + genWhadCandidates_->size();
  theHistograms->fill("nVBs Alternative", "Number of VB per event Alternative", 7,0,7, nVBs_alternative);
  foreach (auto & VB , *genZhadCandidates_)
  {
    theHistograms->fill("ZtoQ mass Alternative", "ZtoQ mass Alternative", 10, 50, 120, VB.mass());
    theHistograms->fill("VB mass Alternative", "VB mass Alternative", 10, 50, 120, VB.mass());
  }
  foreach (auto & VB , *genWhadCandidates_)
  {
    theHistograms->fill("WtoQ mass Alternative", "WtoQ mass Alternative", 10, 50, 120, VB.mass());
    theHistograms->fill("VB mass Alternative", "VB mass Alternative", 10, 50, 120, VB.mass());
  }
  foreach (auto & VB , *genZlepCandidates_)
  {
    theHistograms->fill("ZtoChLep mass Alternative", "ZtoChLep mass Alternative", 10, 50, 120, VB.mass());
    theHistograms->fill("VB mass Alternative", "VB mass Alternative", 10, 50, 120, VB.mass());
  }
  foreach (auto & VB , *genWlepCandidates_)
  {
    theHistograms->fill("WtoLep mass Alternative", "WtoLep mass Alternative", 10, 50, 120, VB.mass());
    theHistograms->fill("VB mass Alternative", "VB mass Alternative", 10, 50, 120, VB.mass());
  }
  // for (auto VB : genVBHelper_.ZtoNeutrinos())
  // {
  //   theHistograms->fill("ZtoNeutrinos mass Alternative", "ZtoNeutrinos mass Alternative", 10, 50, 120, VB.mass());
  //   theHistograms->fill("VB mass Alternative", "VB mass Alternative", 10, 50, 120, VB.mass());
  // }
}


void VZGAnalyzer::Reconstruct(phys::Boson<phys::Jet> *V_JJCandidate)
{
  // Building of every jets pairs combination
  std::vector<phys::Boson<phys::Jet>> DiJets;

  for (uint i = 0; i < jets->size(); i++) // Warning: size can be 0
    for (uint j = i+1; j < jets->size(); j++)
        DiJets.push_back(phys::Boson<phys::Jet>(jets->at(i), jets->at(j)));

  if (jets->size() > 1)
  {
    std::stable_sort(DiJets.begin(), DiJets.end(), phys::Mass2Comparator(phys::ZMASS, phys::WMASS));
    *V_JJCandidate = DiJets.at(0);
  }
}

void VZGAnalyzer::PhotonSelection(std::vector<phys::Photon> *phot)
{
  for (auto p : *photons)
  {
    // Pixel seed and electron veto
    // if (ph.hasPixelSeed() || !ph.passElectronVeto())
    //        continue;

    if (p.id() == 22 && KinematicsOK(p, 20, 2.4) && !p.hasPixelSeed() && p.passElectronVeto() && p.cutBasedIDLoose())
    {
        phot->push_back(p);
    }
  }
  if (verbose==true) std::cout << "Number of selected RECO photons = " << phot->size() << std::endl;
}


void VZGAnalyzer::CompatibilityTest(phys::Boson<phys::Jet> bestCandidate, phys::Boson<phys::Particle> genVB, std::string sample, std::string algorithm)
{
}

void VZGAnalyzer::printHistos(uint i, std::string histoType, phys::Boson<phys::Jet> recoV, std::vector<phys::Photon> selectedphotons)
{

  std::vector<std::string> cuts = {"0", "1", "2", "3", "4", "5", "6", "7"};//, "8", "9", "10", "11"};
  std::vector<std::string> orders = {"0", "1", "2", "3", "4", "5", "6", "7"};//, "8", "9", "10", "11"};
  
  if (i < cuts.size() && cut(i, recoV,selectedphotons))
  {

    theHistograms->fill("#AAA_cut_flow_" + histoType, "Cut flow", 12, 0, 12, i, (theWeight*LumiSF));

    theHistograms->fill("recoVMass_" + histoType + cuts.at(i), "mass of recoV", 40, 0, 200, recoV.mass(), (theWeight*LumiSF));
    //theHistograms->fill("recoVTot_" + histoType + cuts.at(i), "mass of recoV", 1, 50, 120, recoV.mass(), (theWeight*LumiSF));

    theHistograms->fill("recoVDaughter0Pt_" + histoType + cuts.at(i), "pt of recoVDaughter0", 50, 0, 600, recoV.daughter(0).pt(), (theWeight*LumiSF));
    theHistograms->fill("recoVDaughter1Pt_" + histoType + cuts.at(i), "pt of recoVDaughter1", 50, 0, 600, recoV.daughter(1).pt(), (theWeight*LumiSF));
    //PlotJets(recoV.daughter(0),recoV.daughter(1), "", theWeight*LumiSF, histoType + cuts.at(i));

    theHistograms->fill("recoZMass_" + histoType + cuts.at(i), "mass of recoZ", 40, 0, 200, Z->mass(), (theWeight*LumiSF));

    theHistograms->fill("recoVPt_" + histoType + cuts.at(i), "pt of recoV", 30, 0, 300, recoV.pt(), (theWeight*LumiSF));
    //theHistograms->fill("recoVEta_" + histoType + cuts.at(i), "eta of recoV", 30, 0, 3.5, fabs(recoV.eta()), (theWeight*LumiSF));
    //theHistograms->fill("recoVPhi_" + histoType + cuts.at(i), "phi of recoV", 30, -3.2, 3.2, recoV.phi(), (theWeight*LumiSF));
    //theHistograms->fill("recoVEnergy_" + histoType + cuts.at(i), "energy of  recoV", 60, 0, 2400, fabs(recoV.e()), (theWeight*LumiSF));
    theHistograms->fill("recoVDaughtersDeltaPhi_" + histoType + cuts.at(i), "dPhi of recoVDaughters", 20, 0, 3.2, fabs(physmath::deltaPhi(recoV.daughter(0).phi(), recoV.daughter(1).phi())), (theWeight*LumiSF));

    theHistograms->fill("recoZPt_" + histoType + cuts.at(i), "pt of recoZ", 50, 0, 600, Z->pt(), (theWeight*LumiSF));
    theHistograms->fill("recoZEta_" + histoType + cuts.at(i), "eta of recoZ", 35, 0, 3.5, fabs(Z->eta()), (theWeight*LumiSF));
    //theHistograms->fill("recoZEnergy_" + histoType + cuts.at(i), "energy of  recoZ", 120, 0, 400, fabs(Z->e()), (theWeight*LumiSF));
    theHistograms->fill("recoZDeltaPhi_" + histoType + cuts.at(i), "dPhi of recoZ", 30, 0, 3.2, fabs(physmath::deltaPhi(Z->daughter(0).phi(), Z->daughter(1).phi())), (theWeight*LumiSF));




    
  // std::vector<phys::Jet> selectedRECOjets;
  // foreach (const phys::Jet &jet, *jets)
  // {
  //   if (KinematicsOK(jet, ptcut, etacut)) // KinematicsOK(jet)
  //   {
  //     selectedRECOjets.push_back(jet);
  //   }
  // }

  // phys::Photon mostEnergeticPhoton;
  //  if (selectedphotons.size() > 0)
  // {
  //   std::stable_sort(selectedphotons.begin(), selectedphotons.end(), phys::EComparator());
  //   mostEnergeticPhoton = selectedphotons[0];
  // }

  // theHistograms->fill("pt_mostenergeticphoton"+histoType + cuts.at(i), "pt_mostenergeticphoton"+histoType + cuts.at(i), 30, 0, 300, mostEnergeticPhoton.pt(), (theWeight*LumiSF));
  // std::vector<std::pair<phys::Photon, phys::Jet>> nearestRECOjetstoPhoton;

  // if (selectedRECOjets.size() > 0)
  // {
  //   std::stable_sort(selectedRECOjets.begin(), selectedRECOjets.end(), phys::DeltaRComparator(mostEnergeticPhoton));
  //   phys::Jet nearestRECOjet = selectedRECOjets.at(0);
  //   nearestRECOjetstoPhoton.push_back({mostEnergeticPhoton, nearestRECOjet});
  // }
  // for (auto pair : nearestRECOjetstoPhoton)
  // {
  //   ResolutionPlots(pair.first, pair.second, "Photon_vs_jet_", theWeight*LumiSF,"_"+ histoType + cuts.at(i));
  //   theHistograms->fill("DeltaR_mostEnergeticPhoton_vs_BestMatchedRECOJet_"+histoType + cuts.at(i), "DeltaR_mostEnergeticPhoton_vs_BestMatchedRECOJet"+histoType + cuts.at(i)+"; #DeltaR", 50, 0, 5, abs(physmath::deltaR(pair.first, pair.second)), theWeight*LumiSF);
  //   theHistograms->fill("DeltaR_vs_Deltapt_"+histoType + cuts.at(i), "DeltaR_vs_Deltapt"+histoType + cuts.at(i)+";#Delta pt [GeV/c] ; #DeltaR", 20, -100, 100, 20, 0, 2, pair.first.pt()-pair.second.pt(),abs(physmath::deltaR(pair.first, pair.second)), theWeight*LumiSF);
  // }

    phys::Photon mostEnergeticPhoton;
    if (selectedphotons.size() > 0)
      {
	std::stable_sort(selectedphotons.begin(), selectedphotons.end(), phys::EComparator());
	mostEnergeticPhoton = selectedphotons[0];
      }

    theHistograms->fill("pt_mostenergeticphoton"+histoType + cuts.at(i), "pt_mostenergeticphoton"+histoType + cuts.at(i), 30, 0, 300, mostEnergeticPhoton.pt(), (theWeight*LumiSF));

    std::pair<phys::Photon, phys::Jet> nearestRECOjetstoPhoton;

    if(abs(physmath::deltaR(recoV.daughter(0), mostEnergeticPhoton))<abs(physmath::deltaR(recoV.daughter(1), mostEnergeticPhoton)) )
      nearestRECOjetstoPhoton={mostEnergeticPhoton, recoV.daughter(0)};
    else if (abs(physmath::deltaR(recoV.daughter(0), mostEnergeticPhoton))>abs(physmath::deltaR(recoV.daughter(1), mostEnergeticPhoton)) )
      nearestRECOjetstoPhoton={mostEnergeticPhoton, recoV.daughter(1)};

    theHistograms->fill("DR_gammaClosestJet_"+histoType + cuts.at(i), "DR_gammaClosestJet_"+histoType + cuts.at(i)+"; #DeltaR", 50, 0, 5, abs(physmath::deltaR(nearestRECOjetstoPhoton.first, nearestRECOjetstoPhoton.second)), theWeight*LumiSF);
    theHistograms->fill("DeltaR_vs_Deltapt_gammaJet"+histoType + cuts.at(i), "DeltaR_vs_Deltapt_gammaJet"+histoType + cuts.at(i)+";#Delta pt [GeV/c] ; #DeltaR", 20, -100, 100, 50, 0, 5, nearestRECOjetstoPhoton.first.pt()-nearestRECOjetstoPhoton.second.pt(),abs(physmath::deltaR(nearestRECOjetstoPhoton.first, nearestRECOjetstoPhoton.second)), theWeight*LumiSF);

    
    if(selectedphotons.size()>0)
      {
	std::vector<TLorentzVector> jjG;
	jjG.push_back(recoV.daughter(0).p4());
	jjG.push_back(recoV.daughter(1).p4());
	jjG.push_back(selectedphotons.at(0).p4());

	std::vector<TLorentzVector> lljjG;
	lljjG.push_back(Z->daughter(0).p4());
	lljjG.push_back(Z->daughter(1).p4());
	lljjG.push_back(jjG.at(0));
	lljjG.push_back(jjG.at(1));
	lljjG.push_back(jjG.at(2));

	for(int l = 0; l<3; l++)
	  {
    
	    theHistograms->fill("FWM_T"+orders.at(l)+"_jets_"+histoType+cuts.at(i), "FWM_T"+orders.at(l)+"_jets_"+histoType+cuts.at(i)+"; H_"+orders.at(l)+"^T jets + gamma", 40, -3, 3,
				SumFWM(l, 't', jjG), theWeight*LumiSF);

	    theHistograms->fill("FWM_T"+orders.at(l)+"_fullSyst_"+histoType+cuts.at(i), "FWM_T"+orders.at(l)+"_fullSyst_"+histoType+cuts.at(i)+"; H_"+orders.at(l)+"^T jets and gamma", 40, 0, 4,
				SumFWM(l, 't', lljjG), theWeight*LumiSF);
	  }
      }

    printHistos(++i, histoType, recoV, selectedphotons); 
  }
  return;
}



void VZGAnalyzer::genEventSetup(){
  genQuarks_->clear();
  genChLeptons_->clear();
  genNeutrinos_->clear();
  genPhotons_->clear();
  genPhotonsPrompt_->clear();
	
  genZlepCandidates_->clear();
  genWlepCandidates_->clear();
  genZhadCandidates_->clear();
  genWhadCandidates_->clear();
	
  genZZ_ = DiBoson<Particle, Particle>();
  genWZ_ = DiBoson<Particle, Particle>();
	
  // Sort gen particles
  for(auto p : *genParticles){
    unsigned int aPID = abs(p.id());
    if(aPID < 9)
      genQuarks_->push_back(p);
    else if(aPID == 11 || aPID == 13){
      genChLeptons_->push_back(p);
    }
    else if(aPID == 12 || aPID == 14)
      genNeutrinos_->push_back(p);
    else if(p.id() == 22){
      genPhotons_->push_back(p);
      if(p.genStatusFlags().test(phys::isPrompt))
	genPhotonsPrompt_->push_back(p);
    }
  }
	
  // Gen W --> l nu
  if(genNeutrinos_->size() > 0 && genChLeptons_->size() > 0){
    for(auto l : *genChLeptons_){
      for(auto v : *genNeutrinos_){
	if( abs(l.id() + v.id()) == 1 ){
	  Boson<Particle> Wcand(l,v);
	  if(GenWBosonDefinition(Wcand))
	    genWlepCandidates_->push_back(Wcand);
	}
      }
    }
  }
  
  // Gen Z --> l lbar
  if(genChLeptons_->size() >= 2){
    for(size_t i = 0 ; i < genChLeptons_->size(); ++i){
      Particle& l1 = genChLeptons_->at(i);
      for(size_t j = i+1; j < genChLeptons_->size(); ++j){
	Particle& l2 = genChLeptons_->at(j);
	
	if( l1.id() + l2.id() == 0 ){
	  Boson<Particle> Zcand(l1,l2);
	  if(ZBosonDefinition(Zcand))
	    genZlepCandidates_->push_back(Zcand);
	}
      }
    }
  }
  
  if(genQuarks_->size() >= 2){
    for(size_t i = 0  ; i < genQuarks_->size(); ++i){
      Particle& q1 = genQuarks_->at(i);
      if(q1.id() > 5) continue;
      for(size_t j = i+1; j < genQuarks_->size(); ++j){
	Particle& q2 = genQuarks_->at(j);
	if(q2.id() > 5) continue;

	// Gen W --> q q'bar
	if( (q1.id() * q2.id() < 0) && ( abs(q1.id()+q2.id()) % 2 ==1 ) ){
	  Boson<Particle> Wcand(q1,q2);
	  if(GenWBosonDefinition(Wcand))
	    genWhadCandidates_->push_back(Wcand);
	}

	// Gen Z --> q qbar
	if( q1.id() + q2.id() == 0 ){
	  Boson<Particle> Zcand(q1,q2);
	  if(ZBosonDefinition(Zcand))
	    genZhadCandidates_->push_back(Zcand);
	}
      }
    }
  }
  
  // genZZ --> 4l
  if(genChLeptons_->size() >= 4 && genZlepCandidates_->size() >= 2){
    std::sort(genZlepCandidates_->begin(), genZlepCandidates_->end(), MassComparator(phys::ZMASS));
    Boson<Particle>& Z0 = genZlepCandidates_->front();
		
    // Vector containing the rest of the Zll candidates
    std::vector<Boson<Particle>> Zll(genZlepCandidates_->begin()+1, genZlepCandidates_->end());
    std::sort(Zll.begin(), Zll.end(), ScalarSumPtComparator());
    Boson<Particle>* pZ1 = nullptr;
    for(size_t i = 0; i < Zll.size(); ++i){
      if(! haveCommonDaughter(Z0, Zll.at(i))){
	pZ1 = &(Zll.at(i));
	break;
      }
    }
    if(pZ1)
      genZZ_ = DiBoson<Particle, Particle>(Z0, *pZ1);
  }
	
  // genZW --> 3l nu
  if(genChLeptons_->size() >= 3 && genZlepCandidates_->size() >= 1 && genWlepCandidates_->size() >= 1){	
    std::sort(genZlepCandidates_->begin(), genZlepCandidates_->end(), MassComparator(phys::ZMASS));
    Boson<Particle>& Z0 = genZlepCandidates_->front();
		
    std::sort(genWlepCandidates_->begin(), genWlepCandidates_->end(), MassComparator(phys::WMASS));
    Boson<Particle>& W0 = genWlepCandidates_->front();
		
    genWZ_ = DiBoson<Particle, Particle>(Z0, W0);
  }
	
}

























void VZGAnalyzer::genAnalyze()
{
  theHistograms->fill("Signal_fraction_QLG", "Signal_fraction_QLG", 2, 0, 2, IN_GENsignalDef() , theWeight*LumiSF);
  theHistograms->fill("Signal_fraction_Q", "Signal_fraction_Q", 2, 0, 2, HadronicSignalConstraint(), theWeight*LumiSF);
  theHistograms->fill("Signal_fraction_L", "Signal_fraction_L", 2, 0, 2, LeptonicSignalConstraint(), theWeight*LumiSF);
  theHistograms->fill("Signal_fraction_L", "Signal_fraction_L", 2, 0, 2, PhotonSignalConstraint(), theWeight*LumiSF);

  //___________________________________________________________________________________STUDYING_LEPTONIC_CONSTRAINTS_ON_THE_SAMPLE_AT_GEN_LEVEL

  theHistograms->fill("LepSel_ALL_Zexists_ptl1OK_ptl2OK_ptlOK_mllOK", "LepSel_ALL_Zexists_ptl1OK_ptl2OK_ptlOK_mllOK", 6, 0, 6, 0., theWeight*LumiSF);

  if (genVBHelper_.ZtoChLep().size()==1){
    theHistograms->fill("LepSel_ALL_Zexists_ptl1OK_ptl2OK_ptlOK_mllOK", "LepSel_ALL_Zexists_ptl1OK_ptl2OK_ptlOK_mllOK", 6, 0, 6, 1., theWeight*LumiSF);
    if (KinematicsOK(genVBHelper_.ZtoChLep()[0].daughter(0), 5.,2.5))
      theHistograms->fill("LepSel_ALL_Zexists_ptl1OK_ptl2OK_ptlOK_mllOK", "LepSel_ALL_Zexists_ptl1OK_ptl2OK_ptlOK_mllOK", 6, 0, 6, 2., theWeight*LumiSF);
    if (KinematicsOK(genVBHelper_.ZtoChLep()[0].daughter(1), 5.,2.5))
      theHistograms->fill("LepSel_ALL_Zexists_ptl1OK_ptl2OK_ptlOK_mllOK", "LepSel_ALL_Zexists_ptl1OK_ptl2OK_ptlOK_mllOK", 6, 0, 6, 3., theWeight*LumiSF);
    if (KinematicsOK(genVBHelper_.ZtoChLep()[0].daughter(0), 5.,2.5) && KinematicsOK(genVBHelper_.ZtoChLep()[0].daughter(1), 5.,2.5))
      {
	theHistograms->fill("LepSel_ALL_Zexists_ptl1OK_ptl2OK_ptlOK_mllOK", "LepSel_ALL_Zexists_ptl1OK_ptl2OK_ptlOK_mllOK", 6, 0, 6, 4., theWeight*LumiSF);
	if (genVBHelper_.ZtoChLep()[0].mass()>60 && genVBHelper_.ZtoChLep()[0].mass()<120)
	  theHistograms->fill("LepSel_ALL_Zexists_ptl1OK_ptl2OK_ptlOK_mllOK", "LepSel_ALL_Zexists_ptl1OK_ptl2OK_ptlOK_mllOK", 6, 0, 6, 5., theWeight*LumiSF);
      }
  }
  //___________________________________________________________________________________


  
  bool GENsignal=(HadronicSignalConstraint() && LeptonicSignalConstraint() && PhotonSignalConstraint());
  bool IN_RECObaseline=baselineRequirements();

  theHistograms->fill("Signal_fraction_RECO", "Signal_fraction_RECO", 2, 0, 2, IN_RECObaseline, theWeight*LumiSF);

  //_____________________________________________________________________________________________//

  if (GENsignal && IN_RECObaseline)
    theHistograms->fill("GENRECO_trueBkg_sigLoss_fakeSig_trueSig", "GENRECO_trueBkg_sigLoss_fakeSig_trueSig", 4, 0, 4, 3., theWeight*LumiSF);
  else if (!GENsignal && IN_RECObaseline)
    theHistograms->fill("GENRECO_trueBkg_sigLoss_fakeSig_trueSig", "GENRECO_trueBkg_sigLoss_fakeSig_trueSig", 4, 0, 4, 2., theWeight*LumiSF);
  else if (GENsignal && !IN_RECObaseline)
    theHistograms->fill("GENRECO_trueBkg_sigLoss_fakeSig_trueSig", "GENRECO_trueBkg_sigLoss_fakeSig_trueSig", 4, 0, 4, 1., theWeight*LumiSF);
  else
    theHistograms->fill("GENRECO_trueBkg_sigLoss_fakeSig_trueSig", "GENRECO_trueBkg_sigLoss_fakeSig_trueSig", 4, 0, 4, 0., theWeight*LumiSF);
  //_____________________________________________________________________________________________//
  
  theHistograms->fill("GENRECO_11", "GENRECO_11", 2, 0, 2, GENsignal && IN_RECObaseline, theWeight*LumiSF);
  theHistograms->fill("GENRECO_01", "GENRECO_01", 2, 0, 2, !GENsignal && IN_RECObaseline, theWeight*LumiSF);
  theHistograms->fill("GENRECO_10", "GENRECO_10", 2, 0, 2, GENsignal && !IN_RECObaseline, theWeight*LumiSF);
  theHistograms->fill("GENRECO_00", "GENRECO_00", 2, 0, 2, !GENsignal && !IN_RECObaseline, theWeight*LumiSF);
}
//___________________________________________________________________________________
  // genVBAnalyzer();
  /*
  std::vector<phys::Boson<phys::Particle>> genV;
  //std::vector<phys::Boson<phys::Particle>> genV(genVBHelper_.ZtoQ().size()+genVBHelper_.WtoQ().size());
  if(genVBHelper_.ZtoQ().size()>0)
  {
    genV=genVBHelper_.ZtoQ();
    genV.insert(genV.end(), genVBHelper_.WtoQ().begin(), genVBHelper_.WtoQ().end());
  }
  else if(genVBHelper_.WtoQ().size()>0)
  {
         genV=genVBHelper_.WtoQ();
    genV.insert(genV.end(), genVBHelper_.ZtoQ().begin(), genVBHelper_.ZtoQ().end());
  }
    //---------------------------------------- Single q analysis ----------------------------------------//
  std::vector<phys::Particle> genQuarksfromV;
  for (auto VB : genV)
  {
    theHistograms->fill("quarkfromV charge", "quarkfromV charge", 7, -7. / 6., 7. / 6., VB.daughter(0).charge(), theWeight*LumiSF);
    theHistograms->fill("quarkfromV charge", "quarkfromV charge", 7, -7. / 6., 7. / 6., VB.daughter(1).charge(), theWeight*LumiSF);

    theHistograms->fill("quarkfromV pt", "quarkfromV pt", 50, 0, 600, VB.daughter(0).pt(), theWeight*LumiSF);
    theHistograms->fill("quarkfromV pt", "quarkfromV pt", 50, 0, 600, VB.daughter(1).pt(), theWeight*LumiSF);

    genQuarksfromV.push_back(VB.daughter(0));
    genQuarksfromV.push_back(VB.daughter(1));
  }
  theHistograms->fill("0size_GENQuarksfromV_beforecuts", "0size_GENQuarksfromV_beforecuts", 10, -0.5, 9.5, genQuarksfromV.size(), theWeight*LumiSF);
  genQuarksfromV.erase(std::remove_if(genQuarksfromV.begin(), genQuarksfromV.end(), [](phys::Particle p)
                                      { return !KinematicsOK(p, ptcut, etacut); }),
                       genQuarksfromV.end());
  theHistograms->fill("0size_GENQuarksfromV_aftercuts", "0size_GENQuarksfromV_aftercuts", 10, -0.5, 9.5, genQuarksfromV.size(), theWeight*LumiSF);
 
  //----------------------------------------Kinematic Cuts on VB(qq)----------------------------------------//
  std::vector<phys::Boson<phys::Particle>> DiQuarks=genV;

  DiQuarks.erase(std::remove_if(DiQuarks.begin(), DiQuarks.end(), [](phys::Boson<phys::Particle> VB)
                                      { return !(KinematicsOK(VB.daughter(0), ptcut, etacut)&&KinematicsOK(VB.daughter(1), ptcut, etacut)); }),
                       DiQuarks.end());

 //----------------------------------------Kinematic Cuts GEN Jets AK4 & RECO Jets AK4----------------------------------------//
  std::vector<phys::Particle> selectedGENjets;
  foreach (const phys::Particle &jet, *genJets)
  {
    if (KinematicsOK(jet,ptcut,etacut)) // KinematicsOK(jet)
    {
      selectedGENjets.push_back(jet);
    }
  }
  std::vector<phys::Jet> selectedRECOjets;
  foreach (const phys::Jet &jet, *jets)
  {
    if (KinematicsOK(jet,ptcut,etacut)) // KinematicsOK(jet)
    {
      selectedRECOjets.push_back(jet);
    }
  }

  //----------------------------------------Matching efficiency ______ SINGLE QUARK/SINGLE GENJET--------------//
  std::vector<phys::Particle> jetsfromquarks;
  std::vector<std::pair<phys::Particle, phys::Particle>> nearestjetstoquark;

  for (auto quark : genQuarksfromV)
  {
    phys::Particle nearestjet;
    bool makesjet = false;
    theHistograms->fill("Pt_quark_den", " Pt_quark_den; GeV/c", 10, ptcut, 300, quark.pt(), theWeight*LumiSF);

    if (selectedGENjets.size() > 0)
    {
      std::stable_sort(selectedGENjets.begin(), selectedGENjets.end(), phys::DeltaRComparator(quark));
      nearestjet = selectedGENjets.at(0);
      nearestjetstoquark.push_back({quark, nearestjet});
      if (abs(physmath::deltaR(quark, nearestjet)) < 0.4)
      {
        jetsfromquarks.push_back(nearestjet);
        makesjet = true;
      }
    }
    if (makesjet && selectedGENjets.size() > 0)
    {
      theHistograms->fill("#QUARK=>GEN", "#QUARK=>GEN", 2, 0, 2, 1., theWeight*LumiSF);
      theHistograms->fill("Pt_quark_num", " Pt_quark_num; GeV/c", 10, ptcut, 300, quark.pt(), theWeight*LumiSF);
    }
    else
    {
      theHistograms->fill("#QUARK=>GEN", "#QUARK=>GEN", 2, 0, 2, 0., theWeight*LumiSF);
    }
  }
  for (auto pair : nearestjetstoquark)
  {
    ResolutionPlots(pair.first,pair.second,"Hadronization_",theWeight*LumiSF,"");
    theHistograms->fill("DeltaR_quark_vs_BestMatchedGENJet", "DeltaR_quark_vs_BestMatchedGENJet; #DeltaR", 20, 0, 0.5, abs(physmath::deltaR(pair.first, pair.second)), theWeight*LumiSF);
    theHistograms->fill("DeltaR_quark_jet_vs_pt", "DeltaR vs pt;pt [GeV/c] ; #DeltaR", 10, ptcut, 300, 20, 0, 0.2, pair.first.pt(),abs(physmath::deltaR(pair.first, pair.second)), theWeight*LumiSF);

  }
  theHistograms->fill("1size_GENjetsfromquarks", "size_GENjetsfromquarks", 10, -0.5, 9.5, jetsfromquarks.size(), theWeight*LumiSF);

  //----------------------------------------Matching efficiency ______ SINGLE GENJET/SINGLE RECOJET--------------//
  std::vector<phys::Particle> RECOjetsfromGENjets;
  std::vector<std::pair<phys::Particle, phys::Particle>> nearestRECOjetstoGENjets;

  for (auto genJet : selectedGENjets)
  {
    phys::Particle nearestRECOjet;
    bool isreconstructed = false;
    theHistograms->fill("Pt_genJet_den", " Pt_genJet_den; GeV/c", 10, ptcut, 300, genJet.pt(), theWeight*LumiSF);

    if (selectedRECOjets.size() > 0)
    {
      std::stable_sort(selectedRECOjets.begin(), selectedRECOjets.end(), phys::DeltaRComparator(genJet));
      nearestRECOjet = selectedRECOjets.at(0);
      nearestRECOjetstoGENjets.push_back({genJet, nearestRECOjet});
      if (abs(physmath::deltaR(genJet, nearestRECOjet)) < 0.4)
      {
        RECOjetsfromGENjets.push_back(nearestRECOjet);
        isreconstructed = true;
      }
    }
    if (isreconstructed && selectedRECOjets.size() > 0)
    {
      theHistograms->fill("#GEN=>RECO", "#GEN=>RECO", 2, 0, 2, 1., theWeight*LumiSF);
      theHistograms->fill("Pt_genJet_num", " Pt_genJet_num; GeV/c", 10, ptcut, 300, genJet.pt(), theWeight*LumiSF);

    }
    else
    {
      theHistograms->fill("#GEN=>RECO", "#GEN=>RECO", 2, 0, 2, 0., theWeight*LumiSF);
    }
  }
  for (auto pair : nearestRECOjetstoGENjets)
  {
    ResolutionPlots(pair.first,pair.second,"SingleJetsReconstruction_",theWeight*LumiSF,"");
    theHistograms->fill("DeltaR_GENjet_vs_BestMatchedRECOJet", "DeltaR_GENjet_vs_BestMatchedRECOJet; #DeltaR", 20, 0, 0.5, abs(physmath::deltaR(pair.first, pair.second)), theWeight*LumiSF);
    theHistograms->fill("DeltaR_jets_vs_pt", "DeltaR jets vs pt;pt [GeV/c] ; #DeltaR", 10, ptcut, 300, 20, 0, 0.2, pair.first.pt(),abs(physmath::deltaR(pair.first, pair.second)), theWeight*LumiSF);
  }
  theHistograms->fill("2size_GENjetsRECONSTRUCTED", "2size_GENjetsRECONSTRUCTED", 10, -0.5, 9.5, RECOjetsfromGENjets.size(), theWeight*LumiSF);

  //----------------------------------------Matching efficiency ______ QUARKS PAIR------------------------//

  std::vector<std::pair<phys::Particle, phys::Particle>> DijetsmatchedtoDiquark;
  std::vector<phys::Boson<phys::Particle>> DiJetsGEN;
  std::cout << ".................GEN to QUARKS MATCHING..............." << std::endl;

  std::cout << "DiQuarks size: " << DiQuarks.size() << std::endl;
  for (auto Diquark : DiQuarks)
  {

    bool firstmatches = false;
    phys::Particle jetmatchedtoFIRSTquark;

    bool secondmatches = false;
    phys::Particle jetmatchedtoSECONDquark;

    bool atleastonematches = false;
    bool bothmatch = false;
    std::cout << "selectedGENjets size: " << selectedGENjets.size() << std::endl;
    for (auto genJet : selectedGENjets)
    {
      double deltaR1 = abs(physmath::deltaR(Diquark.daughter(0), genJet));
      std::cout << "deltaR1= " << deltaR1 << std::endl;
      double deltaR2 = abs(physmath::deltaR(Diquark.daughter(1), genJet));
      std::cout << "deltaR2= " << deltaR2 << std::endl;

      if (deltaR1 < 0.4)
      {
        std::cout << "first matched" << std::endl;
        jetmatchedtoFIRSTquark = genJet;
        firstmatches = true;
      }
      if (deltaR2 < 0.4)
      {
        jetmatchedtoSECONDquark = genJet;
        std::cout << "second matched" << std::endl;
        secondmatches = true;
      }
    }

    bothmatch = (firstmatches && secondmatches);
    atleastonematches = (firstmatches || secondmatches);

    std::cout << "bothmatch= " << bothmatch << std::endl;
    std::cout << "atleastonematches= " << atleastonematches << std::endl;

    if (bothmatch)
    {
      std::cout << "both matched" << std::endl;
      std::cout << "reconstructing a boson from dijets matched to diquark" << std::endl;
      DiJetsGEN.push_back(phys::Boson<phys::Particle>(jetmatchedtoFIRSTquark, jetmatchedtoSECONDquark));
      DijetsmatchedtoDiquark.push_back({Diquark, phys::Boson<phys::Particle>(jetmatchedtoFIRSTquark, jetmatchedtoSECONDquark)});
      theHistograms->fill("#Bothmatched", "#Bothmatched", 2, 0, 2, 1., theWeight*LumiSF);
    }
    if (!bothmatch)
    {
      theHistograms->fill("#Bothmatched", "#Bothmatched", 2, 0, 2, 0., theWeight*LumiSF);
    }

    if (atleastonematches)
    {
      theHistograms->fill("#AtLeastONEmatches", "#AtLeastONEmatches", 2, 0, 2, 1., theWeight*LumiSF);
    }
    if (!atleastonematches)
    {
      theHistograms->fill("#AtLeastONEmatches", "#AtLeastONEmatches", 2, 0, 2, 0., theWeight*LumiSF);
    }
    if (atleastonematches && bothmatch)
    {
      theHistograms->fill("#Bothmatched|AtLeastONEmatches", "#Bothmatched|AtLeastONEmatches", 2, 0, 2, 1., theWeight*LumiSF);
    }
    if (atleastonematches && !bothmatch)
    {
      theHistograms->fill("#Bothmatched|AtLeastONEmatches", "#Bothmatched|AtLeastONEmatches", 2, 0, 2, 0., theWeight*LumiSF);
    }
  }
  theHistograms->fill("1.1size_GEN_DiJets", "1size_GEN_DiJets", 10, -0.5, 9.5, DiJetsGEN.size(), theWeight*LumiSF);
  for (auto DiJet : DiJetsGEN)
  {
    float mjj = (DiJet.daughter(0).p4() + DiJet.daughter(1).p4()).M();
    theHistograms->fill("mjj_GEN", "mjj_GEN", 10, 50, 120, mjj, theWeight*LumiSF);
  }

  //----------------------------------------Matching efficiency ______ RECO to GEN  PAIR------------------------//
  std::cout << ".................RECO to GEN MATCHING..............." << std::endl;

  // std::vector<std::pair<phys::Particle,phys::Particle>> DiRECOjetsmatchedtoDiGENjets;
  // std::vector<phys::Boson<phys::Particle>> DiJetsRECO;
  std::vector<phys::Boson<phys::Particle>> DiJetsGENreconstructed;
  std::cout << "GEN Dijets size: " << DiJetsGEN.size() << std::endl;
  for (auto DiJet : DiJetsGEN)
  {

    bool firstmatches = false;
    // phys::Particle jetmatchedtoFIRSTgen;

    bool secondmatches = false;
    // phys::Particle jetmatchedtoSECONDgen;

    bool atleastonematches = false;
    bool bothmatch = false;
    std::cout << "selectedRECOjets size: " << selectedRECOjets.size() << std::endl;
    for (auto recoJet : selectedRECOjets)
    {
      double deltaR1 = abs(physmath::deltaR(DiJet.daughter(0), recoJet));
      std::cout << "deltaR1= " << deltaR1 << std::endl;
      double deltaR2 = abs(physmath::deltaR(DiJet.daughter(1), recoJet));
      std::cout << "deltaR2= " << deltaR2 << std::endl;

      if (deltaR1 < 0.4)
      {
        std::cout << "first matched" << std::endl;
        // jetmatchedtoFIRSTgen=recoJet;
        firstmatches = true;
      }
      if (deltaR2 < 0.4)
      {
        std::cout << "second matched" << std::endl;
        // jetmatchedtoSECONDgen=recoJet;
        secondmatches = true;
      }
    }

    bothmatch = (firstmatches && secondmatches);
    atleastonematches = (firstmatches || secondmatches);

    std::cout << "bothmatch= " << bothmatch << std::endl;
    std::cout << "atleastonematches= " << atleastonematches << std::endl;

    if (bothmatch)
    {
      std::cout << "both matched" << std::endl;
      // std::cout << "reconstructing a boson from diRECOjets matched to diGENjets" << std::endl;
      // DiJetsRECO.push_back(phys::Boson<phys::Particle>(jetmatchedtoFIRSTgen, jetmatchedtoSECONDgen));
      DiJetsGENreconstructed.push_back(DiJet);
      // DiRECOjetsmatchedtoDiGENjets.push_back({DiJet,phys::Boson<phys::Particle>(jetmatchedtoFIRSTgen, jetmatchedtoSECONDgen)});

      theHistograms->fill("#RECOGEN_Bothmatched", "#RECOGEN_Bothmatched", 2, 0, 2, 1., theWeight*LumiSF);
    }
    if (!bothmatch)
    {
      theHistograms->fill("#RECOGEN_Bothmatched", "#RECOGEN_Bothmatched", 2, 0, 2, 0., theWeight*LumiSF);
    }

    if (atleastonematches)
    {
      theHistograms->fill("#RECOGEN_AtLeastONEmatches", "#RECOGEN_AtLeastONEmatches", 2, 0, 2, 1., theWeight*LumiSF);
    }
    if (!atleastonematches)
    {
      theHistograms->fill("#RECOGEN_AtLeastONEmatches", "#RECOGEN_AtLeastONEmatches", 2, 0, 2, 0., theWeight*LumiSF);
    }
    if (atleastonematches && bothmatch)
    {
      theHistograms->fill("#RECOGEN_Bothmatched|AtLeastONEmatches", "#RECOGEN_Bothmatched|AtLeastONEmatches", 2, 0, 2, 1., theWeight*LumiSF);
    }
    if (atleastonematches && !bothmatch)
    {
      theHistograms->fill("#RECOGEN_Bothmatched|AtLeastONEmatches", "#RECOGEN_Bothmatched|AtLeastONEmatches", 2, 0, 2, 0., theWeight*LumiSF);
    }
  }
  // theHistograms->fill("2size_RECO_DiJets", "1size_RECO_DiJets", 10, -0.5, 9.5, DiJetsRECO.size(), theWeight*LumiSF);
  //  for (auto DiJet:DiJetsRECO)
  //  {
  //      float mjj = (DiJet.daughter(0).p4() + DiJet.daughter(1).p4()).M();
  //      theHistograms->fill("mjj_RECO", "mjj_RECO", 10, 50, 120, mjj, theWeight*LumiSF);
  //  }

  //----------------------------------------Matching efficiency ______ ALGORITHM------------------------//
  std::cout << ".................Algorithm efficiency..............." << std::endl;

  std::vector<std::pair<phys::Particle, phys::Particle>> DiRECOjetsmatchedtoDiGENjets;
  std::vector<phys::Boson<phys::Particle>> DiJetsRECO;


// reconstructing pairs


    std::map<std::string, Boson<phys::Jet>> Candidates;

    std::vector<phys::Boson<phys::Jet>> JetPairs;

    for (size_t i = 0; i < selectedRECOjets.size(); i++)
    {

      for (size_t j = i + 1; j < selectedRECOjets.size(); j++)
      {
        JetPairs.push_back(phys::Boson<phys::Jet>(selectedRECOjets.at(i), selectedRECOjets.at(j)));
      }
    }

    theHistograms->fill("2size_Jetspairs_RECO", "size_Jetsparis_RECO", 10, -0.5, 9.5, JetPairs.size(), theWeight*LumiSF);
    std::cout << "#reco jets pairs : " << JetPairs.size() << std::endl;

    if (JetPairs.size() > 0)
    {

      // 1st reconstruction model: comparison with WMass
      std::stable_sort(JetPairs.begin(), JetPairs.end(), phys::MassComparator(phys::WMASS));
      Candidates["mW"] = JetPairs.at(0);

      //2nd reconstruction model: comparison with ZMass
      std::stable_sort(JetPairs.begin(), JetPairs.end(), phys::MassComparator(phys::ZMASS));
      Candidates["mZ"] = JetPairs.at(0);

      // 3rd reconstruction model: maximization of candidate Pt
      std::stable_sort(JetPairs.begin(), JetPairs.end(), phys::ScalarSumPtComparator());
      Candidates["maxVPt"] = JetPairs.at(0);

      // // // 4th reconstruction model: minimization of total Pt of Zjj system
      // std::vector<phys::Particle> ZZjj;
      // phys::Particle ZZjjCandidate;
      // for (uint i = 0; i < JetPairs.size(); i++)
      // {
      //   phys::Particle totState(ZZ->p4() + (JetPairs.at(i)).p4());
      //   ZZjj.push_back(totState.p4());
      // }
      // std::stable_sort(ZZjj.begin(), ZZjj.end(), phys::PtComparator());
      // ZZjjCandidate = ZZjj.back();
      // for (uint i = 0; i < JetPairs.size(); i++)
      //   if ((JetPairs.at(i)).p4() == (ZZjjCandidate.p4() - ZZ->p4()))
      //     Candidates["minTotPt"] = JetPairs.at(i);

      // 4th reconstruction model: comparison with a mean value between ZMass and WMass
      std::stable_sort(JetPairs.begin(), JetPairs.end(), phys::MassComparator(0.2 * phys::ZMASS + 0.8 * phys::WMASS));
      Candidates["m8W2Z"] = JetPairs.at(0);
      // 5th reconstruction model: comparison with a mean value between ZMass and WMass

      std::stable_sort(JetPairs.begin(), JetPairs.end(), phys::Mass2Comparator(phys::ZMASS, phys::WMASS));
      Candidates["mWZ"] = JetPairs.at(0);

      for (auto Candidate : Candidates)
      {
        theHistograms->fill("mjj_" + Candidate.first + "_Candidate", "mjj_" + Candidate.first + "_Candidate", 10, 50, 120, Candidate.second.mass(), theWeight*LumiSF);
      }
    }
  
  theHistograms->fill("2size_RECO_JetPairs", "2size_RECO_JetPairs", 10, -0.5, 9.5, JetPairs.size(), theWeight*LumiSF);


  std::cout << "#true gen jets pairs reconstructed: " << DiJetsGENreconstructed.size() << std::endl;

  for (auto DiJet : DiJetsGENreconstructed)
  {
    theHistograms->fill("mjj_den" , " mjj_den; GeV/c^{2}" , 10, 50, 120, DiJet.mass(), theWeight*LumiSF);
    theHistograms->fill("Pt_den" , " Pt_den; GeV/c" , 10, 0, 300, DiJet.pt(), theWeight*LumiSF);
    std::cout<<""<<std::endl;


    //----------------------------------------MATCHED Jets total mass----------------------------------//
    for (auto Candidate : Candidates)
    {

        std::cout<<""<<std::endl;

        bool truepair = false;
        std::cout << "Algorithm: " << Candidate.first << std::endl;
        phys::Jet jetRECOA = Candidate.second.daughter(0);
        phys::Jet jetRECOB = Candidate.second.daughter(1);
        phys::Particle jetGENA = DiJet.daughter(0);
        phys::Particle jetGENB = DiJet.daughter(1);
        double deltaRAA = abs(physmath::deltaR(jetGENA, jetRECOA));
        double deltaRAB = abs(physmath::deltaR(jetGENA, jetRECOB));
        double deltaRBA = abs(physmath::deltaR(jetGENB, jetRECOA));
        double deltaRBB = abs(physmath::deltaR(jetGENB, jetRECOB));
        if ((deltaRAA < 0.4 && deltaRBB < 0.4) || (deltaRAB < 0.4 && deltaRBA < 0.4))
        {
          truepair = true;
        }
        if (truepair)
        {
          std::cout << "the algorithm selected a reco pair matched to a gen pair" << std::endl;
          theHistograms->fill("#Algorithm_" + Candidate.first, "#Algorithm_" + Candidate.first, 2, 0, 2, 1., theWeight*LumiSF);
          theHistograms->fill("PASSED mjj_" + Candidate.first + "_Candidate", " PASSED mjj_" + Candidate.first + "_Candidate", 10, 50, 120, DiJet.mass(), theWeight*LumiSF);
          theHistograms->fill("mjj_" + Candidate.first + "_num", " mjj_" + Candidate.first + "_num; GeV/c^{2}", 10, 50, 120, DiJet.mass(), theWeight*LumiSF);
          theHistograms->fill("Pt_" + Candidate.first + "_num", " Pt_" + Candidate.first + "_num; GeV/c", 10, 0, 300, DiJet.pt(), theWeight*LumiSF);

        }
        else
        {
          theHistograms->fill("FAILED mjj_" + Candidate.first + "_Candidate", " FAILED mjj_" + Candidate.first + "_Candidate", 10, 50, 120, DiJet.mass(), theWeight*LumiSF);
          std::cout << "the algorithm selected a reco pair NOT matched to a gen pair" << std::endl;
          theHistograms->fill("#Algorithm_" + Candidate.first, "#Algorithm_" + Candidate.first, 2, 0, 2, 0., theWeight*LumiSF);
        }
        if (truepair && Candidate.first=="mWZ")
        {
          ResolutionPlots(DiJet,Candidate.second,"VectorBosonReconstruction_",theWeight*LumiSF,"");
        }
    }
  }

}
  */
void VZGAnalyzer::ResolutionPlots(const phys::Particle &gen, const phys::Particle &reco, std::string prename, const float weight, std::string suffix)
{
  std::string where;
  if (fabs(gen.eta()) < 2.4)
  {
    where = "Barrel";
  }
  else if (fabs(gen.eta()) > 2.4 && fabs(gen.eta()) < 4.7)
  {
    where = "Endcap";
  }
  double delta_charge = (reco.charge() - gen.charge());
  double delta_mass = (reco.mass() - gen.mass());
  double delta_trmass = (reco.p4().Mt() - gen.p4().Mt());
  double delta_pt = (reco.pt() - gen.pt());
  double delta_eta = (reco.eta() - gen.eta());
  double delta_phi = (reco.phi() - gen.phi());
  // double JJdeltaR = abs(physmath::deltaR(gen, reco));

  double res_mass = (reco.mass() - gen.mass()) / gen.mass();
  double res_trmass = (reco.p4().Mt() - gen.p4().Mt()) / gen.p4().Mt();
  double res_pt = (reco.pt() - gen.pt()) / gen.pt();

  if (where == "Barrel" || where == "Endcap")
  {
    std::string name = "Delta";
    theHistograms->fill(prename + name + "_charge_" + where + suffix, prename + "#" + name + "_charge_" + where + ";#Delta charge", 9, -4.5, 4.5, delta_charge, weight);
    theHistograms->fill(prename + name + "_mass_" + where + suffix, prename + "#" + name + "_mass_" + where + ";#Delta mass [GeV/c^2]", 10, -20, 20, delta_mass, weight);
    theHistograms->fill(prename + name + "_trmass_" + where + suffix, prename + "#" + name + "_trmass_" + where + ";#Delta Tr mass [GeV/c^2]", 10, -20, 20, delta_trmass, weight);
    theHistograms->fill(prename + name + "_pt_" + where + suffix, prename + "#" + name + "_pt_" + where + ";#Delta pt [GeV/c]", 10, -20, 20, delta_pt, weight);

    theHistograms->fill(prename + name + "_eta_" + where + suffix, prename + "#" + name + "_eta_" + where + ";#Delta eta", 10, -1, 1, delta_eta, weight);
    theHistograms->fill(prename + name + "_phi_" + where + suffix, prename + "#" + name + "_#phi_" + where + ";#Delta phi", 10, -1, 1, delta_phi, weight);

    theHistograms->fill(prename + name + "_mass_vs_mass_" + where + suffix, prename + "#" + name + "_mass_vs_mass_" + where + ";mass[GeV/c^2]; #Delta mass [GeV/c^2]", 20, 0, 200, 10, -20, 20, gen.mass(), delta_mass, weight);
    theHistograms->fill(prename + name + "_trmass_vs_trmass_" + where + suffix, prename + "#" + name + "_trmass_vs_trmass_" + where + ";Trmass[GeV/c^2]; #Delta Trmass [GeV/c^2]", 20, 0, 200, 10, -20, 20, gen.p4().Mt(), delta_trmass, weight);
    theHistograms->fill(prename + name + "_pt_vs_pt_" + where + suffix, prename + "#" + name + "_pt_vs_p_{t}_" + where + ";p_{t} [GeV/c^2]; #Delta p_{t} [GeV/c^2]", 30, 0, 300, 10, -20, 20, gen.pt(), delta_pt, weight);

    name = "Res";
    //  theHistograms->fill(prename + name + "_charge_" + suffix, prename + name + "'s charge", 5, -2.5, 2.5, Jet1.charge(), weight);
    theHistograms->fill(prename + name + "_mass_" + where + suffix, prename + name + "_mass" + where, 10, -1, 1, res_mass, weight);
    theHistograms->fill(prename + name + "_trmass_" + where + suffix, prename + name + "_trmass_" + where, 10, -1, 1, res_trmass, weight);
    theHistograms->fill(prename + name + "_pt_" + where + suffix, prename + name + "_p_{t}_" + where, 10, -1, 1, res_pt, weight);
    //  theHistograms->fill(prename + name + "_Y_" + suffix, prename + name + "'s Y", 50, -5, 5, Jet1.rapidity(), weight);
    //  theHistograms->fill(prename + name + "_eta_" + suffix, prename + name + "'s #eta", 50, -9, 9, Jet1.eta(), weight);
    //  theHistograms->fill(prename + name + "_phi_" + suffix, prename + name + "'s #phi", 50, -3.5, 3.5, Jet1.phi(), weight);
  }
}




void VZGAnalyzer::QuarksToJets()
{
  std::vector<phys::Particle> genQuarks;
  foreach (const Particle &p, *genParticles)
    {
      if  (abs(p.id()) < 10) // Is it a quark? 
	{
	  theHistograms->fill("quark charge", "quark charge", 7, -7. / 6., 7. / 6., p.charge(), theWeight*LumiSF);
	  theHistograms->fill("quark pt", "quark pt", 50, 0, 600, p.pt(), theWeight*LumiSF);
	  genQuarks.push_back(Particle(p));
		     

	}
    }
    //---------------------------------------- Single q analysis and cuts ----------------------------------------//
  // std::vector<phys::Particle> genQuarksfromV;
  // for (auto VB : genV)
  // {
  //   theHistograms->fill("quarkfromV charge", "quarkfromV charge", 7, -7. / 6., 7. / 6., VB.daughter(0).charge(), theWeight*LumiSF);
  //   theHistograms->fill("quarkfromV charge", "quarkfromV charge", 7, -7. / 6., 7. / 6., VB.daughter(1).charge(), theWeight*LumiSF);

  //   theHistograms->fill("quarkfromV pt", "quarkfromV pt", 50, 0, 600, VB.daughter(0).pt(), theWeight*LumiSF);
  //   theHistograms->fill("quarkfromV pt", "quarkfromV pt", 50, 0, 600, VB.daughter(1).pt(), theWeight*LumiSF);

  //   genQuarksfromV.push_back(VB.daughter(0));
  //   genQuarksfromV.push_back(VB.daughter(1));
  // }
  // theHistograms->fill("0size_GENQuarksfromV_beforecuts", "0size_GENQuarksfromV_beforecuts", 10, -0.5, 9.5, genQuarksfromV.size(), theWeight*LumiSF);
  // genQuarksfromV.erase(std::remove_if(genQuarksfromV.begin(), genQuarksfromV.end(), [](phys::Particle p)
  //                                     { return !KinematicsOK(p, ptcut, etacut); }),
  //                      genQuarksfromV.end());
  // theHistograms->fill("0size_GENQuarksfromV_aftercuts", "0size_GENQuarksfromV_aftercuts", 10, -0.5, 9.5, genQuarksfromV.size(), theWeight*LumiSF);
 

 //----------------------------------------Kinematic Cuts GEN Jets AK4 & RECO Jets AK4----------------------------------------//
  std::vector<phys::Particle> selectedGENjets;

  theHistograms->fill("#genJets_overall", "#genJets_overall", 8, 0, 8, genJets->size(), theWeight*LumiSF);

  foreach (const phys::Particle &jet, *genJets)
    {
      if (KinematicsOK(jet,ptcut,etacut))
	selectedGENjets.push_back(jet);
    }
  theHistograms->fill("#genJets_selected", "#genJets_selected", 8, 0, 8, selectedGENjets.size(), theWeight*LumiSF);
  theHistograms->fill("AtLeast_2GenJets", "AtLeast_2GenJets", 2, 0, 2, selectedGENjets.size()>1, theWeight*LumiSF);

  theHistograms->fill("#recoJets_overall", "#recoJets_overall", 8, 0, 8, jets->size(), theWeight*LumiSF);

  std::vector<phys::Jet> selectedRECOjets;
  foreach (const phys::Jet &jet, *jets)
    {
    if (KinematicsOK(jet,ptcut,etacut))
      selectedRECOjets.push_back(jet);
  }
  theHistograms->fill("#recoJets_selected", "#recoJets_selected", 8, 0, 8, selectedRECOjets.size(), theWeight*LumiSF);
  theHistograms->fill("AtLeast_2RecoJets", "AtLeast_2RecoJets", 2, 0, 2, selectedRECOjets.size()>1, theWeight*LumiSF);


  //----------------------------------------Matching efficiency ______ SINGLE QUARK/SINGLE GENJET--------------//
  std::vector<phys::Particle> GENjetsfromquarks;
  std::vector<std::pair<phys::Particle, phys::Particle>> nearestGENjetstoquark;

  phys::Particle firstGENjetMatched;
  
  int quarkMatchingCounter = 0;
  int twoQuarksMatched = 0;
  
  for (auto quark : genQuarks)
  {
    phys::Particle nearestGENjet;
    bool makesGENjet = false;

    theHistograms->fill("Pt_quark_den", " Pt_quark_den; GeV/c", 10, 0, 300, quark.pt(), theWeight*LumiSF);

    if (selectedGENjets.size() > 0)
    {
      std::stable_sort(selectedGENjets.begin(), selectedGENjets.end(), phys::DeltaRComparator(quark));
      nearestGENjet = selectedGENjets.at(0);
      nearestGENjetstoquark.push_back({quark, nearestGENjet});
      if (abs(physmath::deltaR(quark, nearestGENjet)) < 0.4)
      {
        GENjetsfromquarks.push_back(nearestGENjet);
        makesGENjet = true;
	if (quarkMatchingCounter == 0) firstGENjetMatched = nearestGENjet;
        else	theHistograms->fill("overlappedQuarks", "overlappedQuarks", 2, 0, 2, quarkMatchingCounter == 1 && abs(physmath::deltaR(firstGENjetMatched, nearestGENjet))<0.4, theWeight*LumiSF);
	quarkMatchingCounter++;
      }
      else       theHistograms->fill("dR_unmatchedQuarks_closestGENjet", "dR_unmatchedQuarks_closestGENjet", 50, 0, 5, abs(physmath::deltaR(quark, nearestGENjet)), theWeight*LumiSF);


    }
    if (makesGENjet && selectedGENjets.size() > 0)
    {
      theHistograms->fill("#QUARK=>GEN", "#QUARK=>GEN", 2, 0, 2, 1., theWeight*LumiSF);
      theHistograms->fill("Pt_quark_num", " Pt_quark_num; GeV/c", 10, 0, 300, quark.pt(), theWeight*LumiSF);
    }
    else
    {
      theHistograms->fill("#QUARK=>GEN", "#QUARK=>GEN", 2, 0, 2, 0., theWeight*LumiSF);
    }
  }

  if(quarkMatchingCounter==2) twoQuarksMatched =1;
  else   theHistograms->fill("#events with 1 quark matching", "#events with 1 quark matching", 2, 0, 2, quarkMatchingCounter, theWeight*LumiSF);
  theHistograms->fill("#events with 2 quarks matching", "#events with 2 quarks matching", 2, 0, 2, twoQuarksMatched, theWeight*LumiSF);
  theHistograms->fill("#quarks matching", "#quarks matching", 3, 0, 3, quarkMatchingCounter, theWeight*LumiSF);
      
  
  for (auto pair : nearestGENjetstoquark)
  {
    ResolutionPlots(pair.first,pair.second,"Hadronization_",theWeight*LumiSF,"");
    theHistograms->fill("DeltaR_quark_vs_BestMatchedGENJet", "DeltaR_quark_vs_BestMatchedGENJet; #DeltaR", 20, 0, 0.5, abs(physmath::deltaR(pair.first, pair.second)), theWeight*LumiSF);
    theHistograms->fill("DeltaR_quark_GENjet_vs_pt", "DeltaR vs pt;pt [GeV/c] ; #DeltaR", 10, 0, 300, 20, 0, 0.2, pair.first.pt(),abs(physmath::deltaR(pair.first, pair.second)), theWeight*LumiSF);

  }
  theHistograms->fill("1size_GENjetsfromquarks", "size_GENjetsfromquarks", 10, -0.5, 9.5, GENjetsfromquarks.size(), theWeight*LumiSF);


  //----------------------------------------Matching efficiency ______ SINGLE GENJET/SINGLE RECOJET--------------//
  std::vector<phys::Particle> RECOjetsfromGENjets;
  std::vector<std::pair<phys::Particle, phys::Particle>> nearestRECOjetstoGENjets;

  int GENtoRECOjetsCounter=0;

  for (auto genJet : GENjetsfromquarks)
  {
    phys::Particle nearestRECOjet;
    bool isreconstructed = false;
    theHistograms->fill("Pt_genJet_den", " Pt_genJet_den; GeV/c", 10, 0, 300, genJet.pt(), theWeight*LumiSF);

    if (selectedRECOjets.size() > 0)
    {
      std::stable_sort(selectedRECOjets.begin(), selectedRECOjets.end(), phys::DeltaRComparator(genJet));
      nearestRECOjet = selectedRECOjets.at(0);
      nearestRECOjetstoGENjets.push_back({genJet, nearestRECOjet});
      theHistograms->fill("dR_GENclosestRECOjet", "dR_GENclosestRECOjet", 500, 0, 5, abs(physmath::deltaR(genJet, nearestRECOjet)), theWeight*LumiSF);
      if (abs(physmath::deltaR(genJet, nearestRECOjet)) < 0.4)
      {
        RECOjetsfromGENjets.push_back(nearestRECOjet);
        isreconstructed = true;
	GENtoRECOjetsCounter++;
      }
      else       theHistograms->fill("dR_unmatchedGENclosestRECOjet", "dR_unmatchedGENclosestRECOjet", 50, 0, 5, abs(physmath::deltaR(genJet, nearestRECOjet)), theWeight*LumiSF);

    }
    if (isreconstructed && selectedRECOjets.size() > 0)       theHistograms->fill("Pt_genJet_num", " Pt_genJet_num; GeV/c", 10, 0, 300, genJet.pt(), theWeight*LumiSF);
    theHistograms->fill("#GEN=>RECO", "#GEN=>RECO", 2, 0, 2, isreconstructed && selectedRECOjets.size() > 0, theWeight*LumiSF);
  }
  theHistograms->fill("2GENJetsToRECO", "2GENJetsToRECO", 2, 0, 2, GENtoRECOjetsCounter==2, theWeight*LumiSF);
  theHistograms->fill("AtLeast1GENJetToRECO", "AtLeast1GENJetToRECO", 2, 0, 2, GENtoRECOjetsCounter>0, theWeight*LumiSF);

  
  for (auto pair : nearestRECOjetstoGENjets)
  {
    ResolutionPlots(pair.first,pair.second,"SingleJetsReconstruction_",theWeight*LumiSF,"");
    theHistograms->fill("DeltaR_GENjet_vs_BestMatchedRECOJet", "DeltaR_GENjet_vs_BestMatchedRECOJet; #DeltaR", 20, 0, 0.5, abs(physmath::deltaR(pair.first, pair.second)), theWeight*LumiSF);
    theHistograms->fill("DeltaR_jets_vs_pt", "DeltaR jets vs pt;pt [GeV/c] ; #DeltaR", 10, 0, 300, 20, 0, 0.2, pair.first.pt(),abs(physmath::deltaR(pair.first, pair.second)), theWeight*LumiSF);
  }
  theHistograms->fill("2size_GENjetsRECONSTRUCTED", "2size_GENjetsRECONSTRUCTED", 10, -0.5, 9.5, RECOjetsfromGENjets.size(), theWeight*LumiSF);
}

void VZGAnalyzer::PlotJets(const phys::Particle &Jet0, const phys::Particle &Jet1, std::string prename, const float weight, std::string suffix)
{
       std::string name = "J0";
       theHistograms->fill(prename + name + "_charge_" + suffix, prename + name + "'s charge", 5, -2.5, 2.5, Jet0.charge(), weight);
       theHistograms->fill(prename + name + "_mass_" + suffix, prename + name + "'s mass", 63, 0, 252, Jet0.mass(), weight);
       theHistograms->fill(prename + name + "_trmass_" + suffix, prename + name + "'s trmass", 50, 0, 400, Jet0.p4().Mt(), weight);
       theHistograms->fill(prename + name + "_pt_" + suffix, prename + name + "'s p_{t}", 50, ptcut, 600, Jet0.pt(), weight);
       theHistograms->fill(prename + name + "_Y_" + suffix, prename + name + "'s Y", 50, -5, 5, Jet0.rapidity(), weight);
       theHistograms->fill(prename + name + "_eta_" + suffix, prename + name + "'s #eta", 50, -9, 9, Jet0.eta(), weight);
       theHistograms->fill(prename + name + "_phi_" + suffix, prename + name + "'s #phi", 50, -3.5, 3.5, Jet0.phi(), weight);

       name = "J1";
       theHistograms->fill(prename + name + "_charge_" + suffix, prename + name + "'s charge", 5, -2.5, 2.5, Jet1.charge(), weight);
       theHistograms->fill(prename + name + "_mass_" + suffix, prename + name + "'s mass", 63, 0, 252, Jet1.mass(), weight);
       theHistograms->fill(prename + name + "_trmass_" + suffix, prename + name + "'s trmass", 50, 0, 400, Jet1.p4().Mt(), weight);
       theHistograms->fill(prename + name + "_pt_" + suffix, prename + name + "'s p_{t}", 50, 30, 600, Jet1.pt(), weight);
       theHistograms->fill(prename + name + "_Y_" + suffix, prename + name + "'s Y", 50, -5, 5, Jet1.rapidity(), weight);
       theHistograms->fill(prename + name + "_eta_" + suffix, prename + name + "'s #eta", 50, -9, 9, Jet1.eta(), weight);
       theHistograms->fill(prename + name + "_phi_" + suffix, prename + name + "'s #phi", 50, -3.5, 3.5, Jet1.phi(), weight);

       name = "JJ";
       TLorentzVector JJp4 = Jet0.p4() + Jet1.p4();
       double JJdeltaEta = Jet0.eta() - Jet1.eta();
       double JJdeltaPhi = physmath::deltaPhi(Jet0.phi(), Jet1.phi());
       double JJdeltaR = abs(physmath::deltaR(Jet0, Jet1));

       theHistograms->fill(prename + name + "_mass_" + suffix, " Jets' mass", 10, 50, 120, JJp4.M(), weight);
       theHistograms->fill(prename + name + "_trmass_" + suffix, " Jets' trmass", 50, 0, 200, JJp4.Mt(), weight);
       theHistograms->fill(prename + name + "_pt_" + suffix, " Jets' p_{t}", 50, 0, 600, JJp4.Pt(), weight);
       theHistograms->fill(prename + name + "_deltaEta_" + suffix, " Jets' #Delta#eta", 50, -9, 9, JJdeltaEta, weight);
       theHistograms->fill(prename + name + "_deltaEtaabs_" + suffix, " Jets' |#Delta#eta|", 25, 0, 9, abs(JJdeltaEta), weight);
       theHistograms->fill(prename + name + "_deltaR_" + suffix, " Jets' #DeltaR", 25, -0.5, 9, JJdeltaR, weight);
       theHistograms->fill(prename + name + "_deltaPhi_" + suffix, " Jets' #Delta#phi", 50, -3.5, 3.5, JJdeltaPhi, weight);

       theHistograms->fill(prename + name + "_massvsdeltaEta_" + suffix, prename + name + "'s mass(x) vs #Delta#eta(y)", 12, 160, 1780, 10, -6.5, 6.5, JJp4.M(), JJdeltaEta, weight);
       theHistograms->fill(prename + name + "_massvsdeltaEtaabs_" + suffix, prename + name + "'s mass(x) vs |#Delta#eta|(y)", 12, 160, 1780, 10, -6.5, 6.5, JJp4.M(), abs(JJdeltaEta), weight);
}

void VZGAnalyzer::PlotJet(const phys::Particle &Jet, std::string prename, const float weight, std::string suffix)
{
       std::string name = " ";
       theHistograms->fill(prename + name + "_charge_" + suffix, prename + name + "'s charge", 5, -2.5, 2.5, Jet.charge(), weight);
       theHistograms->fill(prename + name + "_mass_" + suffix, prename + name + "'s mass", 63, 50, 120, Jet.mass(), weight);
       theHistograms->fill(prename + name + "_trmass_" + suffix, prename + name + "'s trmass", 50, 0, 200, Jet.p4().Mt(), weight);
       theHistograms->fill(prename + name + "_pt_" + suffix, prename + name + "'s p_{t}", 50, ptcut, 600, Jet.pt(), weight);
       theHistograms->fill(prename + name + "_Y_" + suffix, prename + name + "'s Y", 50, -5, 5, Jet.rapidity(), weight);
       theHistograms->fill(prename + name + "_eta_" + suffix, prename + name + "'s #eta", 50, -9, 9, Jet.eta(), weight);
       theHistograms->fill(prename + name + "_phi_" + suffix, prename + name + "'s #phi", 50, -3.5, 3.5, Jet.phi(), weight);
}


