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

bool IsARunForMVAFeat=true;
bool verbose = false;
bool signalSample = false;
bool fiducial_run =true;
double etacut=4.7;
double ptcut=30;
double LumiSF=3.32;//7.035;//--> 2016post||3.32-->2017||2.29-->2018||8.16--> 2016post||7.035--> 2016pre
//FOR DY PART1
//double LumiSF=7.42;
double dR_jetRatio_cut = 0.4;
double dR_FJRatio_cut = 0.8;
int cutsToApply=7;

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
  bool haveGoodFJCand=false;
  bool haveGoodDiJetCand=false;

  std::vector<phys::Particle> FJCand;

  foreach (const phys::Particle &fatJet, *genJetsAK8)
    {
      if (KinematicsOK(fatJet,ptcut,etacut)) // KinematicsOK(jet)
	FJCand.push_back(fatJet);
    }

  if (FJCand.size() > 0){
    std::stable_sort(FJCand.begin(), FJCand.end(), phys::Mass2Comparator(phys::ZMASS, phys::WMASS));
    haveGoodFJCand=(FJCand[0].mass()>50  &&   FJCand[0].mass()<120);
  }
  
  std::vector<phys::Particle> selectedGENjets;
  std::vector<phys::Boson<phys::Particle>> DiJetsCand;

  foreach (const phys::Particle &jet, *genJets)
  {
    if (KinematicsOK(jet,ptcut,etacut)) // KinematicsOK(jet)
      selectedGENjets.push_back(jet);
  }

  if (selectedGENjets.size() > 1){
    for (uint i = 0; i < selectedGENjets.size() - 1; i++) // Warning: size can be 0
      for (uint j = i+1; j < selectedGENjets.size(); j++)
        DiJetsCand.push_back(phys::Boson<phys::Particle>(selectedGENjets.at(i), selectedGENjets.at(j)));
  }

  if (DiJetsCand.size() > 0){
    std::stable_sort(DiJetsCand.begin(), DiJetsCand.end(), phys::Mass2Comparator(phys::ZMASS, phys::WMASS));
  //*V_JJCandidate = DiJetsCand.at(0);
    haveGoodDiJetCand=(DiJetsCand[0].mass()>50  &&   DiJetsCand[0].mass()<120);
  }
  //if(haveGoodFJCand) std::cout<<"FJ selection accomplished: FJ MASS = "<<FJCand[0].mass()<<std::endl;
  //if(haveGoodDiJetCand) std::cout<<"2AK4 selection accomplished: 2J MASS = "<<DiJetsCand[0].mass()<<std::endl;
  //if(haveGoodFJCand || haveGoodDiJetCand) std::cout<<"hadronic selection accomplished"<<std::endl;
  return (haveGoodFJCand || haveGoodDiJetCand);

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
    {
      std::vector<phys::Particle> selectedGENphotons;
      for (auto p : *genParticles)
	if (p.id() == 22 && KinematicsOK(p, 20, 2.4) && p.genStatusFlags().test(phys::isPrompt) &&  p.genStatusFlags().test(phys::fromHardProcess))
	  selectedGENphotons.push_back(p);
      
  
      if (selectedGENphotons.size()>=1)
	{
	  TLorentzVector  GEN_llPh ;
	  double mGEN_llPh;
	  double mGenZ;
	  if(genVBHelper_.ZtoChLep().size()>=1)
	    {
	      
	      //	      std::cout<< "entered ZToL"<<std::endl;
	      GEN_llPh = genVBHelper_.ZtoChLep()[0].daughter(0).p4()+genVBHelper_.ZtoChLep()[0].daughter(0).p4()+selectedGENphotons.at(0).p4();
	      //	      std::cout<< "GEN_llPh implemented"<<std::endl;
	      mGEN_llPh=GEN_llPh.M();
	      //	      std::cout<< "mass of GEN_llPh implemented"<<std::endl;
	      	      mGenZ=genVBHelper_.ZtoChLep()[0].mass();
	      //	      std::cout<< "GEN_Z mass implemented"<<std::endl;
	      
	      //theHistograms->fill("GEN mjj_vs_mjjG", "GEN mjj_vs_mjjG; mjj [GeV] ; mjj#gamma [GeV]", 35, 50, 120, 30, 50, 350, mGenV, mGEN_jjPh, theWeight*LumiSF);
	      theHistograms->fill("GEN mll_vs_mllG", "GEN mll_vs_mllG; mll [GeV] ; mll#gamma [GeV]", 30, 60, 120, 80, 50, 450, mGenZ, mGEN_llPh, theWeight*LumiSF);
	      
	      //	      std::cout<< "Histos filled"<<std::endl;
	      
	      //	      if(!fiducial_run || mGEN_llPh>95)std::cout<< "SigDef accomplished"<<std::endl;
	      
	      if(!fiducial_run) return true;
	      else if(mGEN_llPh>95) return true;//additional sig def fiducial requirement 
	    }
	}
    
    }
  //    return true;
  return false;
}
/*
bool VZGAnalyzer::IN_GENsignalDef_excludingFSR()
{
  if (LeptonicSignalConstraint() && HadronicSignalConstraint() && PhotonSignalConstraint())
    {
      
      return true;
    }
  return false;
}
*/

bool VZGAnalyzer::baselineRequirements()
{
  //----------------------------------------Building jj pairs ----------------------------------------//
  int topo=0;
  phys::Boson<phys::Jet> recoV;
  phys::Jet recoFJ;
  bool haveGoodRECODiJetCand=false;
  bool haveGoodRECOFJCand=false;
  //topo =  Reconstruct(&recoV,&recoFJ,&haveGoodRECODiJetCand,&haveGoodRECOFJCand);

  /*
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
  */
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
     
  if ( (haveGoodRECODiJetCand || haveGoodRECOFJCand) && selectedphotons.size() >0 && goodZ )
    return true;
  return false;
}


Bool_t VZGAnalyzer::cut(Int_t n, phys::Boson<phys::Jet> recoV, phys::Jet recoFJ, std::vector<phys::Photon> selectedPhotons, int VBTopo)
{ // returns false if the event has to be cut
  //  std::cout<<"entering cut"<<std::endl;
  if(selectedPhotons.size()<1) return false;

  TLorentzVector llGamma = Z->daughter(0).p4()+Z->daughter(1).p4()+selectedPhotons.at(0).p4();
  double mllGamma=llGamma.M();

  phys::Photon mostEnergeticPhoton;
  std::stable_sort(selectedPhotons.begin(), selectedPhotons.end(), phys::EComparator());

  mostEnergeticPhoton = selectedPhotons[0];
    
  std::pair<phys::Photon, phys::Jet> nearestRECOjetstoPhoton;

  if(fabs(physmath::deltaR(recoV.daughter(0), mostEnergeticPhoton))<fabs(physmath::deltaR(recoV.daughter(1), mostEnergeticPhoton)) )
    nearestRECOjetstoPhoton={mostEnergeticPhoton, recoV.daughter(0)};
  else if (fabs(physmath::deltaR(recoV.daughter(0), mostEnergeticPhoton))>fabs(physmath::deltaR(recoV.daughter(1), mostEnergeticPhoton)) )
    nearestRECOjetstoPhoton={mostEnergeticPhoton, recoV.daughter(1)};

  std::pair<phys::Photon, phys::Lepton> nearestLepToPhoton;

  if(fabs(physmath::deltaR(Z->daughter(0), mostEnergeticPhoton))<fabs(physmath::deltaR(Z->daughter(1), mostEnergeticPhoton)) )
    nearestLepToPhoton={mostEnergeticPhoton, Z->daughter(0)};
  else if (fabs(physmath::deltaR(Z->daughter(0), mostEnergeticPhoton))>fabs(physmath::deltaR(Z->daughter(1), mostEnergeticPhoton)) )
    nearestLepToPhoton={mostEnergeticPhoton, Z->daughter(1)};
  
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
  /*
  bool baseline = (VBTopo!=0 //(jets->size() > 1 || jetsAK8->size() > 0)
	&& ((recoV.mass() > 50 && recoV.mass() < 120)||(recoFJ.mass() > 50 && recoFJ.mass() < 120))
	&& KinematicsOK(recoV.daughter(0), ptcut,etacut)
	&& KinematicsOK(recoV.daughter(1), ptcut,etacut)
	&& Z->mass() > 60 && Z->mass() < 120
	&& KinematicsOK(Z->daughter(0), 5.,etacut)
	&& KinematicsOK(Z->daughter(1), 5.,etacut)
	&& fabs(physmath::deltaR(nearestRECOjetstoPhoton.first, nearestRECOjetstoPhoton.second))> dR_jetRatio_cut
	&& selectedPhotons.size()>0
	&& (!fiducial_run || mllGamma>95));
  */
  bool baseline = (VBTopo!=0 //(jets->size() > 1 || jetsAK8->size() > 0)
		   && Z->mass() > 60 && Z->mass() < 120
		   && KinematicsOK(Z->daughter(0), 5.,etacut)
		   && KinematicsOK(Z->daughter(1), 5.,etacut)
		   && selectedPhotons.size()>0
		   && (!fiducial_run || mllGamma>95));
  
  switch (n)
  {
  case 1://baseline (with dRjG) 
    if (baseline)
      return true;
    break;

  case 2://high mllG
    if (baseline
	&& mllGamma>150)
      return true;
    break;

  case 3://dRlG>0.5
    if (baseline
	&& fabs(physmath::deltaR(nearestLepToPhoton.first, nearestLepToPhoton.second))> 0.5
	&& mllGamma>150)
      return true;
    break;

  case 4://mjj=mV+-15
    if (baseline
	&& ((VBTopo == 1 && recoV.mass() > 65 && recoV.mass() < 115)||(VBTopo == -1 && recoFJ.mass() > 65 && recoFJ.mass() < 115))
	&& fabs(physmath::deltaR(nearestLepToPhoton.first, nearestLepToPhoton.second))> 0.5
	&& mllGamma>150)
      return true;
    break;

  case 5://ptj0>40
    if (baseline
	&& (VBTopo==-1 || KinematicsOK(recoV.daughter(0), 40, etacut))
	&& ((VBTopo == 1 && recoV.mass() > 65 && recoV.mass() < 115)||(VBTopo == -1 && recoFJ.mass() > 65 && recoFJ.mass() < 115))
	&& fabs(physmath::deltaR(nearestLepToPhoton.first, nearestLepToPhoton.second))> 0.5
	&& mllGamma>150)
      return true;
    break;

  case 6://pt V > 90
    if (baseline
	&& ((VBTopo==1 && recoV.pt()>90) || (VBTopo==-1 && recoFJ.pt()>90))
	//&& SumFWM(0, 't', lljjG) >2.10 //CT
	&& Z->mass() > 75 && Z->mass() < 110
	&& (VBTopo==-1 || KinematicsOK(recoV.daughter(0), 40, etacut))
	&& ((VBTopo == 1 && recoV.mass() > 65 && recoV.mass() < 115)||(VBTopo == -1 && recoFJ.mass() > 65 && recoFJ.mass() < 115))
	&& fabs(physmath::deltaR(nearestLepToPhoton.first, nearestLepToPhoton.second))> 0.1
	&& mllGamma>150)
      return true;
    break;

  case 7://mll=mZ+-15/20
    if (baseline
        && Z->mass() > 75 && Z->mass() < 110
	&& (VBTopo==-1 || KinematicsOK(recoV.daughter(0), 40, etacut))
	&& ((VBTopo == 1 && recoV.mass() > 65 && recoV.mass() < 115)||(VBTopo == -1 && recoFJ.mass() > 65 && recoFJ.mass() < 115))
	&& fabs(physmath::deltaR(nearestLepToPhoton.first, nearestLepToPhoton.second))> 0.5
	&& mllGamma>150)
      return true;
    break;

  case 8://ptj1>40
    if (baseline
	&& (VBTopo==-1 || KinematicsOK(recoV.daughter(1), 40, etacut) )
	&& ((VBTopo==1 && recoV.pt()>90) || (VBTopo==-1 && recoFJ.pt()>90))
	&& Z->mass() > 75 && Z->mass() < 110
	&& (VBTopo==-1 || KinematicsOK(recoV.daughter(0), 40, etacut))
	&& ((VBTopo == 1 && recoV.mass() > 65 && recoV.mass() < 115)||(VBTopo == -1 && recoFJ.mass() > 65 && recoFJ.mass() < 115))
	&& mllGamma>150)
      return true;
    break;

  case 9://ptj0 > 70
    if (baseline
	&& (VBTopo==-1 || (KinematicsOK(recoV.daughter(0), 70, etacut) && KinematicsOK(recoV.daughter(1), 40, etacut)))
	&& ((VBTopo==1 && recoV.pt()>90) || (VBTopo==-1 && recoFJ.pt()>90))
	//&& SumFWM(0, 't', lljjG) >2.10 //CT
	&& Z->mass() > 75 && Z->mass() < 110
	&& ((VBTopo == 1 && recoV.mass() > 65 && recoV.mass() < 115)||(VBTopo == -1 && recoFJ.mass() > 65 && recoFJ.mass() < 115))
	&& mllGamma>150)
      return true;
    break;

  case 10://dRLG>1.2
    if (baseline
	&& fabs(physmath::deltaR(nearestLepToPhoton.first, nearestLepToPhoton.second))> 1.2
	&& (VBTopo==-1 || (KinematicsOK(recoV.daughter(0), 70, etacut) && KinematicsOK(recoV.daughter(1), 40, etacut)))
	&& ((VBTopo==1 && recoV.pt()>90) || (VBTopo==-1 && recoFJ.pt()>90))
	//	&& SumFWM(0, 't', lljjG) >2.10 //CT
	&& Z->mass() > 75 && Z->mass() < 110
	&& ((VBTopo == 1 && recoV.mass() > 65 && recoV.mass() < 115)||(VBTopo == -1 && recoFJ.mass() > 65 && recoFJ.mass() < 115))
	&& mllGamma>150)
      return true;
    break;

  case 11://dRjj< 2.4
    if (baseline
	&& fabs(physmath::deltaR(Z->daughter(0), Z->daughter(1)))< 2.4
	&& (VBTopo==-1 || (KinematicsOK(recoV.daughter(0), 70, etacut) && KinematicsOK(recoV.daughter(1), 40, etacut)))
	&& fabs(physmath::deltaR(nearestLepToPhoton.first, nearestLepToPhoton.second))> 1.2
	&& ((VBTopo==1 && recoV.pt()>90) || (VBTopo==-1 && recoFJ.pt()>90))
	//&& SumFWM(0, 't', lljjG) >2.10 //CT
	&& Z->mass() > 75 && Z->mass() < 110
	&& ((VBTopo == 1 && recoV.mass() > 65 && recoV.mass() < 115)||(VBTopo == -1 && recoFJ.mass() > 65 && recoFJ.mass() < 115))
	&& mllGamma>150)
      return true;
    break;

  case 12://baseline+mjj=mV+-15 & dRjj 2.4 & FWM T0 
    if (baseline
	&& (VBTopo==-1 || fabs(physmath::deltaR(recoV.daughter(0),recoV.daughter(1)))< 2.4)
	&& fabs(physmath::deltaR(Z->daughter(0), Z->daughter(1)))< 2.4
	&& (VBTopo==-1 || (KinematicsOK(recoV.daughter(0), 70, etacut) && KinematicsOK(recoV.daughter(1), 40, etacut)))
	&& fabs(physmath::deltaR(nearestLepToPhoton.first, nearestLepToPhoton.second))> 1.2
	&& ((VBTopo==1 && recoV.pt()>90) || (VBTopo==-1 && recoFJ.pt()>90))
	//&& SumFWM(0, 't', lljjG) >2.10 //CT
	&& Z->mass() > 75 && Z->mass() < 110
	&& ((VBTopo == 1 && recoV.mass() > 65 && recoV.mass() < 115)||(VBTopo == -1 && recoFJ.mass() > 65 && recoFJ.mass() < 115))
	&& mllGamma>150)
      return true;
    break;

  case 13://baseline+mjj=mV+-15 & dRjj 2.4 & FWM T0 
    if (baseline
	&& (VBTopo==-1 || fabs(physmath::deltaR(recoV.daughter(0),recoV.daughter(1)))< 2.4)
	&& fabs(physmath::deltaR(Z->daughter(0), Z->daughter(1)))< 2.4
	&& (VBTopo==-1 || (KinematicsOK(recoV.daughter(0), 70, etacut) && KinematicsOK(recoV.daughter(1), 40, etacut)))
	&& fabs(physmath::deltaR(nearestLepToPhoton.first, nearestLepToPhoton.second))> 1.2
	&& ((VBTopo==1 && recoV.pt()>90) || (VBTopo==-1 && recoFJ.pt()>90))
	//&& SumFWM(0, 't', lljjG) >2.10 //CT
	&& Z->mass() > 75 && Z->mass() < 110
	&& ((VBTopo == 1 && recoV.mass() > 65 && recoV.mass() < 115)||(VBTopo == -1 && recoFJ.mass() > 65 && recoFJ.mass() < 115))
	&& mllGamma>150)
      return true;
    break;

  default:
    return true;
  }
  return false;

  /*
  switch (n)
  {
  case 13:
  case 12://baseline+mjj=mV+-15 & dRjj 2.4 & FWM T0 
    if (!(VBTopo==-1 && fabs(physmath::deltaR(recoV.daughter(0),recoV.daughter(1)))< 2.4))
      return false;
  case 11:
    if(!fabs(physmath::deltaR(Z->daughter(0), Z->daughter(1)))< 2.4)
      return false;
  case 10:
    if(!fabs(physmath::deltaR(nearestLepToPhoton.first, nearestLepToPhoton.second))> 1.2)
      return false;
  case 9:
    if(!(VBTopo==-1 || KinematicsOK(recoV.daughter(0), 70, etacut))   )
      return false;
  case 8:
    if(!(VBTopo==-1 || KinematicsOK(recoV.daughter(1), 40, etacut)) )
      return false;
  case 7:
    if(!( (VBTopo==-1 && recoV.pt()>90)  || (VBTopo==-1 && recoFJ.pt()>90) ) )
      return false;
  case 6:
    if(!(Z->mass() > 75 && Z->mass() < 110) )
      return false;
  case 5:
    if(!(VBTopo==-1 || KinematicsOK(recoV.daughter(0), 40, etacut)) )
      return false;
  case 4:
    if(!((VBTopo == 1 && recoV.mass() > 65 && recoV.mass() < 115)||(VBTopo == -1 && recoFJ.mass() > 65 && recoFJ.mass() < 115)) )
      return false;
  case 3:
    if(!fabs(physmath::deltaR(nearestLepToPhoton.first, nearestLepToPhoton.second))> 0.5)
      return false;
  case 2://high mllG
    if (!mllGamma>150)    return false;
  case 1://baseline (with dRjG) 
    if (!baseline)    return false;

  default:
    return true;

  }
  */  



}

void VZGAnalyzer::analyze()
{ // It's the only member function running each event.
  if(IsARunForMVAFeat) return;
  if(verbose==true)
    {
      cout << "----------------------------------------------------------------" << endl;
      cout << "Run: " << run << " event: " << event << endl;
    }

  //  genAnalyze();

  //if (IN_GENsignalDef())     QuarksToJets();

  //  if (IN_GENsignalDef()) std::cout<<"passing sig def"<<std::endl;
  int VBTopo = 0;
  phys::Boson<phys::Jet> recoV;
  phys::Jet recoFJ;
  bool haveGoodRECODiJetCand=false;
  bool haveGoodRECOFJCand=false;

  
  std::vector<phys::Photon> selectedphotons;
  PhotonSelection(&selectedphotons);
  //  std::cout<<"selected photons size "<<selectedphotons.size()<<std::endl;
  if(selectedphotons.size()<1){
    if (IN_GENsignalDef())
      {
	printHistos(0, "sign", recoV, recoFJ, selectedphotons, VBTopo);
      } // not signalConstraint anymore


    else{    printHistos(0, "bckg", recoV, recoFJ, selectedphotons, VBTopo);}

    return;
  }
  //  std::cout<<"PASSING selected photons size "<<selectedphotons.size()<<std::endl;

  phys::Photon mostEnergeticPhoton;
  //    if (selectedphotons.size() > 0)
  //  {
  std::stable_sort(selectedphotons.begin(), selectedphotons.end(), phys::EComparator());
  mostEnergeticPhoton = selectedphotons[0];

  
  VBTopo=Reconstruct(&recoV,&recoFJ,&haveGoodRECODiJetCand,&haveGoodRECOFJCand,&mostEnergeticPhoton);

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
  //   theHistograms->fill("DeltaR_mostEnergeticPhoton_vs_BestMatchedRECOJet", "DeltaR_mostEnergeticPhoton_vs_BestMatchedRECOJet; #DeltaR", 50, 0, 5, fabs(physmath::deltaR(pair.first, pair.second)), theWeight*LumiSF);
  //   theHistograms->fill("DeltaR_vs_Deltapt", "DeltaR_vs_Deltapt;#Delta pt [GeV/c] ; #DeltaR", 20, -100, 100, 20, 0, 2, pair.first.pt()-pair.second.pt(),fabs(physmath::deltaR(pair.first, pair.second)), theWeight*LumiSF);
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
    theHistograms->fill("DeltaR_mostEnergeticPhoton_vs_BestMatchedRECOJet", "DeltaR_mostEnergeticPhoton_vs_BestMatchedRECOJet; #DeltaR", 50, 0, 5, fabs(physmath::deltaR(pair.first, pair.second)), theWeight*LumiSF);
    theHistograms->fill("DeltaR_vs_Deltapt", "DeltaR_vs_Deltapt;#Delta pt [GeV/c] ; #DeltaR", 20, -100, 100, 20, 0, 2, pair.first.pt()-pair.second.pt(),fabs(physmath::deltaR(pair.first, pair.second)), theWeight*LumiSF);
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

    printHistos(0, "sign", recoV, recoFJ, selectedphotons, VBTopo);
  } // not signalConstraint anymore

  //---------------BACKGROUND-EVENTS---------------//

  else    printHistos(0, "bckg", recoV, recoFJ, selectedphotons, VBTopo);

}

void VZGAnalyzer::fillFeatTree(FeatList &list, bool &passingPresel )
{
  passingPresel = false;
  if(!IsARunForMVAFeat)  return;
  if(!IN_GENsignalDef()) return;
  //std::cout<<"0: entering fillFeatTree "<<std::endl;


  int nbOfCutsPassed=0;
  
  int VBTopo = 0;
  phys::Boson<phys::Jet> recoV;
  phys::Jet recoFJ;
  bool haveGoodRECODiJetCand=false;
  bool haveGoodRECOFJCand=false;

  
  std::vector<phys::Photon> selectedphotons;
  PhotonSelection(&selectedphotons);

  //std::cout<<"1: passing photon selection "<<std::endl;
  //  std::cout<<"---------------------------------------------------------"<<std::endl;//
  if(selectedphotons.size()<1) {
    list.f_nbOfCutsPassed = 0;  
    return;
  }
  
  phys::Photon mostEnergeticPhoton;

  std::stable_sort(selectedphotons.begin(), selectedphotons.end(), phys::EComparator());
  mostEnergeticPhoton = selectedphotons[0];

  
  VBTopo=Reconstruct(&recoV,&recoFJ,&haveGoodRECODiJetCand,&haveGoodRECOFJCand,&mostEnergeticPhoton);

  //while(cut(nbOfCutsPassed, recoV, recoFJ, selectedphotons, VBTopo))    nbOfCutsPassed++;
  
  //  list.f_nbOfCutsPassed = nbOfCutsPassed;  

  if(VBTopo!=1) return; 
  //  std::cout<<"2: VBTopo "<<VBTopo<<std::endl;

  if (!cut(2, recoV, recoFJ, selectedphotons, VBTopo)) return;

  //  std::cout<<"3: passing cuts "<<std::endl;

  TLorentzVector llPh = Z->daughter(0).p4()+Z->daughter(1).p4()+selectedphotons.at(0).p4();
  TLorentzVector l0Ph = Z->daughter(0).p4()+selectedphotons.at(0).p4();
  TLorentzVector l1Ph = Z->daughter(1).p4()+selectedphotons.at(0).p4();

  double mllPh=llPh.M();
  double m2llPh=llPh.M2();
  double m2l0Ph=l0Ph.M2();
  double m2l1Ph=l1Ph.M2();
  double mll=Z->mass();

  //  std::cout<<"4: first vars filled "<<std::endl;

  
  double dPhiL0G, dPhiL1G, dPhiLL, recoVMass, ptl0, ptl1, FWMT0, ptGamma, ptJ0, ptJ1, etaJ0, etaJ1, etaL0, etaL1, FWMT1, FWMT2, FWMT3, FWMT4, FWMT5, FWMT6, dPhiJ0G, dPhiJ1G, dPhiJJ, dPhiL0J0, dPhiL1J1, dPhiL0J1, dPhiL1J0;  

  dPhiLL=fabs(physmath::deltaPhi(Z->daughter(0).phi(), Z->daughter(1).phi()));
  dPhiL0G=fabs(physmath::deltaPhi(Z->daughter(0).phi(), selectedphotons.at(0).phi() ));
  dPhiL1G=fabs(physmath::deltaPhi(Z->daughter(1).phi(), selectedphotons.at(0).phi()));

  dPhiJJ=fabs(physmath::deltaPhi(recoV.daughter(0).phi(), recoV.daughter(1).phi()));
  dPhiJ0G=fabs(physmath::deltaPhi(recoV.daughter(0).phi(), selectedphotons.at(0).phi() ));
  dPhiJ1G=fabs(physmath::deltaPhi(recoV.daughter(1).phi(), selectedphotons.at(0).phi()));

  dPhiL0J0=fabs(physmath::deltaPhi(Z->daughter(0).phi(), recoV.daughter(0).phi() ));
  dPhiL1J0=fabs(physmath::deltaPhi(Z->daughter(1).phi(), recoV.daughter(0).phi() ));
  dPhiL0J1=fabs(physmath::deltaPhi(Z->daughter(0).phi(), recoV.daughter(1).phi() ));
  dPhiL1J1=fabs(physmath::deltaPhi(Z->daughter(1).phi(), recoV.daughter(1).phi() ));


  std::vector<TLorentzVector> lljjG;
  lljjG.push_back(mostEnergeticPhoton.p4());
  lljjG.push_back(Z->daughter(0).p4());
  lljjG.push_back(Z->daughter(1).p4());
  if (VBTopo==1){
    lljjG.push_back(recoV.daughter(0).p4());
    lljjG.push_back(recoV.daughter(1).p4());
    recoVMass=recoV.mass();
  }
  else if (VBTopo==-1){
    lljjG.push_back(recoFJ.p4());
    recoVMass=recoFJ.mass();
  }
  std::vector<TLorentzVector> llG;
  llG.push_back(Z->daughter(0).p4());
  llG.push_back(Z->daughter(1).p4());
  llG.push_back(selectedphotons.at(0).p4());


    
  std::pair<phys::Photon, phys::Lepton> nearestChLeptToPhoton;

  if(fabs(physmath::deltaR(Z->daughter(0), mostEnergeticPhoton))<fabs(physmath::deltaR(Z->daughter(1), mostEnergeticPhoton)) )
    nearestChLeptToPhoton={mostEnergeticPhoton, Z->daughter(0)};
  else if (fabs(physmath::deltaR(Z->daughter(0), mostEnergeticPhoton))>fabs(physmath::deltaR(Z->daughter(1), mostEnergeticPhoton)) )
    nearestChLeptToPhoton={mostEnergeticPhoton, Z->daughter(1)};

  //  deltaR_LGamma=fabs(physmath::deltaR(Z->daughter(0), selectedphotons.at(0)));

  //std::cout<<"dRLPH "<<deltaR_LGamma<<std::endl;

  //  if (deltaR_LGamma>10) deltaR_LGamma=11;

  
  FWMT0 =SumFWM(0, 't', lljjG);
  FWMT1 =SumFWM(1, 't', lljjG);
  FWMT2 =SumFWM(2, 't', lljjG);
  FWMT3 =SumFWM(3, 't', lljjG);
  FWMT4 =SumFWM(4, 't', lljjG);
  FWMT5 =SumFWM(5, 't', lljjG);
  FWMT6 =SumFWM(6, 't', lljjG);
  //  std::cout<<"----------------------------------------"<<std::endl;    
  list.f_weight = theWeight;
  //  std::cout<<"weight: "<<theWeight<<std::endl;//  theHistograms->fill("the Weight", "the Weight", 50, -5, 5, theWeight, 1);

  list.f_mll  = Z->mass();
  //  std::cout<<"ZCandMass: "<<Z->mass()<<std::endl;//    theHistograms->fill("ZCand mass", "ZCand mass", 30, 60, 120, Z->mass(), 1);

  list.f_ptl1 = Z->daughter(0).pt();
  //  std::cout<<"ptl1: "<<Z->daughter(0).pt()<<std::endl;//      theHistograms->fill("ptl1", "ptl1", 50,  20, 120, Z->daughter(0).pt(), 1);

  list.f_ptl2 = Z->daughter(1).pt();
  //  std::cout<<"ptl2: "<<Z->daughter(1).pt()<<std::endl;//  theHistograms->fill("ptl2", "ptl2", 50,  20, 120, Z->daughter(1).pt(), 1);

  list.f_ptGamma = selectedphotons.at(0).pt();
  //  std::cout<<"ptGamma: "<<selectedphotons.at(0).pt()<<std::endl;//  theHistograms->fill("ptGamma", "ptGamma", 50,  20, 120, selectedphotons.at(0).pt(), 1);
  
  list.f_ptJ0 = recoV.daughter(0).pt();
  //  std::cout<<"ptJ0: "<<recoV.daughter(0).pt()<<std::endl;//  theHistograms->fill("ptJ0", "ptJ0", 50,  30, 130, recoV.daughter(0).pt(), 1);

  list.f_ptJ1 = recoV.daughter(1).pt();
  //  std::cout<<"ptJ1: "<<recoV.daughter(1).pt()<<std::endl;//  theHistograms->fill("ptJ1", "ptJ1", 50,  30, 130, recoV.daughter(1).pt(), 1);

  list.f_etaJ0 = recoV.daughter(0).eta();
  list.f_etaJ1 = recoV.daughter(1).eta();
  list.f_etaL0 = Z->daughter(0).eta();
  list.f_etaL1 = Z->daughter(1).eta();
  
  /*
  std::cout<<"etaJ0: "<<recoV.daughter(0).eta()<<std::endl;//  theHistograms->fill("etaJ0", "etaJ0", 50,  -4, 4, recoV.daughter(0).eta(), 1);
  std::cout<<"etaJ1: "<<recoV.daughter(1).eta()<<std::endl;//  theHistograms->fill("etaJ1", "etaJ1", 50,  -4, 4, recoV.daughter(1).eta(), 1);
  std::cout<<"etaL0: "<<Z->daughter(0).eta()<<std::endl;//  theHistograms->fill("etaL0", "etaL0", 25,  -2.5, 2.5, Z->daughter(0).eta(), 1);
  std::cout<<"etaL1: "<<Z->daughter(1).eta()<<std::endl;//  theHistograms->fill("etaL1", "etaL1", 25,  -2.5, 2.5, Z->daughter(1).eta(), 1);
  */
  
  list.f_dPhiL0G = dPhiL0G;
  list.f_dPhiL1G = dPhiL1G;
  list.f_dPhiLL = dPhiLL;
  /*
  theHistograms->fill("dPhiL0G", "dPhiL0G", 150,  0, 3, dPhiL0G, 1);
  theHistograms->fill("dPhiL1G", "dPhiL1G", 150,  0, 3, dPhiL1G, 1);
  theHistograms->fill("dPhiLL", "dPhiLL", 150,  0, 3, dPhiLL, 1);
  */
  
  list.f_dPhiJ0G = dPhiJ0G;
  list.f_dPhiJ1G = dPhiJ1G;
  list.f_dPhiJJ = dPhiJJ;
  /*
  theHistograms->fill("dPhiJ0G", "dPhiJ0G", 150,  0, 3, dPhiL0G, 1);
  theHistograms->fill("dPhiJ1G", "dPhiJ1G", 150,  0, 3, dPhiL1G, 1);
  theHistograms->fill("dPhiJJ", "dPhiJJ", 150,  0, 3, dPhiLL, 1);
  */
  list.f_dPhiL0J0 = dPhiL0J0;
  list.f_dPhiL0J1 = dPhiL0J1;
  list.f_dPhiL1J0 = dPhiL1J0;
  list.f_dPhiL1J1 = dPhiL1J1;
  /*
  theHistograms->fill("dPhiL0J0", "dPhiL0J0", 150,  0, 3, dPhiL0J0, 1);
  theHistograms->fill("dPhiL0J1", "dPhiL0J1", 150,  0, 3, dPhiL0J1, 1);
  theHistograms->fill("dPhiL1J0", "dPhiL1J0", 150,  0, 3, dPhiL1J0, 1);
  theHistograms->fill("dPhiL1J1", "dPhiL1J1", 150,  0, 3, dPhiL1J1, 1);
  */
  //  list.f_deltaR_LGamma = deltaR_LGamma;
  list.f_recoVMass=recoVMass;
  theHistograms->fill("recoVMass", "recoVMass", 35, 50, 120, recoVMass, 1);

  list.f_FWMT0=FWMT0;
  list.f_FWMT1=FWMT1;
  list.f_FWMT2=FWMT2;
  list.f_FWMT3=FWMT3;
  list.f_FWMT4=FWMT4;
  /*
  theHistograms->fill("FWMT0", "FWMT0", 25, 0, 5, FWMT0, 1);
  theHistograms->fill("FWMT1", "FWMT1", 25, 0, 5, FWMT1, 1);
  theHistograms->fill("FWMT2", "FWMT2", 50, 0, 1, FWMT2, 1);
  theHistograms->fill("FWMT3", "FWMT3", 25, 0, 5, FWMT3, 1);
  theHistograms->fill("FWMT4", "FWMT4", 50, 0, 0.25, FWMT4, 1);
  */
  //  list.f_FWMT5=FWMT5;
  //  list.f_FWMT6=FWMT6;

  //  std::cout<<"5: full list filled "<<std::endl;
  //featureTree.Fill();
  list.f_nbOfCutsPassed = nbOfCutsPassed;  
  passingPresel=true;
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


int VZGAnalyzer::Reconstruct(phys::Boson<phys::Jet> *V_JJCandidate, phys::Jet *V_FJCandidate, bool *haveGoodRECODiJetCand, bool *haveGoodRECOFJCand, phys::Photon *gamma)
{
  int hadrTopo = 0;
  double peakDist_mFJCand=40.;
  double peakDist_mDJCand=40.;
  
  /*
  bool haveGoodRECOFJCand=false;
  bool haveGoodRECODiJetCand=false;
  */
  std::vector<phys::Jet> FJCand;

  foreach (const phys::Jet &fatJet, *jetsAK8)
    {
      if (KinematicsOK(fatJet,ptcut,etacut) && fabs(physmath::deltaR(fatJet,*gamma))> dR_FJRatio_cut) // KinematicsOK(jet)
	FJCand.push_back(fatJet);
    }

  if (FJCand.size() > 0){
    std::stable_sort(FJCand.begin(), FJCand.end(), phys::Mass2Comparator(phys::ZMASS, phys::WMASS));
    *haveGoodRECOFJCand = (FJCand[0].mass()>50  &&   FJCand[0].mass()<120);
  }
  if(*haveGoodRECOFJCand){
    *V_FJCandidate = FJCand.at(0);
    peakDist_mFJCand=fabs(V_FJCandidate->mass()-phys::ZMASS);
    if(fabs(V_FJCandidate->mass()-phys::WMASS)<peakDist_mFJCand)    peakDist_mFJCand=fabs(V_FJCandidate->mass()-phys::WMASS);
  }

  
  std::vector<phys::Jet> selectedJets;
  std::vector<phys::Boson<phys::Jet>> DiJetsCand;

  foreach (const phys::Jet &jet, *jets)
  {
    if (KinematicsOK(jet,ptcut,etacut) && fabs(physmath::deltaR(jet,*gamma))> dR_jetRatio_cut) // KinematicsOK(jet)
      selectedJets.push_back(jet);
  }

  if (selectedJets.size() > 1){
    for (uint i = 0; i < selectedJets.size() - 1; i++) // Warning: size can be 0
      for (uint j = i+1; j < selectedJets.size(); j++)
        DiJetsCand.push_back(phys::Boson<phys::Jet>(selectedJets.at(i), selectedJets.at(j)));
  }

  if (DiJetsCand.size() > 0){
    std::stable_sort(DiJetsCand.begin(), DiJetsCand.end(), phys::Mass2Comparator(phys::ZMASS, phys::WMASS));
  //*V_JJCandidate = DiJetsCand.at(0);
    *haveGoodRECODiJetCand=(DiJetsCand[0].mass()>50  &&   DiJetsCand[0].mass()<120);
  }
  if(*haveGoodRECODiJetCand){
    *V_JJCandidate = DiJetsCand.at(0);
    peakDist_mDJCand=fabs(V_JJCandidate->mass()-phys::ZMASS);
    if(fabs(V_JJCandidate->mass()-phys::WMASS)<peakDist_mDJCand)    peakDist_mDJCand=fabs(V_JJCandidate->mass()-phys::WMASS);
  }

  if(!(*haveGoodRECOFJCand || *haveGoodRECODiJetCand) ){
    hadrTopo=0;
  }else{
    if(*haveGoodRECODiJetCand) hadrTopo=1;//    if(peakDist_mDJCand<peakDist_mFJCand) hadrTopo=1;
    else hadrTopo=0;//    else hadrTopo=-1;
  }
  //_________________________________
  /*
  if(*haveGoodRECODiJetCand){
    theHistograms->fill(" HAD RECO 2J cand mass", " HAD RECO 2J cand mass ; mjj Cand. [GeV]", 28, 50, 120, V_JJCandidate->mass());
    theHistograms->fill(" HAD RECO 2J/FJ cand mass", " HAD RECO 2J/FJ cand mass ; mVB Cand. [GeV]", 28, 50, 120, V_JJCandidate->mass());
  }else if(*haveGoodRECOFJCand){
    theHistograms->fill(" EXTRA HAD RECO FJ cand mass", " EXTRA HAD RECO 2J cand mass ; mFJ Cand. [GeV]", 28, 50, 120, V_FJCandidate->mass());
    theHistograms->fill(" HAD RECO 2J/FJ cand mass", " HAD RECO 2J/FJ cand mass ; mVB Cand. [GeV]", 28, 50, 120, V_FJCandidate->mass());    
  }
    
  if(*haveGoodRECOFJCand)
    theHistograms->fill(" HAD RECO FJ cand mass", " HAD RECO FJ cand mass ; mFJ Cand. [GeV]", 28, 50, 120, V_FJCandidate->mass());
  if(*haveGoodRECOFJCand && *haveGoodRECODiJetCand)
    theHistograms->fill(" HAD RECO 2J vs FJ", " HAD RECO HAD RECO 2J vs FJ ; m2J Cand. [GeV]; mFJ Cand. [GeV]", 28, 50, 120, 28, 50, 120, V_JJCandidate->mass(), V_FJCandidate->mass());

  //_________________________________
  */
  
  //  if(haveGoodRECODiJetCand || haveGoodRECOFJCand) std::cout<<"jets reconstruction worked"<<std::endl;
  // return (haveGoodRECOFJCand || haveGoodRECODiJetCand);


  /*
  
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
  */
  if(IsARunForMVAFeat) hadrTopo=*haveGoodRECODiJetCand;
  
  return hadrTopo;
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
  //std::cout << "Number of selected RECO photons = " << phot->size() << std::endl;
}


void VZGAnalyzer::CompatibilityTest(phys::Boson<phys::Jet> bestCandidate, phys::Boson<phys::Particle> genVB, std::string sample, std::string algorithm)
{
}

void VZGAnalyzer::printHistos(uint i, std::string histoType, phys::Boson<phys::Jet> recoV, phys::Jet recoFJ, std::vector<phys::Photon> selectedphotons, int VBTopo)
{

  std::vector<std::string> cuts = {"0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13"};
  std::vector<std::string> orders = {"0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13"};

  std::string hadTopo="no_jets";
  if(VBTopo==1)
    hadTopo="diJet_Cand";
  else if (VBTopo==-1)
    hadTopo="FJCand";

  if(selectedphotons.size()<1 || VBTopo==0){
    if(i==0)     theHistograms->fill("#AAA_cut_flow_" + histoType, "Cut flow", cutsToApply, 0, cutsToApply, i, (theWeight*LumiSF));//theHistograms->fill("#AAA_cut_flow_" + histoType, "Cut flow", cuts.size(), 0, cuts.size(), i, (theWeight*LumiSF));
    return;
  }
  TLorentzVector jjPh = recoV.daughter(0).p4()+recoV.daughter(1).p4()+selectedphotons.at(0).p4();
  double mjjPh=jjPh.M();
  double m2jjPh=jjPh.M2();

  TLorentzVector llPh = Z->daughter(0).p4()+Z->daughter(1).p4()+selectedphotons.at(0).p4();
  TLorentzVector l0Ph = Z->daughter(0).p4()+selectedphotons.at(0).p4();
  TLorentzVector l1Ph = Z->daughter(1).p4()+selectedphotons.at(0).p4();

  double mllPh=llPh.M();
  double m2llPh=llPh.M2();
  double m2l0Ph=l0Ph.M2();
  double m2l1Ph=l1Ph.M2();
  double mll=Z->mass();
  double mjj=recoV.mass();


  std::vector<TLorentzVector> lljj;
  lljj.push_back(Z->daughter(0).p4());
  lljj.push_back(Z->daughter(1).p4());
  lljj.push_back(recoV.daughter(0).p4());
  lljj.push_back(recoV.daughter(1).p4());

  
  if (i <= cutsToApply && cut(i, recoV, recoFJ, selectedphotons, VBTopo))
  {

    theHistograms->fill("#AAA_cut_flow_" + histoType, "Cut flow", cutsToApply, 0, cutsToApply, i, (theWeight*LumiSF));
    
    if(VBTopo==1){
      theHistograms->fill("V_vs_Z_pt_" + histoType + cuts.at(i), "V_vs_Z_pt_" + histoType + cuts.at(i) +";V p_{t} [GeV/c]; #Z p_{t} [GeV/c]", 30, 0, 300, 30, 0, 300, recoV.pt(), Z->pt(), theWeight*LumiSF);
      //      printHistos(1, "sign", recoV, recoFJ, selectedphotons, VBTopo);
      theHistograms->fill("recoVMass_" + histoType + cuts.at(i), "mass of recoV", 40, 0, 200, recoV.mass(), (theWeight*LumiSF));
      theHistograms->fill("recoVDaughter0Pt_" + histoType + cuts.at(i), "pt of recoVDaughter0", 50, 0, 600, recoV.daughter(0).pt(), (theWeight*LumiSF));
      theHistograms->fill("recoVDaughter1Pt_" + histoType + cuts.at(i), "pt of recoVDaughter1", 50, 0, 600, recoV.daughter(1).pt(), (theWeight*LumiSF));
      /*
      std::cout<<"-----------------------------------------------------------------" <<endl;
      std::cout<<"probb ="<< recoV.daughter(0).deepFlavour().probb <<endl;
      std::cout<<"probc ="<< recoV.daughter(0).deepFlavour().probc <<endl;
      std::cout<<"probg ="<< recoV.daughter(0).deepFlavour().probg <<endl;
      std::cout<<"problepb ="<< recoV.daughter(0).deepFlavour().problepb <<endl;
      std::cout<<"probbb ="<< recoV.daughter(0).deepFlavour().probbb <<endl;
      std::cout<<"probuds ="<< recoV.daughter(0).deepFlavour().probuds <<endl;
      */

      /*
      theHistograms->fill("j0_PN_ZvsQCD"+histoType + cuts.at(i), "j0_PN_ZvsQCD"+histoType + cuts.at(i)+"; lead. jet ParticleNet ZvsQCD", 20, 0, 1, recoV.daughter(0).particleNet().ZvsQCD, theWeight*LumiSF);
      theHistograms->fill("j1_PN_ZvsQCD"+histoType + cuts.at(i), "j1_PN_ZvsQCD"+histoType + cuts.at(i)+"; sublead. jet ParticleNet ZvsQCD", 20, 0, 1, recoV.daughter(1).particleNet().ZvsQCD, theWeight*LumiSF);
      theHistograms->fill("jj_PN_ZvsQCDSum"+histoType + cuts.at(i), "jj_PN_ZvsQCDSum"+histoType + cuts.at(i)+"; DiJet ParticleNet ZvsQCD", 20, 0, 1,  0.5*(recoV.daughter(0).particleNet().ZvsQCD + recoV.daughter(1).particleNet().ZvsQCD), theWeight*LumiSF);
      */
      //____________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________
      
      theHistograms->fill("j0_p(b)"+histoType + cuts.at(i), "j0_p(b)"+histoType + cuts.at(i)+"; lead. jet DeepFlavour p(b)", 20, 0, 1, recoV.daughter(0).deepFlavour().probb, theWeight*LumiSF);
      theHistograms->fill("j1_p(b)"+histoType + cuts.at(i), "j1_p(b)"+histoType + cuts.at(i)+"; sublead. jet DeepFlavour p(b)", 20, 0, 1, recoV.daughter(1).deepFlavour().probb, theWeight*LumiSF);
      theHistograms->fill("jj_p(b)Sum"+histoType + cuts.at(i), "jj_p(b)Sum"+histoType + cuts.at(i)+"; DiJet DeepFlavour p(b)", 20, 0, 1,  0.5*(recoV.daughter(0).deepFlavour().probb + recoV.daughter(1).deepFlavour().probb), theWeight*LumiSF);

      theHistograms->fill("j0_p(c)"+histoType + cuts.at(i), "j0_p(c)"+histoType + cuts.at(i)+"; lead. jet DeepFlavour p(c)", 20, 0, 1, recoV.daughter(0).deepFlavour().probc, theWeight*LumiSF);
      theHistograms->fill("j1_p(c)"+histoType + cuts.at(i), "j1_p(c)"+histoType + cuts.at(i)+"; sublead. jet DeepFlavour p(c)", 20, 0, 1, recoV.daughter(1).deepFlavour().probc, theWeight*LumiSF);
      theHistograms->fill("jj_p(c)Sum"+histoType + cuts.at(i), "jj_p(c)Sum"+histoType + cuts.at(i)+"; DiJet DeepFlavour p(c)", 20, 0, 1,  0.5*(recoV.daughter(0).deepFlavour().probc + recoV.daughter(1).deepFlavour().probc), theWeight*LumiSF);


      theHistograms->fill("j0_p(g)"+histoType + cuts.at(i), "j0_p(g)"+histoType + cuts.at(i)+"; lead. jet DeepFlavour p(g)", 20, 0, 1, recoV.daughter(0).deepFlavour().probg, theWeight*LumiSF);
      theHistograms->fill("j1_p(g)"+histoType + cuts.at(i), "j1_p(g)"+histoType + cuts.at(i)+"; sublead. jet DeepFlavour p(g)", 20, 0, 1, recoV.daughter(1).deepFlavour().probg, theWeight*LumiSF);
      theHistograms->fill("jj_p(g)Sum"+histoType + cuts.at(i), "jj_p(g)Sum"+histoType + cuts.at(i)+"; DiJet DeepFlavour p(g)", 20, 0, 1,  0.5*(recoV.daughter(0).deepFlavour().probg + recoV.daughter(1).deepFlavour().probg), theWeight*LumiSF);

      theHistograms->fill("j0_p(lepb)"+histoType + cuts.at(i), "j0_p(lepb)"+histoType + cuts.at(i)+"; lead. jet DeepFlavour p(lepb)", 20, 0, 1, recoV.daughter(0).deepFlavour().problepb, theWeight*LumiSF);
      theHistograms->fill("j1_p(lepb)"+histoType + cuts.at(i), "j1_p(lepb)"+histoType + cuts.at(i)+"; sublead. jet DeepFlavour p(lepb)", 20, 0, 1, recoV.daughter(1).deepFlavour().problepb, theWeight*LumiSF);
      theHistograms->fill("jj_p(lepb)Sum"+histoType + cuts.at(i), "jj_p(lepb)Sum"+histoType + cuts.at(i)+"; DiJet DeepFlavour p(lepb)", 20, 0, 1,  0.5*(recoV.daughter(0).deepFlavour().problepb + recoV.daughter(1).deepFlavour().problepb), theWeight*LumiSF);

      theHistograms->fill("j0_p(bb)"+histoType + cuts.at(i), "j0_p(bb)"+histoType + cuts.at(i)+"; lead. jet DeepFlavour p(bb)", 20, 0, 1, recoV.daughter(0).deepFlavour().probbb, theWeight*LumiSF);
      theHistograms->fill("j1_p(bb)"+histoType + cuts.at(i), "j1_p(bb)"+histoType + cuts.at(i)+"; sublead. jet DeepFlavour p(bb)", 20, 0, 1, recoV.daughter(1).deepFlavour().probbb, theWeight*LumiSF);
      theHistograms->fill("jj_p(bb)Sum"+histoType + cuts.at(i), "jj_p(bb)Sum"+histoType + cuts.at(i)+"; DiJet DeepFlavour p(bb)", 20, 0, 1,  0.5*(recoV.daughter(0).deepFlavour().probbb + recoV.daughter(1).deepFlavour().probbb), theWeight*LumiSF);

      theHistograms->fill("j0_p(uds)"+histoType + cuts.at(i), "j0_p(uds)"+histoType + cuts.at(i)+"; lead. jet DeepFlavour p(uds)", 20, 0, 1, recoV.daughter(0).deepFlavour().probuds, theWeight*LumiSF);
      theHistograms->fill("j1_p(uds)"+histoType + cuts.at(i), "j1_p(uds)"+histoType + cuts.at(i)+"; sublead. jet DeepFlavour p(uds)", 20, 0, 1, recoV.daughter(1).deepFlavour().probuds, theWeight*LumiSF);
      theHistograms->fill("jj_p(uds)Sum"+histoType + cuts.at(i), "jj_p(uds)Sum"+histoType + cuts.at(i)+"; DiJet DeepFlavour p(uds)", 20, 0, 1,  0.5*(recoV.daughter(0).deepFlavour().probuds + recoV.daughter(1).deepFlavour().probuds), theWeight*LumiSF);


      
      theHistograms->fill("j0_p(udsg)"+histoType + cuts.at(i), "j0_p(udsg)"+histoType + cuts.at(i)+"; lead. jet DeepFlavour p(udsg)", 20, 0, 1, recoV.daughter(0).deepFlavour().probuds + recoV.daughter(0).deepFlavour().probg, theWeight*LumiSF);
      theHistograms->fill("j1_p(udsg)"+histoType + cuts.at(i), "j1_p(udsg)"+histoType + cuts.at(i)+"; sublead. jet DeepFlavour p(udsg)", 20, 0, 1, recoV.daughter(1).deepFlavour().probuds + recoV.daughter(1).deepFlavour().probg, theWeight*LumiSF);
      theHistograms->fill("jj_p(udsg)Sum"+histoType + cuts.at(i), "jj_p(udsg)Sum"+histoType + cuts.at(i)+"; DiJet DeepFlavour p(udsg)", 20, 0, 1,  0.5*( (recoV.daughter(0).deepFlavour().probuds + recoV.daughter(1).deepFlavour().probuds) + (recoV.daughter(0).deepFlavour().probg + recoV.daughter(1).deepFlavour().probg) ), theWeight*LumiSF);

      theHistograms->fill("j0_p(uds)_vs_p(g)"+histoType + cuts.at(i), "j0_p(uds)_p(g)"+histoType + cuts.at(i)+"; lead. jet DeepFlavour p(uds); lead. jet DeepFlavour p(g)" ,  5, 0, 1, 5, 0, 1, recoV.daughter(0).deepFlavour().probuds, recoV.daughter(0).deepFlavour().probg, theWeight*LumiSF);
      theHistograms->fill("j1_p(uds)_vs_p(g)"+histoType + cuts.at(i), "j1_p(uds)_p(g)"+histoType + cuts.at(i)+"; sublead. jet DeepFlavour p(uds); sublead. jet DeepFlavour p(g)" ,  5, 0, 1, 5, 0, 1, recoV.daughter(1).deepFlavour().probuds, recoV.daughter(1).deepFlavour().probg, theWeight*LumiSF);

      
      
    }else if(VBTopo==-1){      
      theHistograms->fill("V_vs_Z_pt_" + histoType + cuts.at(i), "V_vs_Z_pt_" + histoType + cuts.at(i) +";V p_{t} [GeV/c]; #Z p_{t} [GeV/c]", 30, 0, 300, 30, 0, 300, recoFJ.pt(), Z->pt(), theWeight*LumiSF);
      theHistograms->fill("recoVMass_" + histoType + cuts.at(i), "mass of recoV", 40, 0, 200, recoFJ.mass(), (theWeight*LumiSF));
      //      std::cout<<recoFJ.deepFlavour().probb<<endl;
      //      printHistos(1, "sign", recoV, recoFJ, selectedphotons, VBTopo);
      //      theHistograms->fill("FJ_probb"+histoType + cuts.at(i), "FJ_probb"+histoType + cuts.at(i)+"; fat jet DeepFlavour", 16, -0.1, 0.3, recoFJ.deepFlavour().probb, theWeight*LumiSF);


      theHistograms->fill("recoVPt_" + histoType + cuts.at(i), "pt of recoV", 30, 0, 300, recoV.pt(), (theWeight*LumiSF));
      //theHistograms->fill("recoVEta_" + histoType + cuts.at(i), "eta of recoV", 30, 0, 3.5, fabs(recoV.eta()), (theWeight*LumiSF));
      //theHistograms->fill("recoVPhi_" + histoType + cuts.at(i), "phi of recoV", 30, -3.2, 3.2, recoV.phi(), (theWeight*LumiSF));
      //theHistograms->fill("recoVEnergy_" + histoType + cuts.at(i), "energy of  recoV", 60, 0, 2400, fabs(recoV.e()), (theWeight*LumiSF));
      theHistograms->fill("recoVDaughtersDeltaPhi_" + histoType + cuts.at(i), "dPhi of recoVDaughters", 20, 0, 3.2, fabs(physmath::deltaPhi(recoV.daughter(0).phi(), recoV.daughter(1).phi())), (theWeight*LumiSF));

    }

    theHistograms->fill("recoZDaughter0Pt_" + histoType + cuts.at(i), "pt of recoZDaughter0", 50, 0, 600, Z->daughter(0).pt(), (theWeight*LumiSF));
    theHistograms->fill("recoZDaughter1Pt_" + histoType + cuts.at(i), "pt of recoZDaughter1", 50, 0, 600, Z->daughter(1).pt(), (theWeight*LumiSF));
    


    /*
    theHistograms->fill("VZMass_" + histoType + cuts.at(i), "mass of reco VZ system", 40, 0, 400, lljj.mass(), (theWeight*LumiSF));
    theHistograms->fill("VZ_Pt_" + histoType + cuts.at(i), "pt of reco VZ system", 50, 0, 500, lljj.Pt(), (theWeight*LumiSF));
    */
    std::vector<phys::Jet> kinRECOjets;
    std::vector<phys::Jet> kinRECOfatJets;
    foreach (const phys::Jet &jet, *jets)
      {
	if (KinematicsOK(jet,ptcut,etacut))
	  kinRECOjets.push_back(jet);
      }
    theHistograms->fill("#recoJetsAK4_kinAcc_" + histoType + cuts.at(i), "#recoJetsAK4_kinematically_OK", 8, 0, 8, kinRECOjets.size(), theWeight*LumiSF);
    foreach (const phys::Jet &FJ, *jetsAK8)
      {
	if (KinematicsOK(FJ,ptcut,etacut))
	  kinRECOfatJets.push_back(FJ);
      }
    theHistograms->fill("#recoJetsAK8_kinAcc_" + histoType + cuts.at(i), "#recoJetsAK8_kinematically_OK", 8, 0, 8, kinRECOfatJets.size(), theWeight*LumiSF);

    if(VBTopo==1){
      theHistograms->fill("#recoJetsAK4_hadrTopo_"+ hadTopo +"_"+ histoType + cuts.at(i), "#recoJetsAK4_hadrTopo_"+ hadTopo +"_"+ histoType + cuts.at(i), 10, 0, 10, kinRECOjets.size(), theWeight*LumiSF);
      theHistograms->fill("#recoJetsAK8_hadrTopo_"+ hadTopo +"_"+ histoType + cuts.at(i), "#recoJetsAK8_hadrTopo_"+ hadTopo +"_"+ histoType + cuts.at(i), 10, 0, 10, kinRECOfatJets.size(), theWeight*LumiSF);
    }else if(VBTopo==-1){
      theHistograms->fill("#recoJetsAK4_hadrTopo_"+ hadTopo +"_"+ histoType + cuts.at(i), "#recoJetsAK4_hadrTopo_"+ hadTopo +"_"+ histoType + cuts.at(i), 10, 0, 10, kinRECOjets.size(), theWeight*LumiSF);
      theHistograms->fill("#recoJetsAK8_hadrTopo_"+ hadTopo +"_"+ histoType + cuts.at(i), "#recoJetsAK8_hadrTopo_"+ hadTopo +"_"+ histoType + cuts.at(i), 10, 0, 10, kinRECOfatJets.size(), theWeight*LumiSF);      
    }

    theHistograms->fill("#recoJets_total_"+ histoType + cuts.at(i), "#recoJets_total_"+ histoType + cuts.at(i), 10, 0, 10, kinRECOfatJets.size()+kinRECOfatJets.size(), theWeight*LumiSF);      

    //PlotJets(recoV.daughter(0),recoV.daughter(1), "", theWeight*LumiSF, histoType + cuts.at(i));

    theHistograms->fill("recoZMass_" + histoType + cuts.at(i), "mass of recoZ", 40, 0, 200, Z->mass(), (theWeight*LumiSF));

    theHistograms->fill("recoZPt_" + histoType + cuts.at(i), "pt of recoZ", 50, 0, 600, Z->pt(), (theWeight*LumiSF));
    theHistograms->fill("recoZEta_" + histoType + cuts.at(i), "eta of recoZ", 35, 0, 3.5, fabs(Z->eta()), (theWeight*LumiSF));
    //theHistograms->fill("recoZEnergy_" + histoType + cuts.at(i), "energy of  recoZ", 120, 0, 400, fabs(Z->e()), (theWeight*LumiSF));
    theHistograms->fill("recoZDeltaPhi_" + histoType + cuts.at(i), "dPhi of recoZ", 30, 0, 3.2, fabs(physmath::deltaPhi(Z->daughter(0).phi(), Z->daughter(1).phi())), (theWeight*LumiSF));


    if(VBTopo==1){
      theHistograms->fill(" HAD TOPO DJ cand mass_" + histoType + cuts.at(i), " HAD TOPO DJ cand mass cand mass_" + histoType + cuts.at(i)+" ; mjj Cand. [GeV]", 28, 50, 120, recoV.mass());
      theHistograms->fill(" HAD TOPO 2J/FJ cand mass_" + histoType + cuts.at(i), " HAD TOPO 2J/FJ cand mass_" + histoType + cuts.at(i)+" ; mVB Cand. [GeV]", 28, 50, 120, recoV.mass());
    }else if (VBTopo==-1){
      theHistograms->fill(" HAD TOPO FJ cand mass_" + histoType + cuts.at(i), " HAD TOPO FJ cand mass cand mass_" + histoType + cuts.at(i)+" ; mFj Cand. [GeV]", 28, 50, 120, recoFJ.mass());
      theHistograms->fill(" HAD TOPO 2J/FJ cand mass_" + histoType + cuts.at(i), " HAD TOPO 2J/FJ cand mass_" + histoType + cuts.at(i)+" ; mVB Cand. [GeV]", 28, 50, 120, recoFJ.mass());
    }

    theHistograms->fill(" HAD TOPO_" + histoType + cuts.at(i), " HAD TOPO_" + histoType + cuts.at(i), 3, -1.5, 1.5, VBTopo-0.5);


    
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
  //   theHistograms->fill("DeltaR_mostEnergeticPhoton_vs_BestMatchedRECOJet_"+histoType + cuts.at(i), "DeltaR_mostEnergeticPhoton_vs_BestMatchedRECOJet"+histoType + cuts.at(i)+"; #DeltaR", 50, 0, 5, fabs(physmath::deltaR(pair.first, pair.second)), theWeight*LumiSF);
  //   theHistograms->fill("DeltaR_vs_Deltapt_"+histoType + cuts.at(i), "DeltaR_vs_Deltapt"+histoType + cuts.at(i)+";#Delta pt [GeV/c] ; #DeltaR", 20, -100, 100, 20, 0, 2, pair.first.pt()-pair.second.pt(),fabs(physmath::deltaR(pair.first, pair.second)), theWeight*LumiSF);
  // }



    phys::Photon mostEnergeticPhoton;
    //    if (selectedphotons.size() > 0)
    //  {
    std::stable_sort(selectedphotons.begin(), selectedphotons.end(), phys::EComparator());
    mostEnergeticPhoton = selectedphotons[0];
    //  }

    theHistograms->fill("pt_mostenergeticphoton"+histoType + cuts.at(i), "pt_mostenergeticphoton"+histoType + cuts.at(i), 30, 0, 300, mostEnergeticPhoton.pt(), (theWeight*LumiSF));

    std::pair<phys::Photon, phys::Jet> nearestRECOjetstoPhoton;

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

    std::vector<TLorentzVector> llG;
    llG.push_back(Z->daughter(0).p4());
    llG.push_back(Z->daughter(1).p4());
    llG.push_back(selectedphotons.at(0).p4());


    
    if (VBTopo == 1){
      theHistograms->fill("DR_Jets_"+histoType + cuts.at(i), "DR_Jets_"+histoType + cuts.at(i)+"; jets #DeltaR", 50, 0, 5, fabs(physmath::deltaR(recoV.daughter(0),recoV.daughter(1))), theWeight*LumiSF);
      if(fabs(physmath::deltaR(recoV.daughter(0), mostEnergeticPhoton))<fabs(physmath::deltaR(recoV.daughter(1), mostEnergeticPhoton)) )
	nearestRECOjetstoPhoton={mostEnergeticPhoton, recoV.daughter(0)};
      else if (fabs(physmath::deltaR(recoV.daughter(0), mostEnergeticPhoton))>fabs(physmath::deltaR(recoV.daughter(1), mostEnergeticPhoton)) )
	nearestRECOjetstoPhoton={mostEnergeticPhoton, recoV.daughter(1)};


      theHistograms->fill("DR_gammaClosestJet_"+histoType + cuts.at(i), "DR_gammaClosestJet_"+histoType + cuts.at(i)+"; #DeltaR", 50, 0, 5, fabs(physmath::deltaR(nearestRECOjetstoPhoton.first, nearestRECOjetstoPhoton.second)), theWeight*LumiSF);
      theHistograms->fill("DeltaR_vs_Deltapt_gammaJet"+histoType + cuts.at(i), "DeltaR_vs_Deltapt_gammaJet"+histoType + cuts.at(i)+";#Delta pt [GeV/c] ; #DeltaR", 20, -100, 100, 50, 0, 5, nearestRECOjetstoPhoton.first.pt()-nearestRECOjetstoPhoton.second.pt(),fabs(physmath::deltaR(nearestRECOjetstoPhoton.first, nearestRECOjetstoPhoton.second)), theWeight*LumiSF);

    }
      

    std::pair<phys::Photon, phys::Lepton> nearestChLeptToPhoton;

    if(fabs(physmath::deltaR(Z->daughter(0), mostEnergeticPhoton))<fabs(physmath::deltaR(Z->daughter(1), mostEnergeticPhoton)) )
      nearestChLeptToPhoton={mostEnergeticPhoton, Z->daughter(0)};
    else if (fabs(physmath::deltaR(Z->daughter(0), mostEnergeticPhoton))>fabs(physmath::deltaR(Z->daughter(1), mostEnergeticPhoton)) )
      nearestChLeptToPhoton={mostEnergeticPhoton, Z->daughter(1)};

    theHistograms->fill("DR_Lept_"+histoType + cuts.at(i), "DR_Lept_"+histoType + cuts.at(i)+"; leptons #DeltaR", 50, 0, 5, fabs(physmath::deltaR(Z->daughter(0),Z->daughter(1))), theWeight*LumiSF);

    theHistograms->fill("DR_gammaClosestLept_"+histoType + cuts.at(i), "DR_gammaClosestLept_"+histoType + cuts.at(i)+"; #DeltaR", 50, 0, 5, fabs(physmath::deltaR(nearestChLeptToPhoton.first, nearestChLeptToPhoton.second)), theWeight*LumiSF);


    /*
    theHistograms->fill("j0_deepFlavour"+histoType + cuts.at(i), "j0_deepFlavour"+histoType + cuts.at(i)+"; leptons #DeltaR", 50, 0, 5, fabs(physmath::deltaR(Z->daughter(0),Z->daughter(1))), theWeight*LumiSF);
    */
    
    //if(selectedphotons.size()>0)
    //  {


    /*
	theHistograms->fill("DRlGs_vs_mllG_"+histoType + cuts.at(i), "DRjGs_vs_mllG_"+histoType + cuts.at(i)+"; mll#gamma [GeV] ; #DeltaRl#gamma", 20, 0, 200, 50, 0, 5, nearestChLeptToPhoton.first.pt()-nearestChLeptToPhoton.second.pt(),fabs(physmath::deltaR(nearestChLeptToPhoton.first, nearestChLeptToPhoton.second)), theWeight*LumiSF);
	theHistograms->fill("DRjGs_vs_mjjG_"+histoType + cuts.at(i), "DRjGs_vs_mjjG_"+histoType + cuts.at(i)+"; mjj#gamma [GeV] ; #DeltaRj#gamma", 20, -100, 100, 50, 0, 5, nearestRECOjetstoPhoton.first.pt()-nearestRECOjetstoPhoton.second.pt(),fabs(physmath::deltaR(nearestRECOjetstoPhoton.first, nearestRECOjetstoPhoton.second)), theWeight*LumiSF);
	*/

    if(VBTopo==1)
      {
	theHistograms->fill("mjjG_"+histoType + cuts.at(i), 30, 50, 350, mjjPh, theWeight*LumiSF);
	theHistograms->fill("mjj_vs_mjjG_"+histoType + cuts.at(i), "mjj_vs_mjjG_"+histoType + cuts.at(i)+"; mjj [GeV] ; mjj#gamma [GeV]", 35, 50, 120, 30, 50, 350, mjj, mjjPh, theWeight*LumiSF);
	theHistograms->fill("mll_vs_mjj_"+histoType + cuts.at(i), "mll_vs_mjj_"+histoType + cuts.at(i)+"; mll [GeV] ; mjj [GeV]", 30, 60, 120, 35, 50, 120, mll, mjj, theWeight*LumiSF);
	theHistograms->fill("mllG_vs_mjjG_"+histoType + cuts.at(i), "mllG_vs_mjjG_"+histoType + cuts.at(i)+"; mll#gamma [GeV] ; mjj#gamma [GeV]", 80, 50, 450, 30, 50, 350, mllPh, mjjPh, theWeight*LumiSF);
	theHistograms->fill("mllG_vs_mjj_"+histoType + cuts.at(i), "mllG_vs_mjj_"+histoType + cuts.at(i)+"; mll#gamma [GeV] ; mjj [GeV]", 80, 50, 450, 35, 50, 120, mllPh, mjj, theWeight*LumiSF);
	theHistograms->fill("mll_vs_mjjG_"+histoType + cuts.at(i), "mll_vs_mjjG_"+histoType + cuts.at(i)+"; mll [GeV] ; mjj#gamma [GeV]", 30, 60, 120, 30, 50, 350, mllPh, mjjPh, theWeight*LumiSF);


	theHistograms->fill("DRlGs_vs_DRjG_"+histoType + cuts.at(i), "DRlGs_vs_DRjG_"+histoType + cuts.at(i)+"; #DeltaRj#gamma [GeV] ; #DeltaRl#gamma", 50, 0, 5, 50, 0, 5, fabs(physmath::deltaR(nearestRECOjetstoPhoton.first, nearestRECOjetstoPhoton.second)),fabs(physmath::deltaR(nearestChLeptToPhoton.first, nearestChLeptToPhoton.second)), theWeight*LumiSF);

	theHistograms->fill("DALITZ_PLOT_VZG_"+histoType + cuts.at(i), "DALITZ_PLOT_VZG_"+histoType + cuts.at(i)+"; m^2 ll#gamma [GeV] ; m^2  jj#gamma [GeV]", 50, 0, 100000, 50, 0, 100000, m2llPh, m2jjPh, theWeight*LumiSF);

	theHistograms->fill("DALITZ_PLOT_llG_"+histoType + cuts.at(i), "DALITZ_PLOT_llG_"+histoType + cuts.at(i)+"; m^2 l0#gamma [GeV] ; m^2  l1#gamma [GeV]", 50, 0, 50000, 50, 0, 50000, m2l0Ph, m2l1Ph, theWeight*LumiSF);

	theHistograms->fill("relativeLLGmass_vs_cosPhi_ll_"+histoType + cuts.at(i), "relativeLLGmass_vs_cosPhi_ll_"+histoType + cuts.at(i)+"; 2m^2 ll#gamma/(ptl0+ptl1) [GeV]; cos #phi ll", 100, 0, 5000, 10, -1, 1, 2*m2llPh/(Z->daughter(0).pt()+Z->daughter(1).pt()),TMath::Cos(fabs(physmath::deltaPhi(Z->daughter(0).phi(), Z->daughter(1).phi()))), theWeight*LumiSF);

	theHistograms->fill("relativeLLGmass_vs_cosPhi_l0Ph_"+histoType + cuts.at(i), "relativeLLGmass_vs_cosPhi_l0Ph_"+histoType + cuts.at(i)+"; 2m^2 ll#gamma/(ptl0+pt#gamma) [GeV]; cos #phi l0#gamma", 100, 0, 5000, 10, -1, 1, 2*m2llPh/(Z->daughter(0).pt()+selectedphotons.at(0).pt()),TMath::Cos(fabs(physmath::deltaPhi(Z->daughter(0).phi(), selectedphotons.at(0).phi()))), theWeight*LumiSF);
	theHistograms->fill("relativeLLGmass_vs_cosPhi_l1Ph_"+histoType + cuts.at(i), "relativeLLGmass_vs_cosPhi_l1Ph_"+histoType + cuts.at(i)+"; 2m^2 ll#gamma/(ptl1+pt#gamma) [GeV]; cos #phi l1#gamma", 100, 0, 5000, 10, -1, 1, 2*m2llPh/(Z->daughter(1).pt()+selectedphotons.at(0).pt()),TMath::Cos(fabs(physmath::deltaPhi(Z->daughter(1).phi(), selectedphotons.at(0).phi()))), theWeight*LumiSF);

      }    

    theHistograms->fill("mllG_"+histoType + cuts.at(i), 45, 0, 450, mllPh, theWeight*LumiSF);
    theHistograms->fill("mll_vs_mllG_"+histoType + cuts.at(i), "mll_vs_mllG_"+histoType + cuts.at(i)+"; mll [GeV] ; mll#gamma [GeV]", 30, 60, 120, 80, 50, 450, mll, mllPh, theWeight*LumiSF);
    theHistograms->fill("mllG_vs_DRlGs_"+histoType + cuts.at(i), "mllG_vs_DRlGs_"+histoType + cuts.at(i)+"; mll#gamma [GeV] ; #DeltaRl#gamma", 60, 150, 450, 50, 0, 5, mllPh, fabs(physmath::deltaR(nearestChLeptToPhoton.first, nearestChLeptToPhoton.second)), theWeight*LumiSF);



    
    //theHistograms->fill("VZGMass_" + histoType + cuts.at(i), "mass of reco VZG system", 40, 0, 400, lljjG.p4().M(), (theWeight*LumiSF));

    int kinPhotonsCounter=0;
    int kinNoVLPhotonsCounter=0;
    int VLPhotonsCounter=0;
    int VLNoLoosePhotonsCounter=0;
    int LoosePhotonsCounter=0;

    foreach (auto p , *photons)
      {
	if (p.id() == 22 && KinematicsOK(p, 20, 2.4) && !p.hasPixelSeed() && p.passElectronVeto())
	  {
	    kinPhotonsCounter++;

	    if(p.cutBasedID(Photon::IdWp::VeryLoose))
	      {
		VLPhotonsCounter++;
		if (p.cutBasedIDLoose())	LoosePhotonsCounter++;
		else VLNoLoosePhotonsCounter++;
	      }
	    else kinNoVLPhotonsCounter++;
	  }
      }
    theHistograms->fill("#_gamma_kin_" + histoType, "#_gamma_kin_", 12, 0, 12, kinPhotonsCounter, (theWeight*LumiSF));
    theHistograms->fill("#_gamma_kinButNotVL_" + histoType, "#_gamma_kinButNotVL_", 12, 0, 12, kinNoVLPhotonsCounter, (theWeight*LumiSF));
    theHistograms->fill("#_gamma_VL_" + histoType, "#_gamma_VL_", 12, 0, 12, VLPhotonsCounter, (theWeight*LumiSF));
    theHistograms->fill("#_gamma_VLButNotLoose_" + histoType, "#_gamma_VLButNotLoose_", 12, 0, 12, VLNoLoosePhotonsCounter, (theWeight*LumiSF));
    theHistograms->fill("#_gamma_Loose_" + histoType, "#_gamma_Loose_", 12, 0, 12, LoosePhotonsCounter, (theWeight*LumiSF));

    if(VBTopo==1){
      for(int l = 0; l<cutsToApply; l++)
	{
	  theHistograms->fill("FWM_T"+orders.at(l)+"_jets_"+histoType+cuts.at(i), "FWM_T"+orders.at(l)+"_jets_"+histoType+cuts.at(i)+"; H_"+orders.at(l)+"^T jets + gamma", 40, -1, 3,
			      SumFWM(l, 't', jjG), theWeight*LumiSF);

	  theHistograms->fill("FWM_T"+orders.at(l)+"_fullSyst_"+histoType+cuts.at(i), "FWM_T"+orders.at(l)+"_fullSyst_"+histoType+cuts.at(i)+"; H_"+orders.at(l)+"^T jets and gamma", 40, -1, 3,
			      SumFWM(l, 't', lljjG), theWeight*LumiSF);
	}
    }
    printHistos(++i, histoType, recoV, recoFJ, selectedphotons,VBTopo); 
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


  
  bool GENsignal=IN_GENsignalDef();//(HadronicSignalConstraint() && LeptonicSignalConstraint() && PhotonSignalConstraint());
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
  /*
  std::vector<phys::Particle> selectedGENphotons;
  for (auto p : *genParticles)
    if (p.id() == 22 && KinematicsOK(p, 20, 2.4) && p.genStatusFlags().test(phys::isPrompt) &&  p.genStatusFlags().test(phys::fromHardProcess))
      selectedGENphotons.push_back(p);
  if(verbose==true) std::cout<< "Number of selected gen photons = "<<selectedGENphotons.size()<<std::endl;

  
  if (selectedGENphotons.size()>=1 && GENsignal)
    {
      TLorentzVector GEN_jjPh, GEN_llPh ;
      double mGEN_llPh, mGEN_jjPh,  mGenZ, mGenV;

      if(genVBHelper_.WtoQ().size()>=1 || genVBHelper_.ZtoQ().size()>=1)
	{
	  if(genVBHelper_.WtoQ().size()>=1)
	    {
	      GEN_jjPh = genVBHelper_.WtoQ()[0].daughter(0).p4()+genVBHelper_.WtoQ()[0].daughter(1).p4()+selectedGENphotons.at(0).p4();
	      mGenV=genVBHelper_.WtoQ()[0].mass();
	      std::cout<< "GEN_WHad mass implemented"<<std::endl;
	      }
	  if(genVBHelper_.ZtoQ().size()>=1)
	    {
	      GEN_jjPh = genVBHelper_.ZtoQ()[0].daughter(0).p4()+genVBHelper_.WtoQ()[0].daughter(1).p4()+selectedGENphotons.at(0).p4();
	      mGenV=genVBHelper_.ZtoQ()[0].mass();
	      std::cout<< "GEN_ZHad mass implemented"<<std::endl;
	    }
	  mGEN_jjPh=GEN_jjPh.M();
	
	  if(genVBHelper_.ZtoChLep().size()>=1)
	    {
	      std::cout<< "entered ZToL"<<std::endl;
	      GEN_llPh = genVBHelper_.ZtoChLep()[0].daughter(0).p4()+genVBHelper_.ZtoChLep()[0].daughter(0).p4()+selectedGENphotons.at(0).p4();
	      std::cout<< "GEN_llPh implemented"<<std::endl;
	      mGEN_llPh=GEN_llPh.M();
	      std::cout<< "mass of GEN_llPh implemented"<<std::endl;
	      mGenZ=genVBHelper_.ZtoChLep()[0].mass();
	      std::cout<< "GEN_Z mass implemented"<<std::endl;

	      theHistograms->fill("GEN mjj_vs_mjjG", "GEN mjj_vs_mjjG; mjj [GeV] ; mjj#gamma [GeV]", 35, 50, 120, 30, 50, 350, mGenV, mGEN_jjPh, theWeight*LumiSF);
	      theHistograms->fill("GEN mll_vs_mllG", "GEN mll_vs_mllG; mll [GeV] ; mll#gamma [GeV]", 30, 60, 120, 20, 50, 150, mGenZ, mGEN_llPh, theWeight*LumiSF);

	      std::cout<< "Histos filled"<<std::endl;
	      
	    }
	}
    }
  */
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
      if (fabs(physmath::deltaR(quark, nearestjet)) < 0.4)
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
    theHistograms->fill("DeltaR_quark_vs_BestMatchedGENJet", "DeltaR_quark_vs_BestMatchedGENJet; #DeltaR", 20, 0, 0.5, fabs(physmath::deltaR(pair.first, pair.second)), theWeight*LumiSF);
    theHistograms->fill("DeltaR_quark_jet_vs_pt", "DeltaR vs pt;pt [GeV/c] ; #DeltaR", 10, ptcut, 300, 20, 0, 0.2, pair.first.pt(),fabs(physmath::deltaR(pair.first, pair.second)), theWeight*LumiSF);

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
      if (fabs(physmath::deltaR(genJet, nearestRECOjet)) < 0.4)
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
    theHistograms->fill("DeltaR_GENjet_vs_BestMatchedRECOJet", "DeltaR_GENjet_vs_BestMatchedRECOJet; #DeltaR", 20, 0, 0.5, fabs(physmath::deltaR(pair.first, pair.second)), theWeight*LumiSF);
    theHistograms->fill("DeltaR_jets_vs_pt", "DeltaR jets vs pt;pt [GeV/c] ; #DeltaR", 10, ptcut, 300, 20, 0, 0.2, pair.first.pt(),fabs(physmath::deltaR(pair.first, pair.second)), theWeight*LumiSF);
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
      double deltaR1 = fabs(physmath::deltaR(Diquark.daughter(0), genJet));
      std::cout << "deltaR1= " << deltaR1 << std::endl;
      double deltaR2 = fabs(physmath::deltaR(Diquark.daughter(1), genJet));
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
      double deltaR1 = fabs(physmath::deltaR(DiJet.daughter(0), recoJet));
      std::cout << "deltaR1= " << deltaR1 << std::endl;
      double deltaR2 = fabs(physmath::deltaR(DiJet.daughter(1), recoJet));
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
        double deltaRAA = fabs(physmath::deltaR(jetGENA, jetRECOA));
        double deltaRAB = fabs(physmath::deltaR(jetGENA, jetRECOB));
        double deltaRBA = fabs(physmath::deltaR(jetGENB, jetRECOA));
        double deltaRBB = fabs(physmath::deltaR(jetGENB, jetRECOB));
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
  // double JJdeltaR = fabs(physmath::deltaR(gen, reco));

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
      if (fabs(physmath::deltaR(quark, nearestGENjet)) < 0.4)
      {
        GENjetsfromquarks.push_back(nearestGENjet);
        makesGENjet = true;
	if (quarkMatchingCounter == 0) firstGENjetMatched = nearestGENjet;
        else	theHistograms->fill("overlappedQuarks", "overlappedQuarks", 2, 0, 2, quarkMatchingCounter == 1 && fabs(physmath::deltaR(firstGENjetMatched, nearestGENjet))<0.4, theWeight*LumiSF);
	quarkMatchingCounter++;
      }
      else       theHistograms->fill("dR_unmatchedQuarks_closestGENjet", "dR_unmatchedQuarks_closestGENjet", 50, 0, 5, fabs(physmath::deltaR(quark, nearestGENjet)), theWeight*LumiSF);


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
    theHistograms->fill("DeltaR_quark_vs_BestMatchedGENJet", "DeltaR_quark_vs_BestMatchedGENJet; #DeltaR", 20, 0, 0.5, fabs(physmath::deltaR(pair.first, pair.second)), theWeight*LumiSF);
    theHistograms->fill("DeltaR_quark_GENjet_vs_pt", "DeltaR vs pt;pt [GeV/c] ; #DeltaR", 10, 0, 300, 20, 0, 0.2, pair.first.pt(),fabs(physmath::deltaR(pair.first, pair.second)), theWeight*LumiSF);

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
      theHistograms->fill("dR_GENclosestRECOjet", "dR_GENclosestRECOjet", 500, 0, 5, fabs(physmath::deltaR(genJet, nearestRECOjet)), theWeight*LumiSF);
      if (fabs(physmath::deltaR(genJet, nearestRECOjet)) < 0.4)
      {
        RECOjetsfromGENjets.push_back(nearestRECOjet);
        isreconstructed = true;
	GENtoRECOjetsCounter++;
      }
      else       theHistograms->fill("dR_unmatchedGENclosestRECOjet", "dR_unmatchedGENclosestRECOjet", 50, 0, 5, fabs(physmath::deltaR(genJet, nearestRECOjet)), theWeight*LumiSF);

    }
    if (isreconstructed && selectedRECOjets.size() > 0)       theHistograms->fill("Pt_genJet_num", " Pt_genJet_num; GeV/c", 10, 0, 300, genJet.pt(), theWeight*LumiSF);
    theHistograms->fill("#GEN=>RECO", "#GEN=>RECO", 2, 0, 2, isreconstructed && selectedRECOjets.size() > 0, theWeight*LumiSF);
  }
  theHistograms->fill("2GENJetsToRECO", "2GENJetsToRECO", 2, 0, 2, GENtoRECOjetsCounter==2, theWeight*LumiSF);
  theHistograms->fill("AtLeast1GENJetToRECO", "AtLeast1GENJetToRECO", 2, 0, 2, GENtoRECOjetsCounter>0, theWeight*LumiSF);

  
  for (auto pair : nearestRECOjetstoGENjets)
  {
    ResolutionPlots(pair.first,pair.second,"SingleJetsReconstruction_",theWeight*LumiSF,"");
    theHistograms->fill("DeltaR_GENjet_vs_BestMatchedRECOJet", "DeltaR_GENjet_vs_BestMatchedRECOJet; #DeltaR", 20, 0, 0.5, fabs(physmath::deltaR(pair.first, pair.second)), theWeight*LumiSF);
    theHistograms->fill("DeltaR_jets_vs_pt", "DeltaR jets vs pt;pt [GeV/c] ; #DeltaR", 10, 0, 300, 20, 0, 0.2, pair.first.pt(),fabs(physmath::deltaR(pair.first, pair.second)), theWeight*LumiSF);
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
       double JJdeltaR = fabs(physmath::deltaR(Jet0, Jet1));

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


