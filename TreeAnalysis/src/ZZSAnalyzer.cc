#include "VVXAnalysis/TreeAnalysis/interface/ZZSAnalyzer.h"
#include "VVXAnalysis/Commons/interface/SignalDefinitions.h"

#include <boost/foreach.hpp>
#define foreach BOOST_FOREACH

using namespace phys;
//using namespace zz;

void ZZSAnalyzer::analyze(){
  

   foreach(const Boson<Lepton> zmm, *Zmm)
     theHistograms.fill("ZmmMass",100,65,115,zmm.p4().M(),theWeight);
  
  foreach(const Boson<Jet> w, *Wjj)
    if(w.daughter(0).pt() > 40 && w.daughter(1).pt()>40)    
      theHistograms.fill("WjjMass",80,40,120,w.p4().M(),theWeight);
  
  foreach(const Boson<Electron> zee, *Zee)
    theHistograms.fill("ZeeMass",100,65,115,zee.p4().M(),theWeight);
  
  foreach(const Jet jet, *jets)
    theHistograms.fill("jet_Eta",100,-6,6, jet.eta(),theWeight);

  
  
  std::tuple<bool,Boson<Lepton>, Boson<Lepton> > ZZcand = zz::zz4l(*Zmm,*Zee);
  theHistograms.fill("ZZcand",5,-2,2,std::get<0>(ZZcand),theWeight);
  
  if(std::get<0>(ZZcand)){//da mettere poi la condizione sul best candidate ZZ
   
    const int size = jets->size();

    double jet_pt[size];
    double jet_mass = 0;
    double jet_mass_check = 0;
    double jet_deltaeta =0;
    double jet_deltay =0;

    TLorentzVector jet_p4[size];
    TLorentzVector lead_jet1;
    TLorentzVector lead_jet2;

    int i= 0;
    int j =0;
     
    theHistograms.fill("JetMult",10,0,10,size,theWeight);
  
  

  //#events vs M_jj - DeltaEta - DeltaY (eta<4.7, pt>30)
      foreach(const Jet jet, *jets){
      
      jet_pt[i] = jet.pt();
      jet_p4[i] = jet.p4();

      //std::cout << jet_pt[i] << std::endl;

      //std::cout  << " " << jet_p4[i][0]  << " " << jet_p4[i][1] << " " << jet_p4[i][2] << " " << jet_p4[i][3]<< " " << jets->size()<< " " << i << " " << jet_pt[i] << " " << jet.eta() << std::endl;  
      
      i++;

    
    }
    
      
    std::sort(jet_pt,jet_pt + size);  
    
    if(i!=0 && jet_pt[size-1]!=0 && jet_pt[size-2]!=0){
      
      //std::cout << "lead_jet1_pt= " << jet_pt[size-1] <<  " lead_jet2_pt= " << jet_pt[size-2] << " " << std::endl;
      
      for(j = 0; j < size; j++){
	
	if(jet_p4[j].Pt() == jet_pt[size-1]) lead_jet1 = jet_p4[j];
	else if(jet_p4[j].Pt() == jet_pt[size-2]) lead_jet2 = jet_p4[j];
	else continue;
		
      }

      jet_mass = (lead_jet1+lead_jet2).M();
      jet_deltaeta = TMath::Abs(lead_jet1.PseudoRapidity() - lead_jet2.PseudoRapidity());
      jet_deltay = TMath::Abs(lead_jet1.Rapidity() - lead_jet2.Rapidity());

      //jet_mass_check = TMath::Sqrt(lead_jet1.M()*lead_jet1.M()+lead_jet2.M()*lead_jet2.M()+2*(lead_jet1[3]*lead_jet2[3] - lead_jet1[0]*lead_jet2[0] - lead_jet1[1]*lead_jet2[1] - lead_jet1[2]*lead_jet2[2]));
      
      //std::cout<< "M_jj = " << jet_mass << " M1= " << lead_jet1.M() <<  " M2= " << lead_jet2.M() << " " << jet_deltaeta << std::endl; 
      
      
      if(jet_pt[size-2]>30){ 
	
	theHistograms.fill("M_jj",500,0,1000,jet_mass,theWeight);
	theHistograms.fill("DeltaEta_jj",50,0,5,jet_deltaeta,theWeight);
	theHistograms.fill("DeltaY_jj",50,0,5,jet_deltay,theWeight);
     
      }
    }

    
    //#events vs M_jj  - DeltaEta - DeltaY (eta<2.5, pt>30)
    const int csize = centralJets->size();

    // std::cout << "size: " << size << " csize: " << csize << std::endl;

    double cjet_pt[csize];
    double cjet_mass = 0;
    double cjet_mass_check = 0;
    double cjet_deltaeta =0;
    double cjet_deltay =0;

    TLorentzVector cjet_p4[csize];
    TLorentzVector lead_cjet1;
    TLorentzVector lead_cjet2;

    int i_c= 0;
    int j_c =0;
    // int k_c =0;

    theHistograms.fill("JetMultcentral",10,0,10,csize,theWeight);

    foreach(const Jet cjet, *centralJets){
       
      cjet_pt[i_c] = cjet.pt();
      cjet_p4[i_c] = cjet.p4();
      
      //std::cout  << " " << cjet_p4[i_c][0]  << " " << cjet_p4[i_c][1] << " " << cjet_p4[i_c][2] << " " << cjet_p4[i_c][3]<< " " << csize<< " " << i_c << " " << cjet_pt[i_c] << " " << cjet.eta() << std::endl;  
      
      i_c++;
    }
    
    std::sort(cjet_pt,cjet_pt + csize);  
    
    if(i_c!=0 && cjet_pt[csize-1]!=0 && cjet_pt[csize-2]!=0){
      
      //std::cout << "lead_cjet1_pt= " << cjet_pt[csize-1] <<  " lead_cjet2_pt= " << cjet_pt[csize-2] << " " << std::endl;
      
      for(j_c = 0; j_c < csize; j_c++){
	
	if(cjet_p4[j_c].Pt() == cjet_pt[csize-1]) lead_cjet1 = cjet_p4[j_c];
	else if(cjet_p4[j_c].Pt() == cjet_pt[csize-2]) lead_cjet2 = cjet_p4[j_c];
	else continue;
	
      }
      
      cjet_mass = (lead_cjet1+lead_cjet2).M();
      cjet_deltaeta = TMath::Abs(lead_cjet1.PseudoRapidity() - lead_cjet2.PseudoRapidity());
      cjet_deltay = TMath::Abs(lead_cjet1.Rapidity() - lead_cjet2.Rapidity());
      
      cjet_mass_check = TMath::Sqrt(lead_cjet1.M()*lead_cjet1.M()+lead_cjet2.M()*lead_cjet2.M()+2*(lead_cjet1[3]*lead_cjet2[3] - lead_cjet1[0]*lead_cjet2[0] - lead_cjet1[1]*lead_cjet2[1] - lead_cjet1[2]*lead_cjet2[2]));
      
      //std::cout<< "M_jj = " << cjet_mass << " M1= " << lead_cjet1.M() <<  " M2= " << lead_cjet2.M() << " " << cjet_deltaeta << std::endl; 
      
      //#events vs M_jj (eta<4.7, pt>30)
      if(cjet_pt[csize-2]>30) {

	theHistograms.fill("M_jj_central",500,0,1000,cjet_mass,theWeight);
	theHistograms.fill("DeltaEta_jj_central",50,0,5,cjet_deltaeta,theWeight);
	theHistograms.fill("DeltaY_jj_central",50,0,5,cjet_deltay,theWeight);
      }
    }
  }
  























  // foreach(const DiBoson< Boson<Lepton>, Boson<Lepton>> zz4m, *ZZ4m)
//     theHistograms.fill("ZZ4m",100,0,300,(zz4m.first().p4()+zz4m.second().p4()).M())


//   }

//   if(!std::get<0>(ZZ)) theHistograms.fill("ZZCand",5,-2,2,0);
//   else theHistograms.fill("ZZCand",5,-2,2,1);
  
//   if(!std::get<0>(ZZ)) std::cout << "-1" ;
//   else std::cout << "0" ;

  //if(!std::get<0>(ZZ)) std::cout << "-1" ;

  // if(!std::get<0>(ZZ)) std::cout << "-1" <
  // < std::endl;//theHistograms.fill("ZZCand",,150,300,zz.p4().M());
  

    // std::cout << "-1" << std::endl;

 //zz, *ZZ)
    // theHistograms.fill("ZZMass",100,150,300,zz.p4().M());

  //const DiBoson< Boson<Lepton>, Boson<Lepton>> zz4m = *ZZ4m;

//   //foreach(const DiBoson< Boson<Lepton>, Boson<Lepton>> zz4m, *ZZ4m)
// foreach(const DiBoson zz4m, *ZZ4m)
//   //  if(zz4m.passFullSel() == true)    
//   theHistograms.fill("ZZ4m",100,0,300,(zz4m.first().p4()+zz4m.second().p4()).M());
  
}
