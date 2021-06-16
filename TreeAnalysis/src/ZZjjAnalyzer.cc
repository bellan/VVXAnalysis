#include "VVXAnalysis/TreeAnalysis/interface/ZZjjAnalyzer.h"
#include "VVXAnalysis/Commons/interface/SignalDefinitions.h"
#include "VVXAnalysis/Commons/interface/Utilities.h"
#include "VVXAnalysis/Commons/interface/Comparators.h"
#include "VVXAnalysis/Commons/interface/Colours.h"
#include <TSpline.h>


#include <boost/foreach.hpp>
#define foreach BOOST_FOREACH

#include <boost/assign/std/vector.hpp> 
using namespace boost::assign;

using namespace colour;
using namespace physmath;
using namespace phys;
using namespace std;


void ZZjjAnalyzer::GenAnalysis(DiBosonParticle &genZZ, Particle &Jet0, Particle &Jet1){
  // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  // ~~~~~~~~~~~~~~~~~~~~ Begin of gen Analysis ~~~~~~~~~~~~~~~~~~~~
  // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  int cut = 0;  
  theHistograms.fill("GenCuts_wei", "Weighted counters", 10, -0.5, 9.5, cut, theWeight);
  
  // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Jets ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  // ----- jets' number > 2 and |eta| < 4.7
  helper_->FindLeadingJets(genJets, Jet0, Jet1, genParticles);

  if(Jet0.p4().Mt() == 0){
    return;
  }
  cut++;  
  theHistograms.fill("GenCuts_wei", "Weighted counters", 10, -0.5, 9.5, cut, theWeight);

  

  // ~~~~~~~~~~~~~~~~~~~~~~~ ZZ candidates ~~~~~~~~~~~~~~~~~~~~~~~~~
  zz::SignalTopology topology = zz::getSignalTopology(*genParticles, *genJets, *genJetsAK8);
  BosonParticle Z0 = get<1>(topology);
  BosonParticle Z1 = get<2>(topology);

  bitset<16> top = (bitset<16>)get<0>(topology);
  if(top.test(0) == 0 || top.test(1) == 1 || top.test(5) == 0 || top.test(6) == 0 || top.test(7) == 0){
    return;
  }

  eventGen++;
  weightGen += theWeight;
  
  cut++;  
  theHistograms.fill("GenCuts_wei", "Weighted counters", 10, -0.5, 9.5, cut, theWeight);
  theHistograms.fill("ZZ_Events",   "Weighted counters", 10, -0.5, 9.5,   2, theWeight);
  
  
  
  // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ All ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  TLorentzVector llllp4 = Z0.daughter(0).p4() + Z0.daughter(1).p4() + Z1.daughter(0).p4() + Z1.daughter(1).p4();
  TLorentzVector jjp4 = Jet0.p4() + Jet1.p4();
  TLorentzVector ZZjjp4 = Z0.p4() + Z1.p4() + Jet0.p4() + Jet1.p4();


  // Histograms before cuts
  theHistograms.fill("Gen4l_mass_BC",   "m_{4l}",   300, 0, 3000, llllp4.M(), theWeight);
  theHistograms.fill("GenJJ_mass_BC",   "m_{jj}",   400, 0, 4000, jjp4.M(),   theWeight);
  theHistograms.fill("GenZZJJ_mass_BC", "m_{ZZjj}", 500, 0, 5000, ZZjjp4.M(), theWeight);
  
  

  // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Cuts ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  if(llllp4.M() < 180.)
    return;
  
  cut++;  
  theHistograms.fill("GenCuts_wei", "Weighted counters", 10, -0.5, 9.5, cut, theWeight);

  if(jjp4.M() < 100.)
    return;
  
  cut++;  
  theHistograms.fill("GenCuts_wei", "Weighted counters", 10, -0.5, 9.5, cut, theWeight);
  theHistograms.fill("ZZ_Events",   "Weighted counters", 10, -0.5, 9.5,   3, theWeight);
  
  genZZ = DiBosonParticle(Z0, Z1);
  eventGenaftercut++;
  weightGenaftercut += theWeight;
  
  

  // ~~~~~~~~~~~~~~~~~~~~~~~~~~ Histograms ~~~~~~~~~~~~~~~~~~~~~~~~~
  // Histograms after cuts (JJmass is in helper histos)
  theHistograms.fill("Gen4l_mass_AC",   "m_{4l}",   300, 0, 3000, llllp4.M(), theWeight);
  theHistograms.fill("GenZZJJ_mass_AC", "m_{ZZjj}", 500, 0, 5000, ZZjjp4.M(), theWeight);
  
  helper_->PlotBoson(Z0, "GenZ0", theWeight, "AC");
  helper_->PlotBoson(Z1, "GenZ1", theWeight, "AC");
  
  helper_->PlotDiBoson(DiBosonParticle(Z0, Z1), "GenZZ", theWeight, "AC");
  
  helper_->PlotJets(Jet0, Jet1, "Gen", theWeight, "AC");


    
  // ~~~~~~~~~~~~~~~~~~~~~ VBS enriched regions ~~~~~~~~~~~~~~~~~~~~
  // -------------------------- Loose VBS --------------------------
  if(jjp4.M() > 400.){
    theHistograms.fill("Gen4l_mass_AC_Loose",   "m_{4l}",   300, 0, 3000, llllp4.M(), theWeight);
    theHistograms.fill("GenJJ_mass_AC_Loose",   "m_{jj}",   400, 0, 4000, jjp4.M(),   theWeight);
    theHistograms.fill("GenZZJJ_mass_AC_Loose", "m_{ZZjj}", 500, 0, 5000, ZZjjp4.M(), theWeight);
  
    cut++;  
    theHistograms.fill("GenCuts_wei", "Weighted counters", 10, -0.5, 9.5, cut, theWeight);
  }

  
  // -------------------------- Tight VBS -------------------------- 
  if(jjp4.M() > 1000.){
    theHistograms.fill("Gen4l_mass_AC_Tight",   "m_{4l}",   300, 0, 3000, llllp4.M(), theWeight);
    theHistograms.fill("GenJJ_mass_AC_Tight",   "m_{jj}",   400, 0, 4000, jjp4.M(),   theWeight);
    theHistograms.fill("GenZZJJ_mass_AC_Tight", "m_{ZZjj}", 500, 0, 5000, ZZjjp4.M(), theWeight);
  
    cut++;  
    theHistograms.fill("GenCuts_wei", "Weighted counters", 10, -0.5, 9.5, cut, theWeight);
  }

  

  // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  // ~~~~~~~~~~~~~~~~~~~~~ End of gen Analysis ~~~~~~~~~~~~~~~~~~~~~
  // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
}


void ZZjjAnalyzer::RecoAnalysis(DiBosonLepton &recoZZ, Particle &Jet0, Particle &Jet1){
  // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  // ~~~~~~~~~~~~~~~~~~~~ Begin of reco Analysis ~~~~~~~~~~~~~~~~~~~
  // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  int cut = 0;
  theHistograms.fill("RecoCut",     "Reco events after cuts", 14, -0.5, 13.5, cut);
  theHistograms.fill("RecoCut_wei", "Reco events after cuts", 14, -0.5, 13.5, cut, theWeight);

  
  
  // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Jets ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  // ----- jets' number > 2 and |eta| < 4.7
  helper_->FindLeadingJets(jets, Jet0, Jet1);
  
  if(Jet0.p4().Mt() == 0){
    return;
  }
  
  cut++;
  theHistograms.fill("RecoCut",     "Events after cuts", 13, -0.5, 12.5, cut);
  theHistograms.fill("RecoCut_wei", "Events after cuts", 14, -0.5, 13.5, cut, theWeight);

  

  // ~~~~~~~~~~~~~~~~~~~~~~~ ZZ candidates ~~~~~~~~~~~~~~~~~~~~~~~~~
  if(ZZ->mass() == 0){
    return;
  }
  
  cut++;
  theHistograms.fill("RecoCut",     "Events after cuts", 13, -0.5, 12.5, cut);
  theHistograms.fill("RecoCut_wei", "Events after cuts", 14, -0.5, 13.5, cut, theWeight);
  theHistograms.fill("ZZ_Events",   "Weighted counters", 10, -0.5, 9.5,    4, theWeight);
  
  eventReco++;
  weightReco += theWeight;
  
  
  
  // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ All ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  TLorentzVector llllp4 = ZZ->first().daughter(0).p4() + ZZ->first().daughter(1).p4() + ZZ->second().daughter(0).p4() + ZZ->second().daughter(1).p4();
  TLorentzVector jjp4 = Jet0.p4() + Jet1.p4();
  TLorentzVector ZZjjp4 = ZZ->first().p4() + ZZ->first().p4() + Jet0.p4() + Jet1.p4();


  // Histograms before cuts
  theHistograms.fill("Reco4l_mass_BC",   "m_{4l}",   300, 0, 3000, llllp4.M(), theWeight);
  theHistograms.fill("RecoJJ_mass_BC",   "m_{jj}",   400, 0, 4000, jjp4.M(),   theWeight);
  theHistograms.fill("RecoZZJJ_mass_BC", "m_{ZZjj}", 500, 0, 5000, ZZjjp4.M(), theWeight);
  
  

  // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Cuts ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  if(llllp4.M() < 180.)
    return;
  cut++;
  theHistograms.fill("RecoCut",     "Events after cuts", 13, -0.5, 12.5, cut);
  theHistograms.fill("RecoCut_wei", "Events after cuts", 14, -0.5, 13.5, cut, theWeight);

  if(jjp4.M() < 100.)
    return;
  cut++;
  theHistograms.fill("RecoCut",     "Events after cuts", 13, -0.5, 12.5, cut);
  theHistograms.fill("RecoCut_wei", "Events after cuts", 14, -0.5, 13.5, cut, theWeight);

  
  // Histograms after cuts (JJ mass and deltaEta in helper histograms)
  helper_->PlotDiBoson(*ZZ, "ZZ", theWeight, "AC");
  
  helper_->PlotJets(Jet0, Jet1, "Reco", theWeight, "AC");
  
  theHistograms.fill("Reco4l_mass_AC",   "m_{4l}",   300, 0, 3000, llllp4.M(), theWeight);
  theHistograms.fill("RecoZZjj_mass_AC", "m_{ZZjj}", 500, 0, 5000, ZZjjp4.M(), theWeight);



  // ~~~~~~~~~~~~~~~~~~~~~~ MELA Discriminants ~~~~~~~~~~~~~~~~~~~~~
  theHistograms.fill("MELA_JJVBF", "JJVBF Nominal", 202, -2, 40, mela->JJVBF_Nominal(), theWeight);
  theHistograms.fill("MELA_JJQCD", "JJQCD Nominal", 202, -2, 40, mela->JJQCD_Nominal(), theWeight);
  theHistograms.fill("MELA_JJEW",  "JJEW Nominal",  202, -2, 40, mela->JJEW_Nominal(),  theWeight);
  
  TFile* f_ = TFile::Open("../../ZZAnalysis/AnalysisStep/data/cconstants/SmoothKDConstant_m4l_DjjVBF13TeV.root");
  TSpline3* ts = (TSpline3*)(f_->Get("sp_gr_varReco_Constant_Smooth")->Clone());
  f_->Close();
  
  float c = 8.5*ts->Eval(ZZ->mass());  
  float disc_VBF_QCD = mela->JJVBF_Nominal()/(mela->JJVBF_Nominal() + c*mela->JJQCD_Nominal());
  
  theHistograms.fill("MELA_discriminant", "Discriminant #frac{JJVBF}{c JJQCD + JJVBF}", 50, 0, 1, disc_VBF_QCD, theWeight);
  
  
    
  recoZZ = DiBosonLepton(*ZZ);
  eventRecoaftercut++;
  weightRecoaftercut += theWeight;
  theHistograms.fill("ZZ_Events", "Weighted counters", 10, -0.5, 9.5, 5, theWeight);
  
  

  // ~~~~~~~~~~~~~~~~~~~~~~~~~~ Histograms ~~~~~~~~~~~~~~~~~~~~~~~~~
  helper_->PlotDiBoson(recoZZ, "RecoZZ", theWeight, "AC");
  helper_->PlotJets(Jet0, Jet1, "Reco", theWeight, "AC");


  
  // ~~~~~~~~~~~~~~~~~~~~~ VBS enriched regions ~~~~~~~~~~~~~~~~~~~~
  // -------------------------- Loose VBS --------------------------
    cut++;
  if(jjp4.M() > 400.){
    theHistograms.fill("Reco4l_mass_AC_Loose",   "m_{4l}",   300, 0, 3000, llllp4.M(), theWeight);
    theHistograms.fill("Recojj_mass_AC_Loose",   "m_{jj}",   400, 0, 4000, jjp4.M(),   theWeight);
    theHistograms.fill("RecoZZjj_mass_AC_Loose", "m_{ZZjj}", 500, 0, 5000, ZZjjp4.M(), theWeight);
    theHistograms.fill("RecoCut", "Events after cuts", 13, -0.5, 12.5, cut);
    theHistograms.fill("RecoCut_wei", "Events after cuts", 14, -0.5, 13.5, cut, theWeight);
  }

  
  // -------------------------- Tight VBS --------------------------
    cut++; 
  if(jjp4.M() > 1000.){
    theHistograms.fill("Reco4l_mass_AC_Tight",   "m_{4l}",   300, 0, 3000, llllp4.M(), theWeight);
    theHistograms.fill("Recojj_mass_AC_Tight",   "m_{jj}",   400, 0, 4000, jjp4.M(),   theWeight);
    theHistograms.fill("RecoZZjj_mass_AC_Tight", "m_{ZZjj}", 500, 0, 5000, ZZjjp4.M(), theWeight);
    theHistograms.fill("RecoCut", "Events after cuts", 13, -0.5, 12.5, cut);
    theHistograms.fill("RecoCut_wei", "Events after cuts", 14, -0.5, 13.5, cut, theWeight);
  }


  // ---------------------------------------------------------------
    cut++;
  if(disc_VBF_QCD < 0.7){
    theHistograms.fill("RecoCut", "Events after cuts", 13, -0.5, 12.5, cut);
    theHistograms.fill("RecoCut_wei", "Events after cuts", 14, -0.5, 13.5, cut, theWeight);
    theHistograms.fill("Reco4l_mass_AM",   "m_{4l}",   300, 0, 3000, llllp4.M(), theWeight);
    theHistograms.fill("RecoZZjj_mass_AM", "m_{ZZjj}", 500, 0, 5000, ZZjjp4.M(), theWeight);
    helper_->PlotDiBoson(*ZZ, "RecoZZ", theWeight, "AM");
  }

  
  
  // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  // ~~~~~~~~~~~~~~~~~~~~~ End of reco Analysis ~~~~~~~~~~~~~~~~~~~~
  // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
}


void ZZjjAnalyzer::GenRecoAnalysis(const DiBosonParticle genZZ, const Particle genJet0, const Particle genJet1, const DiBosonLepton recoZZ, const Particle recoJet0, const Particle recoJet1){
  // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  // ~~~~~~~~~~~~~~~~~~~~ Begin of Reco vs Gen ~~~~~~~~~~~~~~~~~~~~~
  // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  eventGenReco++;
  weightGenReco += theWeight;
  theHistograms.fill("ZZ_Events", "Weighted counters", 10, -0.5, 9.5, 6, theWeight);

  
  // check if gen and reco IDs are the same
  int genZZID = abs(genZZ.first().daughter(0).id()) + abs(genZZ.second().daughter(0).id());
  int recoZZID = abs(recoZZ.first().daughter(0).id()) + abs(recoZZ.second().daughter(0).id());


  // check if jets are the same
  bool GRJet0 = genJet0.pt() - recoJet0.pt() < 0.1;
  bool GRJet1 = genJet1.pt() - recoJet1.pt() < 0.1;
  theHistograms.fill("GR_JJ_arethesame", "GenJ0 and RecoJ0 (0,1), GenJ1 and RecoJ1 (2,3) are the same", 6, -1.5, 4.5, GRJet0);
  theHistograms.fill("GR_JJ_arethesame", "GenJ0 and RecoJ0 (0,1), GenJ1 and RecoJ1 (2,3) are the same", 6, -1.5, 4.5, GRJet1 + 2);

  
  // ~~~~~~~~~~~~~~~~~~~~~~~~~~ Histograms ~~~~~~~~~~~~~~~~~~~~~~~~~
  // genZ vs recoZ ID
  theHistograms.fill("GR_ID_genZZ_vs_recoZZ", "GenZZ's and RecoZZ's first daughters ID", 5, 19, 29, 5, 19, 29, genZZID, recoZZID);
  if(genZZID != recoZZID)
    theHistograms.fill("GR_ID_genZZ_vs_recoZZ_2", "GenZZ's and RecoZZ's first daughters ID", 5, 19, 29, 5, 19, 29, genZZID, recoZZID);

  
  
  // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  // ~~~~~~~~~~~~~~~~~~~~~ End of Reco vs Gen ~~~~~~~~~~~~~~~~~~~~~~
  // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
}


void ZZjjAnalyzer::GenNoRecoAnalysis(){
  // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  // ~~~~~~~~~~~~~~~~~~~~ Begin of Gen NO Reco ~~~~~~~~~~~~~~~~~~~~~
  // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  eventGenNOReco++;
  weightGenNOReco += theWeight;
  theHistograms.fill("ZZ_Events", "Weighted counters", 10, -0.5, 9.5, 7, theWeight);

  
  // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Jets ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  Particle recoJet0;
  Particle recoJet1;
  helper_->FindLeadingJets(jets, recoJet0, recoJet1);

  if(recoJet0.pt() == 0)
    theHistograms.fill("GNR_Events", "Gen not Reco events", 10, -0.5, 9.5, 0, theWeight);
  
  if(recoJet1.pt() == 0)
    theHistograms.fill("GNR_Events", "Gen not Reco events", 10, -0.5, 9.5, 1, theWeight);

  

  // ~~~~~~~~~~~~~~~~~~~~~~~ ZZ candidates ~~~~~~~~~~~~~~~~~~~~~~~~~
  if(ZZ->pt() == 0)
    theHistograms.fill("GNR_Events", "Gen not Reco events", 10, -0.5, 9.5, 2, theWeight);

  // look for where ZZ is filled
  
  
  
  // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  // ~~~~~~~~~~~~~~~~~~~~~ End of Gen NO Reco ~~~~~~~~~~~~~~~~~~~~~~
  // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
}


void ZZjjAnalyzer::begin(){
  // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  // ~~~~~~~~~~~~~~~~~~~~~~~ Begin function ~~~~~~~~~~~~~~~~~~~~~~~~
  // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  
  cout << "\n--------------------------------------------------------------------------" << endl;
  cout << "\n Starting ZZAnalysis on sample in " << fileName << endl;

  // event counters
  eventSample = 0;
  eventGen = 0;
  eventReco = 0;
  eventGenReco = 0;
  eventGenaftercut = 0;
  eventRecoaftercut = 0;
  eventGenNOReco = 0;
  eventRecoNOGen = 0;
  
  // free counters
  
  
  // time begins
  begintime = ((float)clock())/CLOCKS_PER_SEC;
}


void ZZjjAnalyzer::analyze(){
  // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  // ~~~~~~~~~~~~~~~~~~~~~~~~ Main analysis ~~~~~~~~~~~~~~~~~~~~~~~~
  // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  
  eventSample++;
  theHistograms.fill("ZZ_Events", "Weighted counters", 10, -0.5, 9.5, 1, theWeight);
  
  //gen variables	
  DiBosonParticle genZZ;
  Particle genJet0;
  Particle genJet1;
  
  //reco variables
  DiBosonLepton recoZZ;
  Particle recoJet0;
  Particle recoJet1;
  
  //Gen analysis
  ZZjjAnalyzer::GenAnalysis(genZZ, genJet0, genJet1);
  
  //Reco analysis
  //if(genZZ.pt() != 0){
    ZZjjAnalyzer::RecoAnalysis(recoZZ, recoJet0, recoJet1);
  //}

  ///*
  if(genZZ.pt() != 0. && recoZZ.pt() != 0.){
    //Reco vs Gen analysis
    ZZjjAnalyzer::GenRecoAnalysis(genZZ, genJet0, genJet1, recoZZ, recoJet0, recoJet1);
  }
  //*/

  if(genZZ.pt() != 0. && recoZZ.pt() == 0.){
    //Gen event not reconstructed for reasons
    ZZjjAnalyzer::GenNoRecoAnalysis();    
  }

  if(genZZ.pt() == 0. && recoZZ.pt() != 0.){
    //Reco event not generated for reasons
    eventRecoNOGen++;
    weightRecoNOGen += theWeight;
    theHistograms.fill("ZZ_Events", "Weighted counters", 10, -0.5, 9.5, 8, theWeight);
  }
}


Int_t ZZjjAnalyzer::cut(){  
  return 1;
}


void ZZjjAnalyzer::end(TFile &){
  // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  // ~~~~~~~~~~~~~~~~~~~~~~~~~ End banner ~~~~~~~~~~~~~~~~~~~~~~~~~~
  // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  cout << "\n--------------------------------------------------------------------------" << endl;
  
  cout << "\nEvents of the sample analyzed:                       " << setw(9) << eventSample << endl;
  cout << "Gen events analyzed:                                 " << setw(9) << eventGen << "\t" << eventGen*100./eventSample << "%" << "\t" << eventGen*weightGen*1.0 << endl;
  cout << "Gen events after all cuts:                           " << setw(9) << eventGenaftercut << "\t" << eventGenaftercut*100./eventGen << "%" << "\t" << eventGenaftercut*weightGenaftercut*1.0 << endl;
  cout << "Reco events analyzed:                                " << setw(9) << eventReco << "\t" << eventReco*100./eventSample << "%" << "\t" << eventReco*weightReco*1.0 << endl;
  cout << "Reco events after all cuts:                          " << setw(9) << eventRecoaftercut << "\t" << eventRecoaftercut*100./eventReco << "%" << "\t" << eventRecoaftercut*weightRecoaftercut*1.0 << endl;
  cout << "Gen&Reco events analyzed:                            " << setw(9) << eventGenReco << "\t" << eventGenReco*100./eventSample << "%" << "\t" << eventGenReco*weightGenReco*1.0 << endl;
  cout << "Gen events not reconstructed:                        " << setw(9) << eventGenNOReco << "\t" << eventGenNOReco*100./eventGen << "%" << "\t" << eventGenNOReco*weightGenNOReco*1.0 << endl;
  cout << "Reco events that weren't gen events:                 " << setw(9) << eventRecoNOGen << "\t" << eventRecoNOGen*100./eventReco << "%"<< "\t" << eventRecoNOGen*weightRecoNOGen*1.0  << endl;
    
  // execution time
  endtime = ((float)clock())/CLOCKS_PER_SEC;
  helper_->printTime(begintime, endtime);
  cout << "\n--------------------------------------------------------------------------" << endl;
}
