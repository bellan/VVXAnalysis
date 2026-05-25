#include "VVXAnalysis/TreeAnalysis/interface/VZGAnalyzer.h"
#include "VVXAnalysis/Commons/interface/SignalDefinitions.h"
#include "VVXAnalysis/Commons/interface/Utilities.h"
#include "VVXAnalysis/Commons/interface/Comparators.h"

#include "VVXAnalysis/Commons/interface/GenVBHelper.h"

#include <TFile.h>
#include <TMVA/Reader.h>

#include "correction.h"

#include <fstream>
#include <stdexcept>

#include <boost/foreach.hpp>
#define foreach BOOST_FOREACH

// #include <TString.h>

#include <boost/assign/std/vector.hpp>

#include <functional>
#include <map>

using namespace boost::assign;

using std::cout;
using std::endl;

using namespace phys;
bool UNBLIND=true;
bool isRunForCR=false;
std::string REGION_FOR_FIT = "CR2P_1VL";
int ANALYSIS_CUTs_WP = 2;
int BDT_CUT= 0.95;
double ALPHAS_CAP=1.;
bool verboseControlBlinding= false;
double DFQG_RelVar=0.05;
    
const std::vector<double> binEdges =  {-1.00,-0.85,-0.70,-0.56,-0.43,-0.31,-0.19,-0.07,0.04,0.14,0.24,0.33,0.42,0.51,0.59,0.67,0.74,0.81,0.88,0.94,1.00}; //CT: used for PhD thesis
//const std::vector<double> binEdges =  {-1.00,-0.85,-0.70,-0.56,-0.43,-0.31,-0.19,-0.07,0.05,0.16,0.27,0.38,0.49,0.6,0.7,0.8,0.9,1.00}; //CT: rebin5

//const std::vector<double>   binEdges =  {-1.00,-0.90,-0.75,-0.55,      -0.30,           0.,            0.30,         0.55,       0.75,      0.90,    1.00}; //CT: distributed: rebin4

//const std::vector<double> binEdges =  {-0.31,-0.19,-0.07,0.04,0.14,0.24,0.33,0.42,0.51,0.59,0.67,0.74,0.81,0.88,0.94,1.00}; //CT: for testing constrained unc 
//const std::vector<double> binEdges =  {-1.00, 0.00,0.12,0.23,0.33,0.42,0.52,0.63,0.75,0.875,1.00}; //CT: less aggressive binning: rebin1
//const std::vector<double> binEdges =  {-1.00, -0.4, -0.1, 0.12,0.23,0.33,0.42,0.52,0.63,0.75,0.875,1.00}; //CT: less aggressive binning: rebin2
//const std::vector<double> binEdges =    {-1.00,-0.85,-0.70, -0.50,       -0.25,         0.,         0.25,     0.45,     0.60,    0.74,      0.90,    1.00}; //CT: distributed: rebin3

const std::vector<double> mll_binEdges =  {80., 85., 87., 89., 91., 93., 95., 100., 110., 120.}; 

const std::vector<double> rewgtBinEdges =  {0,10,20,30,40,50,60,70,80,90,100,110,120,130,140,150,160,170,180,190,200,220,240,260,280,300,330,360,400,480,650};

const std::vector<double> DYrewgtBinEdges =  {0,10,20,30,40,50,60,70,80,90,100,110,120,130,140,150,160,170,180,195,210,235,250,270,290,320,350,380,420,490,650};

const std::vector<double> ZGrewgtBinEdges =  {0,10,20,30,40,50,60,70,80,90,100,110,120,130,140,150,160,170,180,400};
/*
const std::vector<double> DYreweights_2016preVFP  = { 2.00917, 1.79445, 1.36992, 1.46493, 1.2469, 1.19764, 1.37287, 1.35927, 1.34544, 1.34675, 1.35912, 1.31268, 1.15761, 1.21491, 1.5606, 1.46089, 1.47166, 1.16549, 1.27846, 1.76397, 1.52727, 1.32305, 1.17898, 2.16624, 1.58863, 0.870507, 1.29542, 1.62954, 1.17145, 3.00147 };
const std::vector<double> DYreweights_2016postVFP = { 1.30377, 1.47856, 1.35739, 1.31193, 1.27247, 1.44967, 1.60628, 1.56198, 1.17672, 1.53046, 1.45697, 1.3789, 1.17725, 1.1753, 1.4612, 1.44108, 1.51864, 1.43546, 1.16643, 1.25587, 1.3805, 1.34821, 0.988318, 1.54942, 1.09622, 1.22009, 2.53514, 2.17051, 1.56687, 0.799289  };
const std::vector<double> DYreweights_2017 = { 1.08715, 1.18149, 1.2471, 1.19377, 1.19708, 1.18842, 1.14534, 1.11701, 1.21096, 1.28702, 1.50582, 1.56286, 1.25695, 1.11851, 1.06644, 1.11402, 1.25952, 1.19327, 1.41701, 1.15879, 1.37471, 1.09833, 1.13833, 1.47355, 1.15392, 1.12907, 1.08345, 0.84334, 1.06811, 2.64004 };
const std::vector<double> DYreweights_2018 = { 1.01513, 1.27004, 1.18847, 1.15435, 1.11292, 1.20551, 1.10328, 1.10396, 1.19841, 1.14451, 1.08021, 1.17151, 1.1553, 1.20533, 1.288, 1.09464, 1.14136, 1.07776, 1.12247, 1.35225, 1.05515, 1.2525, 1.22472, 0.92048, 1.27397, 1.28208, 1.44384, 0.931377, 0.97728, 2.09973 };
*/
const std::vector<double> DYreweights_Run2 = { 1.1425, 1.30874, 1.24618, 1.21632, 1.17348, 1.22305, 1.19032, 1.17641, 1.21754, 1.24948, 1.27542, 1.32013, 1.19056, 1.17448, 1.25021, 1.17724, 1.25261, 1.16321, 1.23544, 1.30837, 1.23083, 1.21792, 1.16034, 1.2197, 1.24126, 1.15815, 1.37485, 1.03262, 1.09324, 1.92615 };

const std::vector<double> DYrewErr_Run2 = { 0.0567861, 0.0530675, 0.0418608, 0.0360185, 0.0328456, 0.0352394, 0.0350258, 0.0355526, 0.0401509, 0.0448921, 0.0492029, 0.0569995, 0.0509339, 0.0533836, 0.0636173, 0.0629042, 0.0753955, 0.0719311, 0.0706459, 0.0917001, 0.0726525, 0.109921, 0.0968275, 0.125411, 0.119484, 0.131386, 0.216469, 0.143534, 0.156537, 0.559408 };
/*
const std::vector<double> ZGreweights_2016preVFP  = { 1.30531, 1.11282, 0.977, 0.970804, 1.21571, 1.48782, 1.27616, 1.27632, 1.30883, 1.31963, 2.02662, 1.3168, 1.32528, 1.2418, 0.830399, 0.783573, -0.0931671, 2.80097, 1.23566};
const std::vector<double> ZGreweights_2016postVFP = { 0.917345, 1.15899, 0.977712, 1.19733, 1.1256, 1.38638, 1.18733, 1.46917, 1.50516, 1.59668, 0.906755, 1.33604, 1.75536, 0.284094, 1.2815, 0.32339, 1.83554, 2.23789, 1.05929  };
const std::vector<double> ZGreweights_2017 = { 0.852612, 0.938297, 0.954081, 1.00565, 0.932689, 1.03828, 1.13545, 0.8911, 1.49594, 1.193, 0.846078, 1.29195, 1.04021, 1.10745, 0.965414, 0.890614, 0.894638, -1.3284, 1.89145 };
const std::vector<double> ZGreweights_2018 = { 0.62344, 0.865009, 0.9259, 0.913941, 0.986357, 1.31534, 1.26514, 1.10043, 1.52301, 1.34909, 1.17461, 1.37038, 1.21519, 1.0367, 0.995551, 1.22477, 1.35664, 1.24588, 0.950638 };
*/
const std::vector<double> ZGreweights_Run2 = { 0.823679, 0.954957, 0.947289, 0.982136, 1.01551, 1.25815, 1.21801, 1.10996, 1.48301, 1.32895, 1.14251, 1.33304, 1.2409, 0.993592, 1.00056, 0.969438, 1.09989, 0.927254, 1.24255 };

const std::vector<double> ZGrewErr_Run2 = { 0.088428, 0.0615822, 0.0524166, 0.0583568, 0.0667823, 0.0797848, 0.0841933, 0.0975099, 0.126066, 0.138374, 0.139703, 0.177405, 0.178318, 0.172954, 0.217512, 0.262971, 0.289531, 0.464014, 0.188374 };


bool IsARunForMVAFeat=false;
bool verbose = false;
bool fullPlotList = false;
bool fiducial_run =true;
bool runningSR2PFJ=false;
double etacut=4.7;
double ptcut=30;
double LumiSF=1.0;//7.035--> 2016post||3.32-->2017||2.29-->2018||8.16--> 2016post||7.035--> 2016pre
//FOR DY PART1
//double LumiSF=7.42;
//double weight=theWeight*LumiSF;
double PhEffSF=1.;
double PhEffSFUnc=0.;

double rewgt=1.;
double rewErr=0.;
double TF=1.;
double TF_uncAbs = 0.18;
double TF_16preVFP  = 1.391;
double TF_16postVFP = 1.395;
double TF_17        = 1.243;
double TF_18        = 1.190;

double subFactor=1.;
double subFactor_16preVFP  = 0.8998324022;
double subFactor_16postVFP = 0.904898959;
double subFactor_17        = 0.9033030263;
double subFactor_18        = 0.8982220908;

double dR_jetRatio_cut = 0.4;
double dR_FJRatio_cut = 0.8;
int cutsToApply=4;

TString BDTmodelPath = "bdtModels/ZZGasOnlySig_VLRetunedAlt__BDT_Xgrad_d3_N050.weights.xml"; // ZZG semilep single sig, WZG semilep as bkg
//TString BDTmodelPath = "bdtModels/postBkgEnriching_kin__BDT_Xgrad_d3_N100.weights.xml";   //vanilla model after retraining SR2P_1k w/bkg enriched in CR2P_1k 
//TString BDTmodelPath = "bdtModels/postQvG_VLRetunedAlt__BDT_Xgrad_d3_N100.weights.xml"; //after QvG SFs
//TString BDTmodelPath = "bdtModels/VLRetunedAlt_noNegWgt_BDT_Xgrad_d3_N030.weights.xml"; //used for PhDthesis

double CRZOFFMassCut = 80;
double CRZONMassCut  = 85;
double CRZinfMassCut = 80;
double CRZsupMassCut = 100;

static TMVA::Reader* reader = nullptr;

static float feat_dRLG, feat_recoVMass, feat_mllG, feat_PhMVAId, feat_J1Girth, feat_ptJ1, feat_ptjj, feat_HT, feat_J0Girth, feat_deltaR_J0Gamma, feat_deltaR_J1Gamma, feat_ptGamma, feat_dPhiZG, feat_j0_btagDeepFlavQG, feat_j1_btagDeepFlavQG;//, feat_J0Puds, feat_J1PG;

std::shared_ptr<correction::CorrectionSet> cset;

void loadQGCorrections(const std::string &jsonFile) {
  cset = correction::CorrectionSet::from_file(jsonFile);
  std::cout << "Loaded QvG SFs from: " << jsonFile << std::endl;
}

struct QGScaleFactor {
  double central;
  double up;
  double down;
};

// wp = "L", "M" o "T"
QGScaleFactor getQGSF(bool isMC, double qvg_score, double eta, double pt, int absFlav, const std::string &wp) {
  QGScaleFactor sf{1.0, 1.0, 1.0};  // default

  if(isMC){
    try {
      auto corr = cset->at("deepJet_shape");  // standard name in the json
      sf.central = corr->evaluate({"central", wp, qvg_score, absFlav, eta, pt});
      sf.up      = corr->evaluate({"up",      wp, qvg_score, absFlav, eta, pt});
      sf.down    = corr->evaluate({"down",    wp, qvg_score, absFlav, eta, pt});
    } catch (std::exception &e) {
      std::cerr << "Error evaluating QG SF: " << e.what() << std::endl;
    }
  }
  return sf;
}

static std::vector<Systematic> jetSysts;
Systematic noJERCsys = {SystType::Nominal, 0, "", ""};

std::string systName(Systematic syst){
  if(syst.type==SystType::Nominal || syst.direction==0){
    std::cout<<"CODING ERROR: systName function wrongly called"<<endl;
    return "";
  }
  std::string jercString = (syst.type==SystType::JES) ? "jes"       : "jer" ;
  std::string srcString  = (syst.type==SystType::JES) ? syst.source : ""    ;
  std::string dirString  = (syst.direction>0)         ? "Up"        : "Down";
  
  return jercString+srcString+syst.yearSys+"_"+dirString;
  
}

using JesFunc = std::function<double(const Jet&)>;

std::map<std::string, JesFunc> jesUncMap = {
    {"Abs",             [](const Jet& j){ return j.jesUnc().Abs           ; }},
    {"BBEC1",           [](const Jet& j){ return j.jesUnc().BBEC1         ; }},
    {"EC2",             [](const Jet& j){ return j.jesUnc().EC2           ; }},
    {"FlavQCD",         [](const Jet& j){ return j.jesUnc().FlavQCD       ; }},
    {"HF",              [](const Jet& j){ return j.jesUnc().HF            ; }},
    {"RelBal",          [](const Jet& j){ return j.jesUnc().RelBal        ; }},

    {"AbsYear",         [](const Jet& j){ return j.jesUnc().Abs_year      ; }},
    {"BBEC1Year",       [](const Jet& j){ return j.jesUnc().BBEC1_year    ; }},
    {"EC2Year",         [](const Jet& j){ return j.jesUnc().EC2_year      ; }},
    {"HFYear",          [](const Jet& j){ return j.jesUnc().HF_year       ; }},
    {"RelSampleYear",   [](const Jet& j){ return j.jesUnc().RelSample_year; }}

};

void ScaleP4(TLorentzVector& p4, double ptScaled, double ptCentral){
  p4.Pz()*ptScaled/ptCentral;
  double px=p4.Px()*ptScaled/ptCentral;
  double py=p4.Py()*ptScaled/ptCentral;
  double pz=p4.Pz()*ptScaled/ptCentral;
  double p0=p4.E() *ptScaled/ptCentral;//TMath::Sqrt(px*px+py*py+pz*pz+p4.M2());
  p4.SetPxPyPzE(px, py, pz, p0);
  
}

double Sigmoid(double x, double x0, double A){
  return 1./(1. + TMath::Exp(-A*(x-x0) ) );
}

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

bool KinematicsOKafterSys(phys::Particle p, float ptVaried, float ptThr,float eta)
{
  if (fabs(p.eta()) < eta && fabs(ptVaried) > ptThr) return true;
  return false;
}

void VZGAnalyzer::begin(){
  
  reader = new TMVA::Reader("Color:Silent");

  reader->AddVariable("ptJ1", &feat_ptJ1);
  reader->AddVariable("ptjj", &feat_ptjj);
  reader->AddVariable("ptGamma",&feat_ptGamma);
  reader->AddVariable("mllPh",&feat_mllG);
  reader->AddVariable("deltaR_J0Gamma",&feat_deltaR_J0Gamma);
  reader->AddVariable("deltaR_J1Gamma",&feat_deltaR_J1Gamma);
  reader->AddVariable("dPhiZG",&feat_dPhiZG);
  reader->AddVariable("recoVMass",&feat_recoVMass);
  reader->AddVariable("dRLG",&feat_dRLG);
  reader->AddVariable("PhMVAId",  &feat_PhMVAId);
  reader->AddVariable("J0Girth",  &feat_J0Girth);
  reader->AddVariable("J1Girth",  &feat_J1Girth);
  //  reader->AddVariable("J1DeepProb_g",  &feat_J1PG);
  //  reader->AddVariable("J0DeepProb_uds",  &feat_J0Puds);
  reader->AddVariable("HT",&feat_HT);

  //redefined features
  reader->AddVariable("J0DeepProb_uds/(J0DeepProb_uds + J0DeepProb_g)",  &feat_j0_btagDeepFlavQG);
  reader->AddVariable("J1DeepProb_uds/(J1DeepProb_uds + J1DeepProb_g)",  &feat_j1_btagDeepFlavQG);
  
  reader->BookMVA("BDT", BDTmodelPath);

  //  std::shared_ptr<correction::CorrectionSet> DFcorrSet;


  
  /*
  cout<<'\n';
  for(char i=0; i<25; ++i) cout<<'-';
  cout<<" Start of VZGAnalyzer ";
  for(char i=0; i<25; ++i) cout<<'-';
  cout<<'\n';
  */
  //  const char* CMSSW_BASE = getenv("CMSSW_BASE");
  std::string VVXAnalysis_dir = "../";

  std::string year_str = std::to_string(year);
  if(year == 2016)      year_str += subEra_;

  //----BLOCK ASSIGNING DY FLAT TF PER YEAR----------//
  if(theSampleInfo.isMC() && theSampleInfo.fileName().find("DY")!=std::string::npos){
    if(year==2016) TF = TF_16preVFP;
    if(year==2017) TF = TF_17;
    if(year==2018) TF = TF_18;
  }
  TF=1.;
  //-------------------------------------------------//
    
  if(!theSampleInfo.isMC()){
    if(year==2016){
      if(year_str.find("preVFP")!=std::string::npos) subFactor = subFactor_16preVFP;
      else subFactor = subFactor_16postVFP;
    }
    if(year==2017) subFactor = subFactor_17;
    if(year==2018) subFactor = subFactor_18;
  }
  subFactor=1.;

  //CT back here
  /*
  std::string year_str;
  if(PHFR_SPLIT){
    year_str = std::to_string(year);
    if(year == 2016)
      year_str += subEra_;
  }
  else{
    year_str = "Run2";
  }
  */
		  
  // Photon efficiency SF for cut-based ID (temporary)
  std::string pathPhotonEffSF(Form("%s/Commons/data/egammaEffi.txt_EGM2D_Pho_Loose_UL%d%s.root", VVXAnalysis_dir.c_str(), year%100, (subEra_.size() > 0 ? ('_'+subEra_).c_str() : "")));
  hPhotonEffSF_           = getHistfromFile(pathPhotonEffSF.c_str(), "EGamma_SF2D");
  hPhotonEffSF_maxPt_     = hPhotonEffSF_->GetYaxis()->GetBinUpEdge(hPhotonEffSF_->GetNbinsY());

  // Photon MVA SF
  mapPhotonMVASF_[Photon::MVAwp::wp80] = getHistfromFile(Form("%s/Commons/data/%d_PhotonsMVAwp80.root", VVXAnalysis_dir.c_str(), year), "EGamma_SF2D");
  mapPhotonMVASF_[Photon::MVAwp::wp90] = getHistfromFile(Form("%s/Commons/data/%d_PhotonsMVAwp90.root", VVXAnalysis_dir.c_str(), year), "EGamma_SF2D");
  for(auto& it: mapPhotonMVASF_)  // std::pair<const Photon::MVAwp, std::unique_ptr<TH2F>>
    mapPhotonMVASF_maxPt_[it.first] = it.second->GetYaxis()->GetBinUpEdge(it.second->GetNbinsY());


  //---block loading QvG DF correction set---//
  std::string QvGcorrJsonPath="qvgJsons/"+year_str+"ULqgtagging.json" ;
  //  QvGcorrJsonPath+=year_str ;
  //  QvGcorrJsonPath+="ULqgtagging.json" ;
  
  if(theSampleInfo.isMC()) loadQGCorrections(QvGcorrJsonPath);

  
  if (jetSysts.empty()) {
    jetSysts.push_back({SystType::Nominal, 0, "",""});
    jetSysts.push_back({SystType::JER, +1, "",""});
    jetSysts.push_back({SystType::JER, -1, "",""});

    std::vector<std::string> jesSourcesRun2   = {"Abs", "BBEC1", "EC2", "FlavQCD", "HF", "RelBal"};

    for (const auto& jesSrc : jesSourcesRun2) {
      jetSysts.push_back({SystType::JES, +1, jesSrc, ""});
      jetSysts.push_back({SystType::JES, -1, jesSrc, ""});
    }

    std::vector<std::string> jesSourcesByYear = {"AbsYear", "BBEC1Year", "EC2Year", "HFYear", "RelSampleYear"};

    std::vector<std::string> yearStrings = {"2016", "2017", "2018"};
    
    for (const auto& jesSrcPerYear : jesSourcesByYear) {
      for (const auto& ys : yearStrings) {
	jetSysts.push_back({SystType::JES, +1, jesSrcPerYear, ys});
	jetSysts.push_back({SystType::JES, -1, jesSrcPerYear, ys});
      }

    }
  }
  
  return;
}

double VZGAnalyzer::getJESUncertainty(const Jet& jet, const std::string& source, const std::string& yearUnc) {
  if (yearUnc!="" && yearUnc!=std::to_string(year) )    return 0.;
  return jesUncMap.at(source)(jet);
}

double VZGAnalyzer::getPhotonEffSF_MVA(const phys::Photon& ph, Photon::MVAwp wp) const{
  const TH2F* hSF = mapPhotonMVASF_.at(wp).get();
  float maxPt = mapPhotonMVASF_maxPt_.at(wp);
  Int_t bin = hSF->FindFixBin(
			      ph.eta(),
			      ph.pt() < maxPt ? ph.pt() : maxPt-0.1
			      );
  return hSF->GetBinContent(bin);
}

double VZGAnalyzer::getPhotonEffSFUnc_MVA(const phys::Photon& ph, Photon::MVAwp wp) const{
  const TH2F* hSF = mapPhotonMVASF_.at(wp).get();
  float maxPt = mapPhotonMVASF_maxPt_.at(wp);
  Int_t bin = hSF->FindFixBin(
			      ph.eta(),
			      ph.pt() < maxPt ? ph.pt() : maxPt-0.1
			      );
  return hSF->GetBinError(bin);
}

std::unique_ptr<TH2F> VZGAnalyzer::getHistfromFile(const char* fname, const char* hname, const char* info){
  TFile fileFR(Form(fname, "READ" ));
  if(!fileFR.IsOpen()){
    cout << colour::Red("WARN") << ": file"<<info<<" not found in \""<<fname<<"\"\n";
    return std::unique_ptr<TH2F> (new TH2F(hname, "DEFAULT", 1,0.,1., 1,0.,1.) );
  }

  TObject* retrievedObj = fileFR.Get(hname);
  if(!retrievedObj){
    cout << colour::Red("WARN") << ": histogram \""<<hname<<"\""<<info<<" not found in \""<<fname<<"\"\n"
	 << "\tavailable keys:" ;
    for(auto key : *fileFR.GetListOfKeys()) cout << ' ' << key->GetName() << ',';
    cout << std::endl;
    fileFR.Close();
    return std::unique_ptr<TH2F> (new TH2F(hname, "DEFAULT", 1,0.,1., 1,0.,1.) );
  }

  std::unique_ptr<TH2F> result((TH2F*) retrievedObj);
  result->SetDirectory(nullptr);  // prevent ROOT from deleting it
  cout << "INFO: retrieved histogram \""<<hname<<"\""<<info<<" from \""<<fname<<"\"\n";
  fileFR.Close();
  return result;
}


bool VZGAnalyzer::LeptonicSignalConstraint()
{
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
  return (haveGoodFJCand || haveGoodDiJetCand);

}

bool VZGAnalyzer::PhotonSignalConstraint()
{
       std::vector<phys::Particle> selectedphotons;
       for (auto p : *genParticles)
              if (p.id() == 22 && KinematicsOK(p, 15, 2.4) && p.genStatusFlags().test(phys::isPrompt) &&  p.genStatusFlags().test(phys::fromHardProcess))
                  selectedphotons.push_back(p);
       if(verbose==true) std::cout<< "Number of selected gen photons = "<<selectedphotons.size()<<std::endl;
       if (selectedphotons.size()>=1) return true;
       return false;
}

bool VZGAnalyzer::IN_GENsignalDef()
{
  if(theSampleInfo.fileName().find("ZH")!=std::string::npos) return true; //CT: by default it's signal
    
  if (LeptonicSignalConstraint() && HadronicSignalConstraint() && PhotonSignalConstraint())
    {
      std::vector<phys::Particle> selectedGENphotons;
      for (auto p : *genParticles)
	if (p.id() == 22 && KinematicsOK(p, 15, 2.4) && p.genStatusFlags().test(phys::isPrompt) &&  p.genStatusFlags().test(phys::fromHardProcess))
	  selectedGENphotons.push_back(p);
      
  
      if (selectedGENphotons.size()>=1)
	{
	  TLorentzVector  GEN_llPh ;
	  double mGEN_llPh;
	  double mGenZ;
	  if(genVBHelper_.ZtoChLep().size()>=1)
	    {
	      GEN_llPh = genVBHelper_.ZtoChLep()[0].p4()+selectedGENphotons.at(0).p4();
	      mGEN_llPh=GEN_llPh.M();
	      mGenZ=genVBHelper_.ZtoChLep()[0].mass();
	      //theHistograms->fill("GEN mjj_vs_mjjG", "GEN mjj_vs_mjjG; mjj [GeV] ; mjj#gamma [GeV]", 35, 50, 120, 30, 50, 350, mGenV, mGEN_jjPh, theWeight*LumiSF);
	      theHistograms->fill("GEN mll_vs_mllG", "GEN mll_vs_mllG; mll [GeV] ; mll#gamma [GeV]", 30, 60, 120, 80, 50, 450, mGenZ, mGEN_llPh, theWeight*LumiSF);
	      if(!fiducial_run)		return true;
	      else if(mGEN_llPh>315-2.5*mGenZ) return true;
	    }
	}
    
    }
  return false;
}

double VZGAnalyzer::VZGMVAScoreEval(phys::Boson<phys::Jet> recoV, phys::Jet recoFJ, std::vector<phys::Photon> selectedphotons, int VBTopo, Systematic jetSystApplied, bool isForDFsys, int DFsysDirection)
{
  double VZGMVAScore=-2.;
  if(!cut(1, recoV, recoFJ, selectedphotons, VBTopo, VZGMVAScore)){
    std::cout<<"BDT building failed"<<endl;
    return -2.;
  }

  TLorentzVector llPh = Z->daughter(0).p4()+Z->daughter(1).p4()+selectedphotons.at(0).p4();
  TLorentzVector l0Ph = Z->daughter(0).p4()+selectedphotons.at(0).p4();
  TLorentzVector l1Ph = Z->daughter(1).p4()+selectedphotons.at(0).p4();

  double mllPh=llPh.M();
  double m2llPh=llPh.M2();
  double m2l0Ph=l0Ph.M2();
  double m2l1Ph=l1Ph.M2();
  //  double mll=Z->mass();
  double mjj=recoV.mass();

  TLorentzVector lljjPh;
  TLorentzVector jjPh;
  double mjjPh=jjPh.M();
  double m2jjPh=jjPh.M2();
  std::vector<TLorentzVector> lljj;

  double mjjScaled_JUncUp, ptjjScaled_JUncUp, ptj0Scaled_JUncUp, ptj1Scaled_JUncUp, mjjScaled_JUncDn, ptjjScaled_JUncDn, ptj0Scaled_JUncDn, ptj1Scaled_JUncDn;
  if (jetSystApplied.type == SystType::JER) {
    ptj0Scaled_JUncUp=recoV.daughter(0).ptJerUp();
    ptj1Scaled_JUncUp=recoV.daughter(1).ptJerUp();
    ptj0Scaled_JUncDn=recoV.daughter(0).ptJerDn();
    ptj1Scaled_JUncDn=recoV.daughter(1).ptJerDn();
  }else if(jetSystApplied.type == SystType::JES){
    ptj0Scaled_JUncUp=recoV.daughter(0).pt()*(1.+getJESUncertainty(recoV.daughter(0), jetSystApplied.source, jetSystApplied.yearSys) ); //recoV.daughter(0).jesUnc().Total );
    ptj1Scaled_JUncUp=recoV.daughter(1).pt()*(1.+getJESUncertainty(recoV.daughter(1), jetSystApplied.source, jetSystApplied.yearSys) ); //recoV.daughter(1).jesUnc().Total );
    ptj0Scaled_JUncDn=recoV.daughter(0).pt()*(1.-getJESUncertainty(recoV.daughter(0), jetSystApplied.source, jetSystApplied.yearSys) ); //recoV.daughter(0).jesUnc().Total );
    ptj1Scaled_JUncDn=recoV.daughter(1).pt()*(1.-getJESUncertainty(recoV.daughter(1), jetSystApplied.source, jetSystApplied.yearSys) ); //recoV.daughter(1).jesUnc().Total );
  }
  if(jetSystApplied.direction<0){
    if(ptj0Scaled_JUncDn < ptj1Scaled_JUncDn) std::swap (ptj0Scaled_JUncDn,ptj1Scaled_JUncDn);
    if(ptj1Scaled_JUncDn < 30) return -2.;
  }else if(jetSystApplied.direction>0){
    if(ptj0Scaled_JUncUp < ptj1Scaled_JUncUp) std::swap (ptj0Scaled_JUncUp,ptj1Scaled_JUncUp);
    if(ptj1Scaled_JUncUp < 30) return -2.;
  }

  

  if(jetSystApplied.direction>0){
    TLorentzVector JUnc_PJ0_scaling_Up = recoV.daughter(0).p4();
    TLorentzVector JUnc_PJ1_scaling_Up = recoV.daughter(1).p4();
    ScaleP4(JUnc_PJ0_scaling_Up, ptj0Scaled_JUncUp, recoV.daughter(0).pt());
    ScaleP4(JUnc_PJ1_scaling_Up, ptj1Scaled_JUncUp, recoV.daughter(1).pt());
    mjjScaled_JUncUp=(JUnc_PJ0_scaling_Up+JUnc_PJ1_scaling_Up).M();
    if(mjjScaled_JUncUp < 50 || mjjScaled_JUncUp > 120) return -2.;

    ptjjScaled_JUncUp=(JUnc_PJ0_scaling_Up+JUnc_PJ1_scaling_Up).Pt();

    lljjPh = Z->daughter(0).p4()+Z->daughter(1).p4()+selectedphotons.at(0).p4()+JUnc_PJ0_scaling_Up+JUnc_PJ1_scaling_Up;
    jjPh = JUnc_PJ0_scaling_Up+JUnc_PJ1_scaling_Up+selectedphotons.at(0).p4();
    mjjPh=jjPh.M();
    m2jjPh=jjPh.M2();

    lljj.push_back(Z->daughter(0).p4());
    lljj.push_back(Z->daughter(1).p4());
    lljj.push_back(JUnc_PJ0_scaling_Up);
    lljj.push_back(JUnc_PJ1_scaling_Up);
    
  }else if(jetSystApplied.direction<0){
    TLorentzVector JUnc_PJ0_scaling_Dn = recoV.daughter(0).p4();
    TLorentzVector JUnc_PJ1_scaling_Dn = recoV.daughter(1).p4();
    ScaleP4(JUnc_PJ0_scaling_Dn, ptj0Scaled_JUncDn, recoV.daughter(0).pt());
    ScaleP4(JUnc_PJ1_scaling_Dn, ptj1Scaled_JUncDn, recoV.daughter(1).pt());
    mjjScaled_JUncDn=(JUnc_PJ0_scaling_Dn+JUnc_PJ1_scaling_Dn).M();
    if(mjjScaled_JUncDn < 50 || mjjScaled_JUncDn > 120) return -2.;
    
    ptjjScaled_JUncDn=(JUnc_PJ0_scaling_Dn+JUnc_PJ1_scaling_Dn).Pt();
    
    lljjPh = Z->daughter(0).p4()+Z->daughter(1).p4()+selectedphotons.at(0).p4()+JUnc_PJ0_scaling_Dn+JUnc_PJ1_scaling_Dn;
    jjPh = JUnc_PJ0_scaling_Dn+JUnc_PJ1_scaling_Dn+selectedphotons.at(0).p4();
    mjjPh=jjPh.M();
    m2jjPh=jjPh.M2();

    lljj.push_back(Z->daughter(0).p4());
    lljj.push_back(Z->daughter(1).p4());
    lljj.push_back(JUnc_PJ0_scaling_Dn);
    lljj.push_back(JUnc_PJ1_scaling_Dn);


  }else if(jetSystApplied.direction==0){
    lljjPh = Z->daughter(0).p4()+Z->daughter(1).p4()+selectedphotons.at(0).p4()+recoV.daughter(0).p4()+recoV.daughter(1).p4();
    jjPh = recoV.daughter(0).p4()+recoV.daughter(1).p4()+selectedphotons.at(0).p4();
    mjjPh=jjPh.M();
    m2jjPh=jjPh.M2();

    lljj.push_back(Z->daughter(0).p4());
    lljj.push_back(Z->daughter(1).p4());
    lljj.push_back(recoV.daughter(0).p4());
    lljj.push_back(recoV.daughter(1).p4());
  }
  /*
  //_____________________________________________BLOCK_TO_CROSS-CHECK_JETScaling___________________________________________//
  if(isJERorJES>0 && isForSysUpDn>0){       //JER UP
    theHistograms->fill("AUX_ptj0_JERup", "AUX_ptj0_JERup" , 40, 0,  400, ptj0Scaled_JUncUp, theWeight*LumiSF);
    theHistograms->fill("AUX_ptj1_JERup", "AUX_ptj1_JERup" , 40, 0,  400, ptj1Scaled_JUncUp, theWeight*LumiSF);
    theHistograms->fill("AUX_ptjj_JERup", "AUX_ptjj_JERup" , 40, 0,  400, ptjjScaled_JUncUp, theWeight*LumiSF);
    theHistograms->fill("AUX_mjj_JERup" , "AUX_mjj_JERup"  , 24, 30, 150, mjjScaled_JUncUp, theWeight*LumiSF);
  }else if(isJERorJES>0 && isForSysUpDn<0){ //JER DN
    theHistograms->fill("AUX_ptj0_JERdn", "AUX_ptj0_JERdn" , 40, 0,  400, ptj0Scaled_JUncDn, theWeight*LumiSF);
    theHistograms->fill("AUX_ptj1_JERdn", "AUX_ptj1_JERdn" , 40, 0,  400, ptj1Scaled_JUncDn, theWeight*LumiSF);
    theHistograms->fill("AUX_ptjj_JERdn", "AUX_ptjj_JERdn" , 40, 0,  400, ptjjScaled_JUncDn, theWeight*LumiSF);
    theHistograms->fill("AUX_mjj_JERdn" , "AUX_mjj_JERdn"  , 24, 30, 150, mjjScaled_JUncDn, theWeight*LumiSF);
  }else if(isJERorJES<0 && isForSysUpDn>0){ //JES UP
    theHistograms->fill("AUX_ptj0_JESup", "AUX_ptj0_JESup" , 40, 0,  400, ptj0Scaled_JUncUp, theWeight*LumiSF);
    theHistograms->fill("AUX_ptj1_JESup", "AUX_ptj1_JESup" , 40, 0,  400, ptj1Scaled_JUncUp, theWeight*LumiSF);
    theHistograms->fill("AUX_ptjj_JESup", "AUX_ptjj_JESup" , 40, 0,  400, ptjjScaled_JUncUp, theWeight*LumiSF);
    theHistograms->fill("AUX_mjj_JESup" , "AUX_mjj_JESup"  , 24, 30, 150, mjjScaled_JUncUp, theWeight*LumiSF);
  }else if(isJERorJES<0 && isForSysUpDn<0){ //JES DN
    theHistograms->fill("AUX_ptj0_JESdn", "AUX_ptj0_JESdn" , 40, 0,  400, ptj0Scaled_JUncDn, theWeight*LumiSF);
    theHistograms->fill("AUX_ptj1_JESdn", "AUX_ptj1_JESdn" , 40, 0,  400, ptj1Scaled_JUncDn, theWeight*LumiSF);
    theHistograms->fill("AUX_ptjj_JESdn", "AUX_ptjj_JESdn" , 40, 0,  400, ptjjScaled_JUncDn, theWeight*LumiSF);
    theHistograms->fill("AUX_mjj_JESdn" , "AUX_mjj_JESdn"  , 24, 30, 150, mjjScaled_JUncDn, theWeight*LumiSF);
  }else if(isJERorJES==0 && isForSysUpDn==0){//central
    theHistograms->fill("AUX_ptj0_central", "AUX_ptj0_central" , 40, 0,  400, recoV.daughter(0).pt(), theWeight*LumiSF);
    theHistograms->fill("AUX_ptj1_central", "AUX_ptj1_central" , 40, 0,  400, recoV.daughter(1).pt(), theWeight*LumiSF);
    theHistograms->fill("AUX_ptjj_central", "AUX_ptjj_central" , 40, 0,  400, recoV.pt(), theWeight*LumiSF);
    theHistograms->fill("AUX_mjj_central" , "AUX_mjj_central"  , 24, 30, 150, mjj, theWeight*LumiSF);
  }
  //______________________________________________________________________________________________________________________//
*/
  
  phys::Photon mostEnergeticPhoton;

  std::stable_sort(selectedphotons.begin(), selectedphotons.end(), phys::EComparator());
  mostEnergeticPhoton = selectedphotons[0];


  std::pair<phys::Photon, phys::Jet> nearestRECOjetstoPhoton;
  std::pair<phys::Photon, phys::Lepton> nearestChLeptToPhoton;

  if(fabs(physmath::deltaR(Z->daughter(0), mostEnergeticPhoton))<fabs(physmath::deltaR(Z->daughter(1), mostEnergeticPhoton)) )
    nearestChLeptToPhoton={mostEnergeticPhoton, Z->daughter(0)};
  else if (fabs(physmath::deltaR(Z->daughter(0), mostEnergeticPhoton))>fabs(physmath::deltaR(Z->daughter(1), mostEnergeticPhoton)) )
    nearestChLeptToPhoton={mostEnergeticPhoton, Z->daughter(1)};
  
  double rawDFqJ0 = recoV.daughter(0).deepFlavour().probuds;
  double rawDFqJ1 = recoV.daughter(1).deepFlavour().probuds;
  double rawDFgJ0 = recoV.daughter(0).deepFlavour().probg;
  double rawDFgJ1 = recoV.daughter(1).deepFlavour().probg; 
  
  feat_dRLG = fabs(physmath::deltaR(nearestChLeptToPhoton.first, nearestChLeptToPhoton.second));
  feat_recoVMass = recoV.mass();
  feat_mllG=  llPh.M();
  //feat_J1PG = DFgJ1;
  feat_PhMVAId=selectedphotons.at(0).MVAvalue();
  feat_J1Girth=recoV.daughter(1).girth();
  feat_ptJ1 =  recoV.daughter(1).pt();
  //feat_J0Puds = DFqJ0;
  feat_ptjj=  recoV.pt();
  feat_HT   =  lljjPh.Pt();    
  feat_J0Girth=recoV.daughter(0).girth();
  feat_deltaR_J0Gamma=fabs(physmath::deltaR(recoV.daughter(0), selectedphotons.at(0)));
  feat_deltaR_J1Gamma=fabs(physmath::deltaR(recoV.daughter(1), selectedphotons.at(0)));
  feat_ptGamma=selectedphotons.at(0).pt();
  feat_dPhiZG=fabs(physmath::deltaPhi(Z->phi(),selectedphotons.at(0).phi()) );

  //refedined features
  feat_j0_btagDeepFlavQG = rawDFqJ0 / ( rawDFqJ0 + rawDFgJ0 );
  feat_j1_btagDeepFlavQG = rawDFqJ1 / ( rawDFqJ1 + rawDFgJ1 );
  
  if(jetSystApplied.direction>0){
    //    feat_ptJ0=ptj0Scaled_JUncUp;
    feat_ptJ1=ptj1Scaled_JUncUp;
    feat_ptjj=ptjjScaled_JUncUp;
    feat_recoVMass = mjjScaled_JUncUp;
  }else if(jetSystApplied.direction<0){
    //    feat_ptJ0=ptj0Scaled_JUncDn;
    feat_ptJ1=ptj1Scaled_JUncDn;
    feat_ptjj=ptjjScaled_JUncDn;
    feat_recoVMass = mjjScaled_JUncDn;
  }

  VZGMVAScore = reader->EvaluateMVA("BDT");

  return VZGMVAScore;/*  if(!isForDFsys)  return VZGMVAScore;  
  //old version of unc on DFscore

  double DFqJ0Up = TMath::Min(1.0, feat_J0Puds * 1. + DFQG_RelVar);
  double DFqJ0Dn = TMath::Max(0.0, feat_J0Puds * 1. - DFQG_RelVar);
  double DFgJ1Up = TMath::Min(1.0, feat_J1PG   * 1. + DFQG_RelVar);
  double DFgJ1Dn = TMath::Max(0.0, feat_J1PG   * 1. - DFQG_RelVar);

  feat_J0Puds    = DFqJ0Up;
  feat_J1PG      = DFgJ1Up;
  double BDT_DF_qJ0up_gJ1up   = reader->EvaluateMVA("BDT");
  feat_J1PG      = DFgJ1Dn;
  double BDT_DF_qJ0up_gJ1dn   = reader->EvaluateMVA("BDT");
  feat_J0Puds    = DFqJ0Dn;
  double BDT_DF_qJ0dn_gJ1dn   = reader->EvaluateMVA("BDT");
  feat_J1PG      = DFgJ1Up;
  double BDT_DF_qJ0dn_gJ1up   = reader->EvaluateMVA("BDT");

  double BDT_DFvaried=-2.;
  
  if(isForDFsys && isJERorJES==0){//DF sys case
    if(DFsysDirection>0){
      BDT_DFvaried = std::max({BDT_DF_qJ0up_gJ1up, BDT_DF_qJ0up_gJ1dn, BDT_DF_qJ0dn_gJ1up, BDT_DF_qJ0dn_gJ1dn});  //up variation
    }else if(DFsysDirection<0){
      BDT_DFvaried = std::min({BDT_DF_qJ0up_gJ1up, BDT_DF_qJ0up_gJ1dn, BDT_DF_qJ0dn_gJ1up, BDT_DF_qJ0dn_gJ1dn});  //dn variation
    }
  }
  if(BDT_DFvaried>-1.)
    return BDT_DFvaried;
  */
}

bool VZGAnalyzer::inSR(phys::Boson<phys::Jet> recoV, phys::Jet recoFJ, std::vector<phys::Photon> selectedphotons, int VBTopo, double VZGMVAScore)
{
  std::vector<std::string> orders = {"0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13"};

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



  phys::Photon mostEnergeticPhoton;

  std::stable_sort(selectedphotons.begin(), selectedphotons.end(), phys::EComparator());
  mostEnergeticPhoton = selectedphotons[0];


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
  std::pair<phys::Photon, phys::Lepton> nearestChLeptToPhoton;

  if(fabs(physmath::deltaR(Z->daughter(0), mostEnergeticPhoton))<fabs(physmath::deltaR(Z->daughter(1), mostEnergeticPhoton)) )
    nearestChLeptToPhoton={mostEnergeticPhoton, Z->daughter(0)};
  else if (fabs(physmath::deltaR(Z->daughter(0), mostEnergeticPhoton))>fabs(physmath::deltaR(Z->daughter(1), mostEnergeticPhoton)) )
    nearestChLeptToPhoton={mostEnergeticPhoton, Z->daughter(1)};


  
  VZGMVAScore= -2.;

  VZGMVAScore=VZGMVAScoreEval(recoV, recoFJ,  selectedphotons, VBTopo, noJERCsys, false, 0);
  
  return cut(7, recoV, recoFJ, selectedphotons, VBTopo, VZGMVAScore);

}

bool VZGAnalyzer::inCRZOFF( phys::Boson<phys::Jet> recoV, phys::Jet recoFJ, std::vector<phys::Photon> selectedphotons, int VBTopo, double VZGMVAScore)
{
  if( !cut(4, recoV, recoFJ, selectedphotons, VBTopo, VZGMVAScore)) return false; //out of common baseline
  return Z->mass()<CRZOFFMassCut;
  //  return Z->mass()<CRZOFFMassCut && cut(4, recoV, recoFJ, selectedphotons, VBTopo, VZGMVAScore);
}

bool VZGAnalyzer::inCRFSRTight( phys::Boson<phys::Jet> recoV, phys::Jet recoFJ, std::vector<phys::Photon> selectedphotons, int VBTopo, double VZGMVAScore)
{
   if( !cut(4, recoV, recoFJ, selectedphotons, VBTopo, VZGMVAScore)) return false; //out of common bas
   return !cut(7, recoV, recoFJ, selectedphotons, VBTopo, VZGMVAScore);
}

bool VZGAnalyzer::inCRZON_FSRTight( phys::Boson<phys::Jet> recoV, phys::Jet recoFJ, std::vector<phys::Photon> selectedphotons, int VBTopo, double VZGMVAScore)
{
  return Z->mass()>CRZONMassCut && !cut(7, recoV, recoFJ, selectedphotons, VBTopo, VZGMVAScore);
  //  return Z->mass()<CRZOFFMassCut && cut(4, recoV, recoFJ, selectedphotons, VBTopo, VZGMVAScore);
}

bool VZGAnalyzer::inCRZOFF_DIB( phys::Boson<phys::Jet> recoV, phys::Jet recoFJ, std::vector<phys::Photon> selectedphotons, int VBTopo, double VZGMVAScore)
{
  return inCRZOFF(recoV, recoFJ, selectedphotons, VBTopo, VZGMVAScore) && cut(9, recoV, recoFJ, selectedphotons, VBTopo, VZGMVAScore);
  //  return Z->mass()<CRZOFFMassCut && cut(4, recoV, recoFJ, selectedphotons, VBTopo, VZGMVAScore);
}

bool VZGAnalyzer::inCRZOFF_FSRTight( phys::Boson<phys::Jet> recoV, phys::Jet recoFJ, std::vector<phys::Photon> selectedphotons, int VBTopo, double VZGMVAScore)
{
  return inCRZOFF(recoV, recoFJ, selectedphotons, VBTopo, VZGMVAScore) && inCRFSRTight(recoV, recoFJ, selectedphotons, VBTopo, VZGMVAScore);
  //  return Z->mass()<CRZOFFMassCut && cut(4, recoV, recoFJ, selectedphotons, VBTopo, VZGMVAScore);
}


bool VZGAnalyzer::inCR2P_1VL(  phys::Boson<phys::Jet> recoV, phys::Jet recoFJ, std::vector<phys::Photon> selectedVLPhotons, int VBTopo, double VZGMVAScore)
{//Dec 25: implemented as CR2P_kinBut!WP90 
  if( !cut(ANALYSIS_CUTs_WP, recoV, recoFJ, selectedVLPhotons, VBTopo, VZGMVAScore)) return false; //out of common bas
  bool looseGammaExists=false;
  bool VLGammaExists=false;
  for (auto p : *photons){
    //    if(!looseGammaExists && p.id() == 22 && KinematicsOK(p, 20, 2.4) && !p.hasPixelSeed() && p.passElectronVeto() && p.cutBasedID(Photon::IdWp::VeryLoose)){
    if(!looseGammaExists && p.id() == 22 && KinematicsOK(p, 20, 2.4) && !p.hasPixelSeed() && p.passElectronVeto()){// && p.cutBasedID(Photon::IdWp::VeryLoose)){
      VLGammaExists=true;
      looseGammaExists = p.passMVA(Photon::MVAwp::wp90);//p.cutBasedIDLoose();
    }
  }
   
  return VLGammaExists && !looseGammaExists;//actually kin but not MVA WP 90 passed
}

bool VZGAnalyzer::inCR2P_1L(  phys::Boson<phys::Jet> recoV, phys::Jet recoFJ, std::vector<phys::Photon> selectedLoosePhotons, int VBTopo, double VZGMVAScore)
{
  if( !cut(4, recoV, recoFJ, selectedLoosePhotons, VBTopo, VZGMVAScore)) return false; //out of common bas
  bool looseGammaExists=false;
  bool mediumGammaExists=false;
  for (auto p : *photons){
    if(!mediumGammaExists && p.id() == 22 && KinematicsOK(p, 20, 2.4) && !p.hasPixelSeed() && p.passElectronVeto() && p.cutBasedIDLoose()){
      looseGammaExists=true;
      mediumGammaExists = p.cutBasedIDMedium();
    }
  }
   
  return looseGammaExists && !mediumGammaExists;
}

bool VZGAnalyzer::inCRVSide(  phys::Boson<phys::Jet> recoV, phys::Jet recoFJ, std::vector<phys::Photon> selectedphotons, int VBTopo, double VZGMVAScore)
{
   if( !cut(4, recoV, recoFJ, selectedphotons, VBTopo, VZGMVAScore)) return false; //out of common bas
   return true;//TO IMPLEMENT
}
bool VZGAnalyzer::inCRZSide( phys::Boson<phys::Jet> recoV, phys::Jet recoFJ, std::vector<phys::Photon> selectedphotons, int VBTopo, double VZGMVAScore)
{
   if( !cut(4, recoV, recoFJ, selectedphotons, VBTopo, VZGMVAScore)) return false; //out of common bas
   return Z->mass()<CRZinfMassCut && Z->mass()>CRZsupMassCut;

}


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


Bool_t VZGAnalyzer::cut(Int_t n, phys::Boson<phys::Jet> recoV, phys::Jet recoFJ, std::vector<phys::Photon> selectedPhotons, int VBTopo, double& VZGMVAScore)
{ // returns false if the event has to be cut
  //  std::cout<<"entering cut"<<std::endl;
  if(n<0)  return cut(8, recoV, recoFJ, selectedPhotons, VBTopo, VZGMVAScore);
  if(n==0) return true;
  if(n==1) return selectedPhotons.size()>=1;

  
  TLorentzVector llGamma = Z->daughter(0).p4()+Z->daughter(1).p4()+selectedPhotons.at(0).p4();
  double mllGamma=llGamma.M();

  phys::Photon mostEnergeticPhoton;
  //  std::stable_sort(selectedPhotons.begin(), selectedPhotons.end(), phys::EComparator());

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
  bool objectsExist = (VBTopo!=0 && selectedPhotons.size()>0);

  //addition for CRZOFF/ON
  bool isZOff= Z->mass() < 81;
  //bool ZMassWindow = isZOff;//Z->mass() > 60 && Z->mass() < 120;
  //bool ZMassWindow = Z->mass() > 60 && Z->mass() < 120 && !isZOff;
  bool ZMassWindow = Z->mass() > 60 && Z->mass() < 120;
  bool areGoodZCand = KinematicsOK(Z->daughter(0), 5.,etacut) && KinematicsOK(Z->daughter(1), 5.,etacut);

  bool JGsolved = fabs(physmath::deltaR(recoV.daughter(0), mostEnergeticPhoton))>0.4 && fabs(physmath::deltaR(recoV.daughter(1), mostEnergeticPhoton) )>0.4;
  bool LGsolved = fabs(physmath::deltaR(Z->daughter(0), mostEnergeticPhoton))>0.4 && fabs(physmath::deltaR(Z->daughter(1), mostEnergeticPhoton) )>0.4;
  //if(VBTopo==-1)  JGsolved = physmath::deltaR(recoFJ, mostEnergeticPhoton)>0.8;
											     
  bool baseline = (objectsExist && VBTopo==1  && ZMassWindow   && areGoodZCand && LGsolved && JGsolved);
  bool MVAScoreWP0 = VZGMVAScore > 0.;
  bool MVAScoreWP7 = VZGMVAScore > 0.7;
  bool MVAScoreWP8 = VZGMVAScore > 0.8;
  bool MVAScoreWP9 = VZGMVAScore > 0.9;
  
  
  bool noFSR_ifFiducial = (!fiducial_run || mllGamma>-2.5*Z->mass()+315);

  switch (n)
  {
  case 2://baseline (with dRjG) 
    if (objectsExist)// && Z->mass() > 80)
      return true;
    break;
  case 3:
    if (objectsExist   && ZMassWindow   && areGoodZCand)
      return true;
    break;
  case 4:
    if (baseline)//	&& MVAScoreWP0)
      return true;
    break;
  case 5://high mllG
    if (baseline
	&& mllGamma>90)
      //	&& mllGamma<95)//CUT 6 REVERTED FOR CRFSRTight 
      return true;
    break;
  case 6://high mllG
    if (baseline
	&& mllGamma>95)
      //	&& mllGamma<95)//CUT 6 REVERTED FOR CRFSRTight 
      return true;
    break;
  case 7://high mllG
    if (baseline
	&& mllGamma>100)
      return true;
    break;
  case 8://high mllG
    if (baseline
	&& mllGamma>120)
      return true;
    break;
  case 9://high mllG
    if (baseline
	&& mllGamma>140)//CUT 9 REVERTED FOR CRFSRMedium 
      return true;
    break;
  case 10://dRlG>0.5
    if (baseline
	&& fabs(physmath::deltaR(nearestLepToPhoton.first, nearestLepToPhoton.second))> 0.5
	&& mllGamma>140)//	&& mllGamma>150)
      return true;
    break;
  case 11:
    if(baseline
       &&       mostEnergeticPhoton.passMVA(Photon::MVAwp::wp80)
       &&       mllGamma>140)
      return true;
    break;
   
    /*    
  case 11://MET<90
    if (baseline
	&& met->pt()<90
	&& mllGamma>140)
      return true;
    break;
  case 12://MET<60
    if (baseline
	&& met->pt()<60
	&& mllGamma>140)
      return true;
    break;
  case 13://MET<50
    if (baseline
	&& met->pt()<50
	&& mllGamma>140)
      return true;
    break;
  case 14://MET<50
    if (baseline
	&& ! std::any_of(jets->begin(), jets->end(), [](const Jet& j){ auto dF = j.deepFlavour(); return dF.probb + dF.probbb + dF.problepb > 0.2770; })    // Note: this is the medium WP for Legacy samples (102X)
	&& met->pt()<50
	&& mllGamma>140)
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
  //if(IsARunForMVAFeat) return;
  //bool isBkg=false;
  
  int nbOfGenQuarks=0;

  foreach (const Particle &p, *genParticles)
    {
      if  (abs(p.id()) < 10) // Is it a quark? 
	{
	  nbOfGenQuarks++;
		     
	}
    }

  //  theHistograms->fill("#######nbOfGenQuarks"  , "########nbOfGenQuarks"               , 5, -0.5, 4.5, nbOfGenQuarks, theWeight);

  
  if(verboseControlBlinding){
    std::cout<<"TF="<<TF<<endl;
    std::cout<<"Starting Eff SF = "<<PhEffSF<<endl;
  }
  PhEffSF=1.;
  PhEffSFUnc=0.;

  if(verboseControlBlinding)
    std::cout<<"For this event Eff SF = "<<PhEffSF<<endl;
  
  bool isSigSample = theSampleInfo.isMC() && (theSampleInfo.fileName().find("WZG")!=std::string::npos || theSampleInfo.fileName().find("ZZG")!=std::string::npos || theSampleInfo.fileName().find("ZH")!=std::string::npos);
  bool isDYSample = theSampleInfo.isMC() && theSampleInfo.fileName().find("DY")!=std::string::npos;
  bool isZGSample = theSampleInfo.isMC() && !isSigSample && (theSampleInfo.fileName().find("ZG")!=std::string::npos);

  double VZGMVAScore = -2.;
  std::string region = "";

  if(verbose==true)    {

    cout << "----------------------------------------------------------------" << endl;
    cout << "Run: " << run << " event: " << event << endl;
    
    if(theSampleInfo.isMC()){
      cout << "----------------------------------------------------------------" << endl;
      cout << "MC sample" << endl;
    }else{
      cout << "----------------------------------------------------------------" << endl;
      cout << "DATA sample" << endl;
    }

    if(theSampleInfo.isMC() && (theSampleInfo.fileName().find("WZG")!=std::string::npos || theSampleInfo.fileName().find("ZZG")!=std::string::npos || theSampleInfo.fileName().find("ZH")!=std::string::npos ) ){
      cout << "----------------------------------------------------------------" << endl;
      cout << "signal sample" << endl;
    }else{
      cout << "----------------------------------------------------------------" << endl;
      cout << "non-sig sample" << endl;
    }
  }
  //  genAnalyze();
  
  //----BLOCK ASSIGNING DY REWGT PER YEAR------------//
  std::string year_str = std::to_string(year);
  rewgt=1.;
  rewErr=0.;
  if(theSampleInfo.isMC() && isDYSample && genVBHelper_.ZtoChLep().size()>0){    
    for(int i = 0; i<DYrewgtBinEdges.size()-1 && rewgt==1.; i++){
      if( genVBHelper_.ZtoChLep()[0].pt() > DYrewgtBinEdges[i] && genVBHelper_.ZtoChLep()[0].pt() < DYrewgtBinEdges[i+1]){
	/*
	if(year==2016){
	  if(year_str.find("preVFP")!=std::string::npos)      rewgt=DYreweights_2016preVFP[i];
	  else rewgt=DYreweights_2016postVFP[i];
	}
	if(year==2017) rewgt=DYreweights_2017[i];
	if(year==2018) rewgt=DYreweights_2018[i];
	*/
	rewgt=DYreweights_Run2[i];
	rewErr=DYrewErr_Run2[i];
      }
    }
  }
  //rewgt=1.; //CT: momentaneously de-activating DY reweighting and leaving ZG weights as they are. TO BE REMOVED AFTERWARDS 
    /*
    if (genVBHelper_.ZtoChLep()[0].pt() > DYrewgtBinEdges[DYrewgtBinEdges.size()-1]){ //exception: overflow
      if(year==2016){
	if(year_str.find("preVFP")!=std::string::npos)      rewgt=DYreweights_2016preVFP[DYreweights_2016preVFP.size()-1];
	else rewgt=DYreweights_2016postVFP[reweights_2016postVFP.size()-1];
      }
      if(year==2017) rewgt=DYreweights_2017[DYreweights_2017.size()-1];
      if(year==2018) rewgt=DYreweights_2018[DYreweights_2018.size()-1];

    }
    }*/
  //----BLOCK ASSIGNING ZG REWGT PER YEAR------------//
  if(theSampleInfo.isMC() && isZGSample && genVBHelper_.ZtoChLep().size()>0){
    for(int i = 0; i<ZGrewgtBinEdges.size()-1 && rewgt==1.; i++){
      if( genVBHelper_.ZtoChLep()[0].pt() > ZGrewgtBinEdges[i] && genVBHelper_.ZtoChLep()[0].pt() < ZGrewgtBinEdges[i+1]){
	/*
	if(year==2016){
	  if(year_str.find("preVFP")!=std::string::npos)      rewgt=ZGreweights_2016preVFP[i];
	  else rewgt=ZGreweights_2016postVFP[i];
	}
	if(year==2017) rewgt=ZGreweights_2017[i];
	if(year==2018) rewgt=ZGreweights_2018[i];
	*/
	rewgt=ZGreweights_Run2[i];
	rewErr=ZGrewErr_Run2[i];
      }
    }
  }
  //  rewgt=1.; //CT: momentaneously de-activating ZG reweighting and leaving ZG weights as they are. TO BE REMOVED AFTERWARDS 

    //-------------------------------------------------//
  
  
  int VBTopo = 0;
  phys::Boson<phys::Jet> recoV;
  phys::Jet recoFJ;
  bool haveGoodRECODiJetCand=false;
  bool haveGoodRECOFJCand=false;

  int VBTopo_2P1VL = 0;
  phys::Boson<phys::Jet> recoV_2P1VL;
  phys::Jet recoFJ_2P1VL;
  bool haveGoodRECODiJetCand_2P1VL=false;
  bool haveGoodRECOFJCand_2P1VL=false;

  //_____________________________________________________________________
  //CT: TEMP block for ZFSR subtractioN from DY
  bool promptPhExists = false;
  bool isFSR=false;
  std::vector<phys::Particle> selectedGENphotons, selectedPROMPTphotons;
  for (auto p : *genParticles)
    if (p.id() == 22 && KinematicsOK(p, 15, 2.4)){// && p.genStatusFlags().test(phys::isPrompt) &&  p.genStatusFlags().test(phys::fromHardProcess))
      selectedGENphotons.push_back(p);
      if (p.genStatusFlags().test(phys::isPrompt)) selectedPROMPTphotons.push_back(p);
    }
  if (selectedGENphotons.size()>=1){
    std::stable_sort(selectedGENphotons.begin(), selectedGENphotons.end(), phys::EComparator());
    promptPhExists=selectedPROMPTphotons.size()>=1;
    TLorentzVector  GEN_llPh ;
    double mGEN_llPh;
    double mGenZ;
    if(genVBHelper_.ZtoChLep().size()>=1){
      GEN_llPh = genVBHelper_.ZtoChLep()[0].p4()+selectedGENphotons.at(0).p4();
      mGEN_llPh=GEN_llPh.M();
      mGenZ=genVBHelper_.ZtoChLep()[0].mass();

      isFSR=mGEN_llPh<315-2.5*mGenZ;
    }
  }
  //______________END OF TEMP BLOCK_______________________________________

  std::vector<phys::Photon> selectedVLPhotons;
  PhotonVLSelection(&selectedVLPhotons, 1);
  
  std::vector<phys::Photon> selectedphotons;
  PhotonSelection(&selectedphotons);
  //  std::cout<<"selected photons size "<<selectedphotons.size()<<std::endl;

  bool isCR = false;
  if(selectedphotons.size()<1){
    if (IN_GENsignalDef()){
      printHistos(0, "sign", recoV, recoFJ, selectedphotons, VBTopo, region, isCR);   
    }else if(!IN_GENsignalDef() ){
      printHistos(0, "bckg", recoV, recoFJ, selectedphotons, VBTopo, region, isCR);
    }
    printHistos(0, "all", recoV, recoFJ, selectedphotons, VBTopo, region, isCR);
    /*
    if(!isFSR){
      printHistos(0, "dib", recoV, recoFJ, selectedphotons, VBTopo, region, isCR);
    }else{
      printHistos(0, "fsr", recoV, recoFJ, selectedphotons, VBTopo, region, isCR);
    }
    */
    if(!promptPhExists) {
      printHistos(0, "nonPrompt", recoV, recoFJ, selectedphotons, VBTopo, region, isCR);
    }else{
      printHistos(0, "prompt", recoV, recoFJ, selectedphotons, VBTopo, region, isCR);
    }
    //    std::cout<<"No Loose Photons "<<selectedphotons.size()<<std::endl;
    isCR = true;
    if(selectedVLPhotons.size()<1){

      //            std::cout<<"No Loose Photons, but no VL Photons either"<<selectedVLPhotons.size()<<std::endl;
      return;//NOTE: it does exist at least 1 VL photon but there are no Loose Photons at all 
    }


    //std::cout<<selectedVLPhotons.size()<<" VL Photons "<<std::endl;
    //std::cout<<selectedphotons.size()<<" Loose Photons "<<std::endl;
    phys::Photon mostEnergeticPhoton_2P1VL;

    mostEnergeticPhoton_2P1VL = selectedVLPhotons[0];
  
    VBTopo_2P1VL=Reconstruct(&recoV_2P1VL,&recoFJ_2P1VL,&haveGoodRECODiJetCand_2P1VL,&haveGoodRECOFJCand_2P1VL,&mostEnergeticPhoton_2P1VL, false, noJERCsys);

    if(inCR2P_1VL( recoV_2P1VL, recoFJ_2P1VL, selectedVLPhotons, VBTopo_2P1VL, VZGMVAScore) ){
      region = "CR2P_1VL";
	  //    if(inCR2P_1L( recoV_2P1VL, recoFJ_2P1VL, selectedVLPhotons, VBTopo_2P1VL, VZGMVAScore) ){
	  //      region = "CR2P_1L";

      //      std::cout<<"region: "<<region<<endl;

      if (IN_GENsignalDef()){
	//	printHistos(ANALYSIS_CUTs_WP, "sign", recoV_2P1VL, recoFJ_2P1VL, selectedVLPhotons, VBTopo_2P1VL, region, isCR);
	printHistos(1, "sign", recoV_2P1VL, recoFJ_2P1VL, selectedVLPhotons, VBTopo_2P1VL, region, isCR);   
      }else if (!IN_GENsignalDef()){
	//	printHistos(ANALYSIS_CUTs_WP, "bckg", recoV_2P1VL, recoFJ_2P1VL, selectedVLPhotons, VBTopo_2P1VL, region, isCR);
        printHistos(1, "bckg", recoV_2P1VL, recoFJ_2P1VL, selectedVLPhotons, VBTopo_2P1VL, region, isCR);
      }
      //      printHistos(ANALYSIS_CUTs_WP, "all", recoV_2P1VL, recoFJ_2P1VL, selectedVLPhotons, VBTopo_2P1VL, region, isCR);
      printHistos(1, "all", recoV_2P1VL, recoFJ_2P1VL, selectedVLPhotons, VBTopo_2P1VL, region, isCR); 
      /*
      if(!isFSR){
	printHistos(ANALYSIS_CUTs_WP, "dib", recoV, recoFJ, selectedVLPhotons, VBTopo_2P1VL, region, isCR);
      }else{
	printHistos(ANALYSIS_CUTs_WP, "fsr", recoV, recoFJ, selectedVLPhotons, VBTopo_2P1VL, region, isCR);      
      }
      */
      if(!promptPhExists){
	//	printHistos(ANALYSIS_CUTs_WP, "nonPrompt", recoV_2P1VL, recoFJ_2P1VL, selectedVLPhotons, VBTopo_2P1VL, region, isCR);
	printHistos(1, "nonPrompt", recoV_2P1VL, recoFJ_2P1VL, selectedVLPhotons, VBTopo_2P1VL, region, isCR);
      }else{
	//	printHistos(1, "ANALYSIS_CUTs_WP", recoV_2P1VL, recoFJ_2P1VL, selectedVLPhotons, VBTopo_2P1VL, region, isCR);
	printHistos(1, "prompt", recoV_2P1VL, recoFJ_2P1VL, selectedVLPhotons, VBTopo_2P1VL, region, isCR);
      }
    }
    return;
  }
  
  //  std::cout<<"PASSING selected photons size "<<selectedphotons.size()<<std::endl;

  //________________________________________________________________________________________________________
  // NO MORE EVENTS W/O PHOTONS PROCESSED FROM HERE ON

  std::vector<phys::Boson<phys::Particle>> genV;
  //std::vector<phys::Boson<phys::Particle>> genV(genVBHelper_.ZtoQ().size()+genVBHelper_.WtoQ().size());
  if(genVBHelper_.ZtoQ().size()>0){
    genV=genVBHelper_.ZtoQ();
    genV.insert(genV.end(), genVBHelper_.WtoQ().begin(), genVBHelper_.WtoQ().end());
  }else if(genVBHelper_.WtoQ().size()>0){
    genV=genVBHelper_.WtoQ();
    genV.insert(genV.end(), genVBHelper_.ZtoQ().begin(), genVBHelper_.ZtoQ().end());
  }
    //---------------------------------------- Single q analysis ----------------------------------------//
  std::vector<phys::Particle> genQuarksfromV;
  for (auto VB : genV){
    /*
    theHistograms->fill("quarkfromV charge", "quarkfromV charge", 7, -7. / 6., 7. / 6., VB.daughter(0).charge(), theWeight*LumiSF);
    theHistograms->fill("quarkfromV charge", "quarkfromV charge", 7, -7. / 6., 7. / 6., VB.daughter(1).charge(), theWeight*LumiSF);
    theHistograms->fill("quarkfromV pt", "quarkfromV pt", 50, 0, 600, VB.daughter(0).pt(), theWeight*LumiSF);
    theHistograms->fill("quarkfromV pt", "quarkfromV pt", 50, 0, 600, VB.daughter(1).pt(), theWeight*LumiSF);
    */
    genQuarksfromV.push_back(VB.daughter(0));
    genQuarksfromV.push_back(VB.daughter(1));
  }
  //theHistograms->fill("0size_GENQuarksfromV_beforecuts", "0size_GENQuarksfromV_beforecuts", 10, -0.5, 9.5, genQuarksfromV.size(), theWeight*LumiSF);
  //genQuarksfromV.erase(std::remove_if(genQuarksfromV.begin(), genQuarksfromV.end(),
  //				      [](phys::Particle p){ return !KinematicsOK(p, ptcut, etacut); }),
  //		       genQuarksfromV.end());
  //theHistograms->fill("0size_GENQuarksfromV_aftercuts", "0size_GENQuarksfromV_aftercuts", 10, -0.5, 9.5, genQuarksfromV.size(), theWeight*LumiSF);
 
  //----------------------------------------Kinematic Cuts on VB(qq)----------------------------------------//
  std::vector<phys::Boson<phys::Particle>> DiQuarks=genV;

  //DiQuarks.erase(std::remove_if(DiQuarks.begin(), DiQuarks.end(), [](phys::Boson<phys::Particle> VB)
  //                                      { return !(KinematicsOK(VB.daughter(0), ptcut, etacut)&&KinematicsOK(VB.daughter(1), ptcut, etacut)); }),
  //                       DiQuarks.end());

 //----------------------------------------Kinematic Cuts GEN Jets AK4 & RECO Jets AK4----------------------------------------//
  std::vector<phys::Particle> selectedGENjets;
  foreach (const phys::Particle &jet, *genJets){
    if (KinematicsOK(jet,ptcut,etacut)){
      selectedGENjets.push_back(jet);
    }
  }
  std::vector<phys::Jet> selectedRECOjets;
  foreach (const phys::Jet &jet, *jets){
    if (KinematicsOK(jet,ptcut,etacut)){
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
      if (fabs(physmath::deltaR(quark, nearestjet)) < 0.4){
        jetsfromquarks.push_back(nearestjet);
        makesjet = true;
      }
    }
    if (makesjet && selectedGENjets.size() > 0){
      theHistograms->fill("#QUARK=>GEN", "#QUARK=>GEN", 2, 0, 2, 1., theWeight*LumiSF);
      theHistograms->fill("Pt_quark_num", " Pt_quark_num; GeV/c", 10, ptcut, 300, quark.pt(), theWeight*LumiSF);
    }
    else{
      theHistograms->fill("#QUARK=>GEN", "#QUARK=>GEN", 2, 0, 2, 0., theWeight*LumiSF);
    }
  }
  for (auto pair : nearestjetstoquark){
    //ResolutionPlots(pair.first,pair.second,"Hadronization_",theWeight*LumiSF,"");
    theHistograms->fill("DeltaR_quark_vs_BestMatchedGENJet", "DeltaR_quark_vs_BestMatchedGENJet; #DeltaR", 20, 0, 0.5, fabs(physmath::deltaR(pair.first, pair.second)), theWeight*LumiSF);
    theHistograms->fill("DeltaR_quark_jet_vs_pt", "DeltaR vs pt;pt [GeV/c] ; #DeltaR", 10, ptcut, 300, 20, 0, 0.2, pair.first.pt(),fabs(physmath::deltaR(pair.first, pair.second)), theWeight*LumiSF);
  }
  //theHistograms->fill("1size_GENjetsfromquarks", "size_GENjetsfromquarks", 10, -0.5, 9.5, jetsfromquarks.size(), theWeight*LumiSF);

  //----------------------------------------Matching efficiency ______ SINGLE GENJET/SINGLE RECOJET--------------//
  std::vector<phys::Particle> RECOjetsfromGENjets;
  std::vector<std::pair<phys::Particle, phys::Particle>> nearestRECOjetstoGENjets;

  for (auto genJet : selectedGENjets)
  {
    phys::Particle nearestRECOjet;
    bool isreconstructed = false;
    theHistograms->fill("Pt_genJet_den", " Pt_genJet_den; GeV/c", 10, ptcut, 300, genJet.pt(), theWeight*LumiSF);

    if (selectedRECOjets.size() > 0){
      std::stable_sort(selectedRECOjets.begin(), selectedRECOjets.end(), phys::DeltaRComparator(genJet));
      nearestRECOjet = selectedRECOjets.at(0);
      nearestRECOjetstoGENjets.push_back({genJet, nearestRECOjet});
      if (fabs(physmath::deltaR(genJet, nearestRECOjet)) < 0.4){
        RECOjetsfromGENjets.push_back(nearestRECOjet);
        isreconstructed = true;
      }
    }
    if (isreconstructed && selectedRECOjets.size() > 0){
      theHistograms->fill("#GEN=>RECO", "#GEN=>RECO", 2, 0, 2, 1., theWeight*LumiSF);
      theHistograms->fill("Pt_genJet_num", " Pt_genJet_num; GeV/c", 10, ptcut, 300, genJet.pt(), theWeight*LumiSF);
    }
    else{
      theHistograms->fill("#GEN=>RECO", "#GEN=>RECO", 2, 0, 2, 0., theWeight*LumiSF);
    }
  }
  for (auto pair : nearestRECOjetstoGENjets){
    //    ResolutionPlots(pair.first,pair.second,"SingleJetsReconstruction_",theWeight*LumiSF,"");
    theHistograms->fill("DeltaR_GENjet_vs_BestMatchedRECOJet", "DeltaR_GENjet_vs_BestMatchedRECOJet; #DeltaR", 20, 0, 0.5, fabs(physmath::deltaR(pair.first, pair.second)), theWeight*LumiSF);
    theHistograms->fill("DeltaR_jets_vs_pt", "DeltaR jets vs pt;pt [GeV/c] ; #DeltaR", 10, ptcut, 300, 20, 0, 0.2, pair.first.pt(),fabs(physmath::deltaR(pair.first, pair.second)), theWeight*LumiSF);
  }
  //theHistograms->fill("2size_GENjetsRECONSTRUCTED", "2size_GENjetsRECONSTRUCTED", 10, -0.5, 9.5, RECOjetsfromGENjets.size(), theWeight*LumiSF);
  
  //----------------------------------------Matching efficiency ______ QUARKS PAIR------------------------//

  std::vector<std::pair<phys::Particle, phys::Particle>> DijetsmatchedtoDiquark;
  std::vector<phys::Boson<phys::Particle>> DiJetsGEN;
  //std::cout << ".................GEN to QUARKS MATCHING..............." << std::endl;

  //std::cout << "DiQuarks size: " << DiQuarks.size() << std::endl;
  for (auto Diquark : DiQuarks){

    bool firstmatches = false;
    phys::Particle jetmatchedtoFIRSTquark;

    bool secondmatches = false;
    phys::Particle jetmatchedtoSECONDquark;

    bool atleastonematches = false;
    bool bothmatch = false;
    //    std::cout << "selectedGENjets size: " << selectedGENjets.size() << std::endl;
    for (auto genJet : selectedGENjets){
      double deltaR1 = fabs(physmath::deltaR(Diquark.daughter(0), genJet));
      //      std::cout << "deltaR1= " << deltaR1 << std::endl;
      double deltaR2 = fabs(physmath::deltaR(Diquark.daughter(1), genJet));
      //      std::cout << "deltaR2= " << deltaR2 << std::endl;

      if (deltaR1 < 0.4){
	//        std::cout << "first matched" << std::endl;
        jetmatchedtoFIRSTquark = genJet;
        firstmatches = true;
      }
      if (deltaR2 < 0.4){
        jetmatchedtoSECONDquark = genJet;
	//        std::cout << "second matched" << std::endl;
        secondmatches = true;
      }
    }

    bothmatch = (firstmatches && secondmatches);
    atleastonematches = (firstmatches || secondmatches);

    //    std::cout << "bothmatch= " << bothmatch << std::endl;
    //    std::cout << "atleastonematches= " << atleastonematches << std::endl;

    if (bothmatch){
      //      std::cout << "both matched" << std::endl;
      //      std::cout << "reconstructing a boson from dijets matched to diquark" << std::endl;
      DiJetsGEN.push_back(phys::Boson<phys::Particle>(jetmatchedtoFIRSTquark, jetmatchedtoSECONDquark));
      DijetsmatchedtoDiquark.push_back({Diquark, phys::Boson<phys::Particle>(jetmatchedtoFIRSTquark, jetmatchedtoSECONDquark)});
      //      theHistograms->fill("#Bothmatched", "#Bothmatched", 2, 0, 2, 1., theWeight*LumiSF);
    }
    theHistograms->fill("#Bothmatched", "#Bothmatched", 2, 0, 2, bothmatch, theWeight*LumiSF);
    theHistograms->fill("#AtLeastONEmatches", "#AtLeastONEmatches", 2, 0, 2, atleastonematches, theWeight*LumiSF);

  }
  //theHistograms->fill("1.1size_GEN_DiJets", "1size_GEN_DiJets", 10, -0.5, 9.5, DiJetsGEN.size(), theWeight*LumiSF);
  for (auto DiJet : DiJetsGEN)
  {
    float mjj = (DiJet.daughter(0).p4() + DiJet.daughter(1).p4()).M();
    theHistograms->fill("mjj_GEN", "mjj_GEN", 10, 50, 120, mjj, theWeight*LumiSF);
  }

  //----------------------------------------Matching efficiency ______ RECO to GEN  PAIR------------------------//
  //  std::cout << ".................RECO to GEN MATCHING..............." << std::endl;

  // std::vector<std::pair<phys::Particle,phys::Particle>> DiRECOjetsmatchedtoDiGENjets;
  // std::vector<phys::Boson<phys::Particle>> DiJetsRECO;
  std::vector<phys::Boson<phys::Particle>> DiJetsGENreconstructed;
  //std::cout << "GEN Dijets size: " << DiJetsGEN.size() << std::endl;
  for (auto DiJet : DiJetsGEN)
  {

    bool firstmatches = false;
    // phys::Particle jetmatchedtoFIRSTgen;

    bool secondmatches = false;
    // phys::Particle jetmatchedtoSECONDgen;

    bool atleastonematches = false;
    bool bothmatch = false;
    //    std::cout << "selectedRECOjets size: " << selectedRECOjets.size() << std::endl;
    for (auto recoJet : selectedRECOjets)
    {
      double deltaR1 = fabs(physmath::deltaR(DiJet.daughter(0), recoJet));
      //      std::cout << "deltaR1= " << deltaR1 << std::endl;
      double deltaR2 = fabs(physmath::deltaR(DiJet.daughter(1), recoJet));
      //      std::cout << "deltaR2= " << deltaR2 << std::endl;

      if (deltaR1 < 0.4)
      {
	//        std::cout << "first matched" << std::endl;
        // jetmatchedtoFIRSTgen=recoJet;
        firstmatches = true;
      }
      if (deltaR2 < 0.4)
      {
	//        std::cout << "second matched" << std::endl;
        // jetmatchedtoSECONDgen=recoJet;
        secondmatches = true;
      }
    }

    bothmatch = (firstmatches && secondmatches);
    atleastonematches = (firstmatches || secondmatches);

    //    std::cout << "bothmatch= " << bothmatch << std::endl;
    //    std::cout << "atleastonematches= " << atleastonematches << std::endl;

    if (bothmatch)
    {
      //      std::cout << "both matched" << std::endl;
      // std::cout << "reconstructing a boson from diRECOjets matched to diGENjets" << std::endl;
      // DiJetsRECO.push_back(phys::Boson<phys::Particle>(jetmatchedtoFIRSTgen, jetmatchedtoSECONDgen));
      DiJetsGENreconstructed.push_back(DiJet);
      // DiRECOjetsmatchedtoDiGENjets.push_back({DiJet,phys::Boson<phys::Particle>(jetmatchedtoFIRSTgen, jetmatchedtoSECONDgen)});

    }
    theHistograms->fill("#RECOGEN_Bothmatched", "#RECOGEN_Bothmatched", 2, 0, 2, bothmatch, theWeight*LumiSF);
    theHistograms->fill("#RECOGEN_AtLeastONEmatches", "#RECOGEN_AtLeastONEmatches", 2, 0, 2, atleastonematches, theWeight*LumiSF);
  }
  // theHistograms->fill("2size_RECO_DiJets", "1size_RECO_DiJets", 10, -0.5, 9.5, DiJetsRECO.size(), theWeight*LumiSF);
  //  for (auto DiJet:DiJetsRECO)
  //  {
  //      float mjj = (DiJet.daughter(0).p4() + DiJet.daughter(1).p4()).M();
  //      theHistograms->fill("mjj_RECO", "mjj_RECO", 10, 50, 120, mjj, theWeight*LumiSF);
  //  }


  
  //________________________________________________________________________________________________________
  
  // NEW BLOCK FOR VHad reconstruction test
  

  //----------------------------------------Matching efficiency ______ ALGORITHM------------------------//
  if(verbose==true)std::cout << ".................Algorithm efficiency..............." << std::endl;


  std::map<std::string, Boson<phys::Jet>> Candidates;
  std::vector<phys::Boson<phys::Jet>> JetPairs;
  /*
  std::vector<phys::Jet> selectedRECOjets;
  foreach (const phys::Jet &jet, *jets)
    if (KinematicsOK(jet,ptcut,etacut) && jet.passLooseJetID()) // KinematicsOK(jet)	
	selectedRECOjets.push_back(jet);
  *///filled before
  for (int i = 0; i < selectedRECOjets.size(); i++)    
    for (int j = i + 1; j < selectedRECOjets.size(); j++)	
	JetPairs.push_back(phys::Boson<phys::Jet>(selectedRECOjets.at(i), selectedRECOjets.at(j)));

  theHistograms->fill("size_JetsPairs_RECO", ">=2size_JetsPairs_RECO", 10, -0.5, 9.5, JetPairs.size(), theWeight*LumiSF);
  if(verbose==true)  std::cout << "#reco jets pairs : " << JetPairs.size() << std::endl;

  if (JetPairs.size() > 0){

    Candidates["mjjBased"] = *std::max_element(JetPairs.begin(), JetPairs.end(),
					       [this](const phys::Boson<phys::Jet>& DJA, const phys::Boson<phys::Jet>& DJB) {
						 return this->VHadScore(DJA, 1, 0, 0) < this->VHadScore(DJB, 1, 0, 0);
					       });

    Candidates["qglBased"] = *std::max_element(JetPairs.begin(), JetPairs.end(),
					       [this](const phys::Boson<phys::Jet>& DJA, const phys::Boson<phys::Jet>& DJB) {
						 return this->VHadScore(DJA, -1, 0, 0) < this->VHadScore(DJB, -1, 0, 0);
					       });

    Candidates["mixSmth01"] = *std::max_element(JetPairs.begin(), JetPairs.end(),
						[this](const phys::Boson<phys::Jet>& DJA, const phys::Boson<phys::Jet>& DJB) {
						  return this->VHadScore(DJA, 0, 0.5, 0.1) < this->VHadScore(DJB, 0, 0.5, 0.1); 	//VHadScore(DiJetsCand[j], 0, inflecPt, smoothness);
						});

    Candidates["mixSmth05"] = *std::max_element(JetPairs.begin(), JetPairs.end(),
						[this](const phys::Boson<phys::Jet>& DJA, const phys::Boson<phys::Jet>& DJB) {
						  return this->VHadScore(DJA, 0, 0.5, 0.5) < this->VHadScore(DJB, 0, 0.5, 0.5); 	//VHadScore(DiJetsCand[j], 0, inflecPt, smoothness);
						});
    
    Candidates["mixSmth10"] = *std::max_element(JetPairs.begin(), JetPairs.end(),
						[this](const phys::Boson<phys::Jet>& DJA, const phys::Boson<phys::Jet>& DJB) {
						  return this->VHadScore(DJA, 0, 0.5, 1.0) < this->VHadScore(DJB, 0, 0.5, 1.0); 	//VHadScore(DiJetsCand[j], 0, inflecPt, smoothness);
						});

    Candidates["mixSmth20"] = *std::max_element(JetPairs.begin(), JetPairs.end(),
						[this](const phys::Boson<phys::Jet>& DJA, const phys::Boson<phys::Jet>& DJB) {
						  return this->VHadScore(DJA, 0, 0.5, 2.0) < this->VHadScore(DJB, 0, 0.5, 2.0); 	//VHadScore(DiJetsCand[j], 0, inflecPt, smoothness);
						});
    Candidates["mixSmth40"] = *std::max_element(JetPairs.begin(), JetPairs.end(),
						[this](const phys::Boson<phys::Jet>& DJA, const phys::Boson<phys::Jet>& DJB) {
						  return this->VHadScore(DJA, 0, 0.5, 4.0) < this->VHadScore(DJB, 0, 0.5, 4.0); 	//VHadScore(DiJetsCand[j], 0, inflecPt, smoothness);
						});
    Candidates["mixSmth80"] = *std::max_element(JetPairs.begin(), JetPairs.end(),
						[this](const phys::Boson<phys::Jet>& DJA, const phys::Boson<phys::Jet>& DJB) {
						  return this->VHadScore(DJA, 0, 0.5, 8.0) < this->VHadScore(DJB, 0, 0.5, 8.0); 	//VHadScore(DiJetsCand[j], 0, inflecPt, smoothness);
						});
    Candidates["mixSmth100"] = *std::max_element(JetPairs.begin(), JetPairs.end(),
						[this](const phys::Boson<phys::Jet>& DJA, const phys::Boson<phys::Jet>& DJB) {
						  return this->VHadScore(DJA, 0, 0.5, 8.0) < this->VHadScore(DJB, 0, 0.5, 8.0); 	//VHadScore(DiJetsCand[j], 0, inflecPt, smoothness);
						});
    
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
      
    for (auto Candidate : Candidates)      {
      theHistograms->fill("mVHCand_" + Candidate.first + "_Candidate", "mjj_" + Candidate.first + "_Candidate", 10, 50, 120, Candidate.second.mass(), theWeight*LumiSF);
      theHistograms->fill("qglVHCand_" + Candidate.first + "_Candidate", "qglV_" + Candidate.first + "_Candidate", 10, 0, 1, TMath::Sqrt(Candidate.second.daughter(0).qgLikelihood()*Candidate.second.daughter(1).qgLikelihood()), theWeight*LumiSF);
    }
  }
  

  //  std::cout << "#true gen jets pairs reconstructed: " << DiJetsGENreconstructed.size() << std::endl;

  for (auto DiJet : DiJetsGENreconstructed)
  {
    theHistograms->fill("mjj_den" , " mjj_den; GeV/c^{2}" , 10, 50, 120, DiJet.mass(), theWeight*LumiSF);
    theHistograms->fill("Pt_den" , " Pt_den; GeV/c" , 10, 0, 300, DiJet.pt(), theWeight*LumiSF);


    //----------------------------------------MATCHED Jets total mass----------------------------------//
    for (auto Candidate : Candidates)
    {

      //        std::cout<<""<<std::endl;

        bool truepair = false;
	//        std::cout << "Algorithm: " << Candidate.first << std::endl;
        phys::Jet jetRECOA = Candidate.second.daughter(0);
        phys::Jet jetRECOB = Candidate.second.daughter(1);
        phys::Particle jetGENA = DiJet.daughter(0);
        phys::Particle jetGENB = DiJet.daughter(1);
        double deltaRAA = fabs(physmath::deltaR(jetGENA, jetRECOA));
        double deltaRAB = fabs(physmath::deltaR(jetGENA, jetRECOB));
        double deltaRBA = fabs(physmath::deltaR(jetGENB, jetRECOA));
        double deltaRBB = fabs(physmath::deltaR(jetGENB, jetRECOB));
        if ((deltaRAA < 0.4 && deltaRBB < 0.4) || (deltaRAB < 0.4 && deltaRBA < 0.4))
          truepair = true;
        if (truepair){
	  //          std::cout << "the algorithm selected a reco pair matched to a gen pair" << std::endl;
          theHistograms->fill("#Algorithm_" + Candidate.first, "#Algorithm_" + Candidate.first, 2, 0, 2, 1., theWeight*LumiSF);
          theHistograms->fill("PASSED mjj_" + Candidate.first + "_Candidate", " PASSED mjj_" + Candidate.first + "_Candidate", 10, 50, 120, DiJet.mass(), theWeight*LumiSF);
          theHistograms->fill("mjj_" + Candidate.first + "_num", " mjj_" + Candidate.first + "_num; GeV/c^{2}", 10, 50, 120, DiJet.mass(), theWeight*LumiSF);
          theHistograms->fill("Pt_" + Candidate.first + "_num", " Pt_" + Candidate.first + "_num; GeV/c", 10, 0, 300, DiJet.pt(), theWeight*LumiSF);
        }
        else{
          theHistograms->fill("FAILED mjj_" + Candidate.first + "_Candidate", " FAILED mjj_" + Candidate.first + "_Candidate", 10, 50, 120, DiJet.mass(), theWeight*LumiSF);
	  //          std::cout << "the algorithm selected a reco pair NOT matched to a gen pair" << std::endl;
          theHistograms->fill("#Algorithm_" + Candidate.first, "#Algorithm_" + Candidate.first, 2, 0, 2, 0., theWeight*LumiSF);
        }
	/*
        if (truepair && Candidate.first=="mWZ")
        {
          ResolutionPlots(DiJet,Candidate.second,"VectorBosonReconstruction_",theWeight*LumiSF,"");
        }
	*/
    }
  }// end of loop over RECO DiJets matching GEN DiJets
  
  
  phys::Photon mostEnergeticPhoton;

  //  std::stable_sort(selectedphotons.begin(), selectedphotons.end(), phys::EComparator());
  mostEnergeticPhoton = selectedphotons[0];

  //___________________BLOCK_FOR_PhMVAEffSF___________________//
  //weight=theWeight*LumiSF;
  PhEffSF=1.;
  PhEffSFUnc=0.;
  if(theSampleInfo.isMC()){
    PhEffSF    = getPhotonEffSF_MVA(    selectedphotons.at(0), Photon::MVAwp::wp90);
    PhEffSFUnc = getPhotonEffSFUnc_MVA( selectedphotons.at(0), Photon::MVAwp::wp90);
    //weight=theWeight*LumiSF*PhEffSF;
  }
  //__________________________________________________________//


  
  VBTopo=Reconstruct(&recoV,&recoFJ,&haveGoodRECODiJetCand,&haveGoodRECOFJCand,&mostEnergeticPhoton, true, noJERCsys);
  //if(VBTopo==-1) cout<<"FJTopo!"<<endl;

  region="";
  isCR=false;
  //_______SR_______//
  if (IN_GENsignalDef())    printHistos(0, "sign", recoV, recoFJ, selectedphotons, VBTopo, region, isCR);
  else    printHistos(0, "bckg", recoV, recoFJ, selectedphotons, VBTopo, region, isCR);
  printHistos(0, "all", recoV, recoFJ, selectedphotons, VBTopo, region, isCR);   
  /*
  if(!isFSR) printHistos(0, "dib", recoV, recoFJ, selectedphotons, VBTopo, region, isCR);
  else printHistos(0, "fsr", recoV, recoFJ, selectedphotons, VBTopo, region, isCR);     
  */
  if(!promptPhExists) printHistos(0, "nonPrompt", recoV, recoFJ, selectedphotons, VBTopo, region, isCR);
  else printHistos(0, "prompt", recoV, recoFJ, selectedphotons, VBTopo, region, isCR);

  theHistograms->fill(" ######### nbOfGenQuarks vs prompt/nonPro"  , "nbOfGenQuarks vs prompt/nonPro"               , 5, -0.5, 4.5, 2, 0, 2, nbOfGenQuarks, promptPhExists, theWeight);

  //_______CRZOFF_fullSideBand_______//
  isCR=true;
  if(inCRZOFF( recoV, recoFJ, selectedphotons, VBTopo, VZGMVAScore) ){
    region = "CRZOFF";
    if (IN_GENsignalDef())	printHistos(1, "sign", recoV, recoFJ, selectedphotons, VBTopo, region, isCR);   
    else	printHistos(1, "bckg", recoV, recoFJ, selectedphotons, VBTopo, region, isCR);
    printHistos(1, "all", recoV, recoFJ, selectedphotons, VBTopo, region, isCR);   
    /*
    if(!isFSR) printHistos(ANALYSIS_CUTs_WP, "dib", recoV, recoFJ, selectedphotons, VBTopo, region, isCR);
    else printHistos(ANALYSIS_CUTs_WP, "fsr", recoV, recoFJ, selectedphotons, VBTopo, region, isCR);   
    */    
    if(!promptPhExists) printHistos(1, "nonPrompt", recoV, recoFJ, selectedphotons, VBTopo, region, isCR);
    else printHistos(1, "prompt", recoV, recoFJ, selectedphotons, VBTopo, region, isCR);
  }

  //_______CRFSRT_fullLowBand_______//
  if(inCRFSRTight( recoV, recoFJ, selectedphotons, VBTopo, VZGMVAScore) ){
    region = "CRFSRT";
    if (IN_GENsignalDef())	printHistos(ANALYSIS_CUTs_WP, "sign", recoV, recoFJ, selectedphotons, VBTopo, region, isCR);   
    else	printHistos(ANALYSIS_CUTs_WP, "bckg", recoV, recoFJ, selectedphotons, VBTopo, region, isCR);
    printHistos(ANALYSIS_CUTs_WP, "all", recoV, recoFJ, selectedphotons, VBTopo, region, isCR);   
    /*
    if(!isFSR) printHistos(4, "dib", recoV, recoFJ, selectedphotons, VBTopo, region, isCR);
    else printHistos(4, "fsr", recoV, recoFJ, selectedphotons, VBTopo, region, isCR);
    */
    if(!promptPhExists) printHistos(ANALYSIS_CUTs_WP, "nonPrompt", recoV, recoFJ, selectedphotons, VBTopo, region, isCR);
    else printHistos(ANALYSIS_CUTs_WP, "prompt", recoV, recoFJ, selectedphotons, VBTopo, region, isCR);

  }

  //_______ABCD_categorization_______//
  if (inCRZOFF_DIB( recoV, recoFJ, selectedphotons, VBTopo, VZGMVAScore) ) region = "CRZOFF_DIB"; //TOP LEFT
  else if (inCRZOFF_FSRTight( recoV, recoFJ, selectedphotons, VBTopo, VZGMVAScore) ) region = "CRZOFF_FSRT";//BOTTOM LEFT
  else if (inCRZON_FSRTight( recoV, recoFJ, selectedphotons, VBTopo, VZGMVAScore) ) region = "CRZON_FSRT";//BOTTOM RIGHT
  else return;//Central CROSS || SR

  if (IN_GENsignalDef())	printHistos(ANALYSIS_CUTs_WP, "sign", recoV, recoFJ, selectedphotons, VBTopo, region, isCR);   
  else	printHistos(ANALYSIS_CUTs_WP, "bckg", recoV, recoFJ, selectedphotons, VBTopo, region, isCR);
  printHistos(ANALYSIS_CUTs_WP, "all", recoV, recoFJ, selectedphotons, VBTopo, region, isCR);   
  /*
  if(!isFSR) printHistos(ANALYSIS_CUTs_WP, "dib", recoV, recoFJ, selectedphotons, VBTopo, region, isCR);
  else printHistos(ANALYSIS_CUTs_WP, "fsr", recoV, recoFJ, selectedphotons, VBTopo, region, isCR);
  */
  if(!promptPhExists) printHistos(ANALYSIS_CUTs_WP, "nonPrompt", recoV, recoFJ, selectedphotons, VBTopo, region, isCR);
  else printHistos(ANALYSIS_CUTs_WP, "prompt", recoV, recoFJ, selectedphotons, VBTopo, region, isCR);

  
  //_______CRFSRT_categorization_ZON/OFF_______//
  /*
    if(inCRZOFF( recoV, recoFJ, selectedphotons, VBTopo, VZGMVAScore) ) region = "CRZOFF_FSRT";
    else       region = "CRZON_FSRT";
    if (IN_GENsignalDef())	printHistos(ANALYSIS_CUTs_WP, "sign", recoV, recoFJ, selectedphotons, VBTopo, region, isCR);   
    else	printHistos(ANALYSIS_CUTs_WP, "bckg", recoV, recoFJ, selectedphotons, VBTopo, region, isCR);
    printHistos(ANALYSIS_CUTs_WP, "all", recoV, recoFJ, selectedphotons, VBTopo, region, isCR);   
    if(!isFSR) printHistos(ANALYSIS_CUTs_WP, "dib", recoV, recoFJ, selectedphotons, VBTopo, region, isCR);
    else printHistos(ANALYSIS_CUTs_WP, "fsr", recoV, recoFJ, selectedphotons, VBTopo, region, isCR);   
  */
  
}


void VZGAnalyzer::fillFeatTree(FeatList &list, bool &passingPresel )
{
  
  bool isSigSample = theSampleInfo.isMC() && (theSampleInfo.fileName().find("WZG")!=std::string::npos || theSampleInfo.fileName().find("ZZG")!=std::string::npos || theSampleInfo.fileName().find("ZH")!=std::string::npos);
  bool isDYSample = theSampleInfo.isMC() && (theSampleInfo.fileName().find("DY")!=std::string::npos);
  bool isZGSample = theSampleInfo.isMC() && !isSigSample && (theSampleInfo.fileName().find("ZG")!=std::string::npos);
  
  passingPresel = false;
  //if(!IsARunForMVAFeat)  return;
  //  if( (isSigSample && !IN_GENsignalDef() ) || !theSampleInfo.isMC()) return;
  if( isSigSample && !IN_GENsignalDef() ) return;//|| !theSampleInfo.isMC()) return;
  //std::cout<<"0: entering fillFeatTree "<<std::endl;

  int nbOfCutsPassed=0;
  
  int VBTopo = 0;
  phys::Boson<phys::Jet> recoV;
  phys::Jet recoFJ;
  bool haveGoodRECODiJetCand=false;
  bool haveGoodRECOFJCand=false;
  std::vector<phys::Photon> selectedphotons;
  PhotonVLSelection(&selectedphotons,1);

  if(selectedphotons.size()<1) {
    //    list.f_nbOfCutsPassed = 0;  
    return;
  }
  
  phys::Photon mostEnergeticPhoton;

  //  std::stable_sort(selectedphotons.begin(), selectedphotons.end(), phys::EComparator());
  mostEnergeticPhoton = selectedphotons[0];

  if(isDYSample && mostEnergeticPhoton.genStatusFlags().test(phys::isPrompt)) return;
  
  VBTopo=Reconstruct(&recoV,&recoFJ,&haveGoodRECODiJetCand,&haveGoodRECOFJCand,&mostEnergeticPhoton, false, noJERCsys);

  //while(cut(nbOfCutsPassed, recoV, recoFJ, selectedphotons, VBTopo))    nbOfCutsPassed++;
  
  //  list.f_nbOfCutsPassed = nbOfCutsPassed;  

  if(VBTopo!=1) return; 
  //  std::cout<<"2: VBTopo "<<VBTopo<<std::endl;
  double temporaryMVAScore = -2.0;
  if (!cut(2, recoV, recoFJ, selectedphotons, VBTopo, temporaryMVAScore)) return;

  //  std::cout<<"3: passing cuts "<<std::endl;

    //---Block for ZX rwgt---//
  std::string year_str = std::to_string(year);
  double ZXrwgt=1.;
  if(theSampleInfo.isMC() && isDYSample && genVBHelper_.ZtoChLep().size()>0){    
    for(int i = 0; i<DYrewgtBinEdges.size()-1 && ZXrwgt==1.; i++){
      if( genVBHelper_.ZtoChLep()[0].pt() > DYrewgtBinEdges[i] && genVBHelper_.ZtoChLep()[0].pt() < DYrewgtBinEdges[i+1]){
	/*
	if(year==2016){
	  if(year_str.find("preVFP")!=std::string::npos)      ZXrwgt=DYreweights_2016preVFP[i];
	  else ZXrwgt=DYreweights_2016postVFP[i];
	}
	if(year==2017) ZXrwgt=DYreweights_2017[i];
	if(year==2018) ZXrwgt=DYreweights_2018[i];
	*/
	ZXrwgt=DYreweights_Run2[i];
      }
    }
  }

  if(theSampleInfo.isMC() && isZGSample && genVBHelper_.ZtoChLep().size()>0){
    for(int i = 0; i<ZGrewgtBinEdges.size()-1 && ZXrwgt==1.; i++){
      if( genVBHelper_.ZtoChLep()[0].pt() > ZGrewgtBinEdges[i] && genVBHelper_.ZtoChLep()[0].pt() < ZGrewgtBinEdges[i+1]){
	/*
	if(year==2016){
	  if(year_str.find("preVFP")!=std::string::npos)      ZXrwgt=ZGreweights_2016preVFP[i];
	  else ZXrwgt=ZGreweights_2016postVFP[i];
	}
	if(year==2017) ZXrwgt=ZGreweights_2017[i];
	if(year==2018) ZXrwgt=ZGreweights_2018[i];
	*/
	ZXrwgt=ZGreweights_Run2[i];
      }
    }
  }
  //---end of ZX rwgt---//

  
  //---miniblock for qvg sf attempt---//
  double J0qvg= recoV.daughter(0).deepFlavour().probuds/(recoV.daughter(0).deepFlavour().probuds + recoV.daughter(0).deepFlavour().probg);
  double J1qvg= recoV.daughter(1).deepFlavour().probuds/(recoV.daughter(1).deepFlavour().probuds + recoV.daughter(1).deepFlavour().probg);
  int absFlavor=0;
      
  if(isSigSample) absFlavor=1;
  else if(isDYSample) absFlavor=21;
  else if(isZGSample) absFlavor=0;
  else absFlavor=5;
      
  QGScaleFactor J0qvgSF = getQGSF(theSampleInfo.isMC(), J0qvg, recoV.daughter(0).eta(), recoV.daughter(0).pt(), absFlavor, "M");
  QGScaleFactor J1qvgSF = getQGSF(theSampleInfo.isMC(), J1qvg, recoV.daughter(1).eta(), recoV.daughter(1).pt(), absFlavor, "M");
  /*  
  double J0qgTagSF    = J0qvgSF.central;
  double J1qgTagSF    = J1qvgSF.central;
  //  double qgTagSF_up = J0qvgSF.up     *J1qvgSF.up     ;
  //  double qgTagSF_dn = J0qvgSF.down   *J1qvgSF.down   ;    
  */
  //---end of block for qvg sf attempt---//


  TLorentzVector llPh = Z->daughter(0).p4()+Z->daughter(1).p4()+selectedphotons.at(0).p4();
  TLorentzVector l0Ph = Z->daughter(0).p4()+selectedphotons.at(0).p4();
  TLorentzVector l1Ph = Z->daughter(1).p4()+selectedphotons.at(0).p4();

  TLorentzVector jjPh = recoV.daughter(0).p4() + recoV.daughter(1).p4() + selectedphotons.at(0).p4();
  TLorentzVector lljjPh = mostEnergeticPhoton.p4() + Z->daughter(0).p4() + Z->daughter(1).p4() + recoV.daughter(0).p4() + recoV.daughter(1).p4();
  TLorentzVector lljj   = Z->daughter(0).p4() + Z->daughter(1).p4() + recoV.daughter(0).p4() + recoV.daughter(1).p4();

  double mllPh=llPh.M();
  double m2llPh=llPh.M2();
  double m2l0Ph=l0Ph.M2();
  double m2l1Ph=l1Ph.M2();
  double mll=Z->mass();

  //  std::cout<<"4: first vars filled "<<std::endl;

  
  double dPhiZG, dPhiL0G, dPhiL1G, dPhiLL, recoVMass, ptl0, ptl1, FWMT0, ptGamma, ptJ0, ptJ1, etaJ0, etaJ1, etaL0, etaL1, FWMT1, FWMT2, FWMT3, FWMT4, FWMT5, FWMT6, dPhiJ0G, dPhiJ1G, dPhiJJ, dPhiL0J0, dPhiL1J1, dPhiL0J1, dPhiL1J0, deltaR_L0Gamma, deltaR_L1Gamma, deltaR_LL, deltaR_JJ, deltaR_J0Gamma, deltaR_J1Gamma, dRLG,       etaG, mjjG, mlljj, mlljjG, ptll, ptjj, HT;
  int nbOfGoodJets=0;
  int nbOfAllJets=0;
  int phIDpassed=0;

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

  dPhiZG = fabs(physmath::deltaPhi(Z->phi(),selectedphotons.at(0).phi()) );
  
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


  
  deltaR_L0Gamma=fabs(physmath::deltaR(Z->daughter(0), selectedphotons.at(0)));
  deltaR_L1Gamma=fabs(physmath::deltaR(Z->daughter(1), selectedphotons.at(0)));
  if(deltaR_L0Gamma<deltaR_L1Gamma) dRLG=deltaR_L0Gamma;
  else dRLG=deltaR_L1Gamma;
  deltaR_LL=fabs(physmath::deltaR(Z->daughter(0), Z->daughter(1) ) );
  deltaR_JJ=fabs(physmath::deltaR(recoV.daughter(0), recoV.daughter(1) ) );
  deltaR_J0Gamma=fabs(physmath::deltaR(recoV.daughter(0), selectedphotons.at(0)));
  deltaR_J1Gamma=fabs(physmath::deltaR(recoV.daughter(1), selectedphotons.at(0)));
  
  foreach (const phys::Jet &jet, *jets)
    {
    
      if (KinematicsOK(jet,ptcut,etacut) && fabs(physmath::deltaR(jet,selectedphotons[0]))> dR_jetRatio_cut) // KinematicsOK(jet)	
	nbOfGoodJets++;
    }
  
  nbOfAllJets=jets->size();
  //std::cout<<"dRLPH "<<deltaR_LGamma<<std::endl;

  //  if (deltaR_LGamma>10) deltaR_LGamma=11;

  
  FWMT0 =SumFWM(0, 't', lljjG);
  FWMT1 =SumFWM(1, 't', lljjG);
  FWMT2 =SumFWM(2, 't', lljjG);
  FWMT3 =SumFWM(3, 't', lljjG);
  FWMT4 =SumFWM(4, 't', lljjG);
  FWMT5 =SumFWM(5, 't', lljjG);
  FWMT6 =SumFWM(6, 't', lljjG);


  if ( selectedphotons[0].id() == 22 && KinematicsOK(selectedphotons[0], 20, 2.4) && !selectedphotons[0].hasPixelSeed() && selectedphotons[0].passElectronVeto()) phIDpassed=1;
  if ( selectedphotons[0].cutBasedID(Photon::IdWp::VeryLoose ) ) phIDpassed=2;
  if ( selectedphotons[0].passMVA(Photon::MVAwp::wp90        ) )  phIDpassed=3; 
  /*
  if ( selectedphotons[0].cutBasedIDLoose()  ) phIDpassed=3;
  if ( selectedphotons[0].cutBasedIDMedium() ) phIDpassed=4;
  if ( selectedphotons[0].cutBasedIDTight()  ) phIDpassed=5; 
  */   
  
  //p.cutBasedIDLoose()

  if(!theSampleInfo.isMC() && phIDpassed==3) return;
  
  list.f_isData = !theSampleInfo.isMC();
  
  //  std::cout<<"----------------------------------------"<<std::endl;    
  list.f_weight = theWeight;

  list.f_posWeight         = theWeight > 0.    ? theWeight : 0.;
  list.f_weightIfCR2P1F    = phIDpassed == 2   ? theWeight : 0.;
  list.f_posWeightIfCR2P1F = (theWeight > 0. && phIDpassed == 2) ? theWeight : 0.;

  list.f_J0qvgSF = J0qvgSF.central;
  list.f_J1qvgSF = J1qvgSF.central;
  list.f_ZXrw = ZXrwgt;
  
  //  std::cout<<"weight: "<<theWeight<<std::endl;//  theHistograms->fill("the Weight", "the Weight", 50, -5, 5, theWeight, 1);

  list.f_mll  = Z->mass();
  //  std::cout<<"ZCandMass: "<<Z->mass()<<std::endl;//    theHistograms->fill("ZCand mass", "ZCand mass", 30, 60, 120, Z->mass(), 1);

  list.f_ptl1 = Z->daughter(0).pt();
  //  std::cout<<"ptl1: "<<Z->daughter(0).pt()<<std::endl;//      theHistograms->fill("ptl1", "ptl1", 50,  20, 120, Z->daughter(0).pt(), 1);

  list.f_ptl2 = Z->daughter(1).pt();
  //  std::cout<<"ptl2: "<<Z->daughter(1).pt()<<std::endl;//  theHistograms->fill("ptl2", "ptl2", 50,  20, 120, Z->daughter(1).pt(), 1);

  list.f_ptGamma = selectedphotons.at(0).pt();
  //  std::cout<<"ptGamma: "<<selectedphotons.at(0).pt()<<std::endl;//  theHistograms->fill("ptGamma", "ptGamma", 50,  20, 120, selectedphotons.at(0).pt(), 1);

  list.f_etaG = selectedphotons.at(0).eta();
  list.f_ptll = Z->pt();
  list.f_ptjj = recoV.pt();
  list.f_mlljjPh = lljjPh.M();
  list.f_mjjG = jjPh.M();    
  list.f_mlljj = lljj.M();
  
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
  
  list.f_dPhiZG = dPhiZG;
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
  list.f_deltaR_L0Gamma = deltaR_L0Gamma;
  list.f_deltaR_L1Gamma = deltaR_L1Gamma;
  list.f_deltaR_LL = deltaR_LL;
  list.f_deltaR_J0Gamma = deltaR_J0Gamma;
  list.f_deltaR_J1Gamma = deltaR_J1Gamma;
  list.f_deltaR_JJ = deltaR_JJ;
  list.f_mllPh = mllPh;
  list.f_recoVMass=recoVMass;
  //  theHistograms->fill("recoVMass", "recoVMass", 35, 50, 120, recoVMass, 1);

  list.f_FWMT0=FWMT0;
  list.f_FWMT1=FWMT1;
  list.f_FWMT2=FWMT2;
  list.f_FWMT3=FWMT3;
  list.f_FWMT4=FWMT4;
  list.f_FWMT5=FWMT5;
  list.f_FWMT6=FWMT6;
  /*
  theHistograms->fill("FWMT0", "FWMT0", 25, 0, 5, FWMT0, 1);
  theHistograms->fill("FWMT1", "FWMT1", 25, 0, 5, FWMT1, 1);
  theHistograms->fill("FWMT2", "FWMT2", 50, 0, 1, FWMT2, 1);
  theHistograms->fill("FWMT3", "FWMT3", 25, 0, 5, FWMT3, 1);
  theHistograms->fill("FWMT4", "FWMT4", 50, 0, 0.25, FWMT4, 1);
  */

  //  std::cout<<"5: full list filled "<<std::endl;
  //featureTree.Fill();
  //  list.f_nbOfCutsPassed = nbOfCutsPassed;
  list.f_nbOfGoodJets = nbOfGoodJets;
  list.f_nbOfAllJets = nbOfAllJets;

  list.f_PhMVAId=selectedphotons.at(0).MVAvalue();

  list.f_phIDpassed = phIDpassed;//int
  list.f_dRLG = dRLG;  
  list.f_J0DeepProb_b=recoV.daughter(0).deepFlavour().probb;
  list.f_J1DeepProb_b=recoV.daughter(1).deepFlavour().probb;
  list.f_J0DeepProb_c=recoV.daughter(0).deepFlavour().probc;
  list.f_J1DeepProb_c=recoV.daughter(1).deepFlavour().probc;
  list.f_J0DeepProb_g=recoV.daughter(0).deepFlavour().probg;
  list.f_J1DeepProb_g=recoV.daughter(1).deepFlavour().probg;
  list.f_J0DeepProb_lepb=recoV.daughter(0).deepFlavour().problepb;
  list.f_J1DeepProb_lepb=recoV.daughter(1).deepFlavour().problepb;
  list.f_J0DeepProb_uds=recoV.daughter(0).deepFlavour().probuds;
  list.f_J1DeepProb_uds=recoV.daughter(1).deepFlavour().probuds;
  list.f_J0ChMult =recoV.daughter(0).chargedMultiplicity() ;//int
  list.f_J1ChMult =recoV.daughter(1).chargedMultiplicity() ;//int
  list.f_J0NeuMult=recoV.daughter(0).neutralMultiplicity() ;//int
  list.f_J1NeuMult=recoV.daughter(1).neutralMultiplicity() ;//int
  list.f_J0ChEmFrac =recoV.daughter(0).chargedEmEnergyFraction() ;//float
  list.f_J1ChEmFrac =recoV.daughter(1).chargedEmEnergyFraction() ;//float
  list.f_J0NeuEmFrac=recoV.daughter(0).neutralEmEnergyFraction() ;//float
  list.f_J1NeuEmFrac=recoV.daughter(1).neutralEmEnergyFraction() ;//float
  list.f_J0MuFrac=recoV.daughter(0).muonEnergyFraction() ;//float
  list.f_J1MuFrac=recoV.daughter(1).muonEnergyFraction() ;//float
  list.f_J0EleFrac=recoV.daughter(0).electronEnergyFraction() ;//float
  list.f_J1EleFrac=recoV.daughter(1).electronEnergyFraction() ;//float
  list.f_J0PhFrac=recoV.daughter(0).photonEnergyFraction() ;//float
  list.f_J1PhFrac=recoV.daughter(1).photonEnergyFraction() ;//float
  list.f_J0Girth=recoV.daughter(0).girth() ;
  list.f_J1Girth=recoV.daughter(1).girth() ;
  list.f_J0GirthCh=recoV.daughter(0).girth_charged() ;
  list.f_J1GirthCh=recoV.daughter(1).girth_charged() ;
  list.f_J0Area=recoV.daughter(0).jetArea()        ;
  list.f_J1Area=recoV.daughter(1).jetArea()        ;
  list.f_J0Loose=recoV.daughter(0).passLooseJetID()        ;//bool
  list.f_J1Loose=recoV.daughter(1).passLooseJetID()        ;//bool
  list.f_J0QGL=recoV.daughter(0).qgLikelihood()   ;
  list.f_J1QGL=recoV.daughter(1).qgLikelihood()   ;

  list.f_HT=lljjPh.Pt();
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


int VZGAnalyzer::Reconstruct(phys::Boson<phys::Jet> *V_JJCandidate, phys::Jet *V_FJCandidate, bool *haveGoodRECODiJetCand, bool *haveGoodRECOFJCand, phys::Photon *gamma, bool doControlPlots, Systematic jetSystApplied)
{
  int hadrTopo = 0;
  double peakDist_mFJCand=400.;
  double peakDist_mDJCand=400.;
  
  *haveGoodRECOFJCand=false;
  *haveGoodRECODiJetCand=false;
  
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

  bool isGoodEventDiscarded = false;
  bool wasGoodEvent = false;
  bool isBadEventRecovered  = false;

  std::vector<phys::Jet> selectedJets;
  std::vector<phys::Boson<phys::Jet>> DiJetsCand;

  foreach (const phys::Jet &jet, *jets){
    
    if(jetSystApplied.type==SystType::Nominal){
      if (KinematicsOK(jet,ptcut,etacut) && fabs(physmath::deltaR(jet,*gamma))> dR_jetRatio_cut && jet.passLooseJetID()){
	selectedJets.push_back(jet);
      }
    }else if(jetSystApplied.type==SystType::JER && jetSystApplied.direction>0 && KinematicsOKafterSys(jet,jet.ptJerUp(),ptcut,etacut) && fabs(physmath::deltaR(jet,*gamma))> dR_jetRatio_cut && jet.passLooseJetID()){
      selectedJets.push_back(jet);
    }else if(jetSystApplied.type==SystType::JER && jetSystApplied.direction<0 && KinematicsOKafterSys(jet,jet.ptJerDn(),ptcut,etacut) && fabs(physmath::deltaR(jet,*gamma))> dR_jetRatio_cut && jet.passLooseJetID()){
      selectedJets.push_back(jet);
    }else if(jetSystApplied.type==SystType::JES && KinematicsOKafterSys(jet,jet.pt()*(1.+jetSystApplied.direction*getJESUncertainty(jet, jetSystApplied.source, jetSystApplied.yearSys) ),ptcut,etacut) && fabs(physmath::deltaR(jet,*gamma))> dR_jetRatio_cut && jet.passLooseJetID()){
      selectedJets.push_back(jet);
    }
    
    //if (KinematicsOK(jet,ptcut,etacut) && fabs(physmath::deltaR(jet,*gamma))> dR_jetRatio_cut && jet.passLooseJetID()) 	selectedJets.push_back(jet);
      
  }//end of loop over jets


  if (selectedJets.size() > 1){
    for (uint i = 0; i < selectedJets.size() - 1; i++) // Warning: size can be 0
      for (uint j = i+1; j < selectedJets.size(); j++)
        DiJetsCand.push_back(phys::Boson<phys::Jet>(selectedJets.at(i), selectedJets.at(j)));
  }
  //std::cout<<"-->DJ Cand size =" <<DiJetsCand.size()<<endl;
  int updatedBestCandIndex = 0;
  if (DiJetsCand.size() > 0){
    std::stable_sort(DiJetsCand.begin(), DiJetsCand.end(), phys::Mass2Comparator(phys::ZMASS, phys::WMASS));

    if(jetSystApplied.type==SystType::Nominal){
      *haveGoodRECODiJetCand=(DiJetsCand[0].mass()>50  &&   DiJetsCand[0].mass()<120);
    }else{//if(isForSysUpDn!=0){
      double currentPeakDist_mDJCand=999.;

      double ptj0Scaled_JUncUp, ptj0Scaled_JUncDn, ptj0Scaled_JES, ptj1Scaled_JES, ptj1Scaled_JUncUp, ptj1Scaled_JUncDn, mjjScaled_JUncUp, mjjScaled_JUncDn, mjjScaled_JES;
      for(int i=0; i<DiJetsCand.size();i++){
	if(jetSystApplied.type==SystType::JER){
	  ptj0Scaled_JUncUp=DiJetsCand[i].daughter(0).ptJerUp();
	  ptj0Scaled_JUncDn=DiJetsCand[i].daughter(0).ptJerDn();
	  ptj1Scaled_JUncUp=DiJetsCand[i].daughter(1).ptJerUp();
	  ptj1Scaled_JUncDn=DiJetsCand[i].daughter(1).ptJerDn();

	  if(jetSystApplied.direction>0){
	    TLorentzVector JUnc_PJ0_scaling_Up = DiJetsCand[i].daughter(0).p4();
	    TLorentzVector JUnc_PJ1_scaling_Up = DiJetsCand[i].daughter(1).p4();
	    ScaleP4(JUnc_PJ0_scaling_Up, ptj0Scaled_JUncUp, DiJetsCand[i].daughter(0).pt());
	    ScaleP4(JUnc_PJ1_scaling_Up, ptj1Scaled_JUncUp, DiJetsCand[i].daughter(1).pt());
	    mjjScaled_JUncUp=(JUnc_PJ0_scaling_Up+JUnc_PJ1_scaling_Up).M();
	    double candPeakDist_mDJCand=fabs(mjjScaled_JUncUp-phys::ZMASS); 
	    if(fabs(mjjScaled_JUncUp-phys::WMASS)<candPeakDist_mDJCand)    candPeakDist_mDJCand=fabs(mjjScaled_JUncUp-phys::WMASS);
	    if(candPeakDist_mDJCand<currentPeakDist_mDJCand){
	      currentPeakDist_mDJCand=candPeakDist_mDJCand;
	      *haveGoodRECODiJetCand=(mjjScaled_JUncUp>50  && mjjScaled_JUncUp<120); //remove this line to avoid threshold effects
	      updatedBestCandIndex=i;
	      
	      if(!isBadEventRecovered && !wasGoodEvent && i!=0 && mjjScaled_JUncUp>50  && mjjScaled_JUncUp<120) isBadEventRecovered=true;
		
	    }
	  }else if(jetSystApplied.direction<0){
	    TLorentzVector JUnc_PJ0_scaling_Dn = DiJetsCand[i].daughter(0).p4();
	    TLorentzVector JUnc_PJ1_scaling_Dn = DiJetsCand[i].daughter(1).p4();
	    ScaleP4(JUnc_PJ0_scaling_Dn, ptj0Scaled_JUncDn, DiJetsCand[i].daughter(0).pt());
	    ScaleP4(JUnc_PJ1_scaling_Dn, ptj1Scaled_JUncDn, DiJetsCand[i].daughter(1).pt());
	    mjjScaled_JUncDn=(JUnc_PJ0_scaling_Dn+JUnc_PJ1_scaling_Dn).M();
	    double candPeakDist_mDJCand=fabs(mjjScaled_JUncDn-phys::ZMASS);
	    if(fabs(mjjScaled_JUncDn-phys::WMASS)<candPeakDist_mDJCand)    candPeakDist_mDJCand=fabs(mjjScaled_JUncDn-phys::WMASS);
	    if(candPeakDist_mDJCand<currentPeakDist_mDJCand){
	      currentPeakDist_mDJCand=candPeakDist_mDJCand;
	      *haveGoodRECODiJetCand=(mjjScaled_JUncDn>50  && mjjScaled_JUncDn<120); //remove this line to avoid threshold effects
	      updatedBestCandIndex=i;
	    }
	  }
	}else if(jetSystApplied.type==SystType::JES){
	  ptj0Scaled_JES=DiJetsCand[i].daughter(0).pt()*(1.+jetSystApplied.direction*getJESUncertainty(DiJetsCand[i].daughter(0), jetSystApplied.source, jetSystApplied.yearSys) );// (DiJetsCand[i].daughter(0).jesUnc().Total));    
	  ptj1Scaled_JES=DiJetsCand[i].daughter(1).pt()*(1.+jetSystApplied.direction*getJESUncertainty(DiJetsCand[i].daughter(1), jetSystApplied.source, jetSystApplied.yearSys) );// (DiJetsCand[i].daughter(1).jesUnc().Total));
	  TLorentzVector JES_PJ0_scaling = DiJetsCand[i].daughter(0).p4();
	  TLorentzVector JES_PJ1_scaling = DiJetsCand[i].daughter(1).p4();
	  ScaleP4(JES_PJ0_scaling, ptj0Scaled_JES, DiJetsCand[i].daughter(0).pt());
	  ScaleP4(JES_PJ1_scaling, ptj1Scaled_JES, DiJetsCand[i].daughter(1).pt());
	  mjjScaled_JES=(JES_PJ0_scaling+JES_PJ1_scaling).M();
	  double candPeakDist_mDJCand=fabs(mjjScaled_JES-phys::ZMASS);
	  if(fabs(mjjScaled_JES-phys::WMASS)<candPeakDist_mDJCand)    candPeakDist_mDJCand=fabs(mjjScaled_JES-phys::WMASS);
	  if(candPeakDist_mDJCand<currentPeakDist_mDJCand){
	    currentPeakDist_mDJCand=candPeakDist_mDJCand;
	    *haveGoodRECODiJetCand=(mjjScaled_JES>50  && mjjScaled_JES<120); //remove this line to avoid threshold effects
	    updatedBestCandIndex=i;
	  }
	}//closing JES case
      }//closing loop over DiJetCand
    }//closing if not for Sys 
  }//closing if at least one DiJetCand


  if(*haveGoodRECODiJetCand){
    *V_JJCandidate = DiJetsCand.at(updatedBestCandIndex);
    peakDist_mDJCand=fabs(V_JJCandidate->mass()-phys::ZMASS);
    if(fabs(V_JJCandidate->mass()-phys::WMASS)<peakDist_mDJCand)    peakDist_mDJCand=fabs(V_JJCandidate->mass()-phys::WMASS);
    //std::cout<<"peakDist "<<peakDist_mDJCand<<endl;  
  }
  //  if(updatedBestCandIndex!=0)std::cout<<"updated best cand index "<<updatedBestCandIndex<<endl;  

  if(!(*haveGoodRECOFJCand || *haveGoodRECODiJetCand) ){
    hadrTopo=0;
  }else{//changed to prioritize the DJ topology
    if(*haveGoodRECODiJetCand) hadrTopo=1;//    if(peakDist_mDJCand<peakDist_mFJCand) hadrTopo=1;
    else hadrTopo=0;// hadrTopo=-1; // to be modified to consider FJ topo
  }
  /*
  if(*haveGoodRECODiJetCand && doControlPlots){
    double VMassDist_mDJCand=50;
    int DJMax=DiJetsCand.size();
    if(DJMax>5) DJMax = 5;
    if(DJMax>1) theHistograms->fill("dQGL_DJCand10 ", "dQGL_DJCand10 ; dQGL DJCand1 - DJCand0", 20, -1, 1, TMath::Sqrt(DiJetsCand[1].daughter(0).qgLikelihood()*DiJetsCand[1].daughter(1).qgLikelihood()) -  TMath::Sqrt(DiJetsCand[0].daughter(0).qgLikelihood()*DiJetsCand[0].daughter(1).qgLikelihood()), theWeight*LumiSF);

    for(int j = 0; j<DJMax; j++){
      VMassDist_mDJCand=fabs(DiJetsCand[j].mass()-phys::ZMASS);
      if(fabs(DiJetsCand[j].mass()-phys::WMASS)<VMassDist_mDJCand)    VMassDist_mDJCand=fabs(DiJetsCand[j].mass()-phys::WMASS);

      theHistograms->fill("DJCand"+std::to_string(j)+"_j0_QGL_vs_NJ", "DJCand"+std::to_string(j)+"_j1_QGL_vs_NJ; DJCand"+std::to_string(j)+" j0 QGL; NJ", 20, 0, 1, 3, 1.5, 4.5, DiJetsCand[j].daughter(0).qgLikelihood(), selectedJets.size(), theWeight*LumiSF);
      theHistograms->fill("DJCand"+std::to_string(j)+"_j1_QGL_vs_NJ", "DJCand"+std::to_string(j)+"_j1_QGL_vs_NJ; DJCand"+std::to_string(j)+" j0 QGL; NJ", 20, 0, 1, 3, 1.5, 4.5, DiJetsCand[j].daughter(1).qgLikelihood(), selectedJets.size(), theWeight*LumiSF);
      theHistograms->fill("DJCand"+std::to_string(j)+"_jj_QGL_vs_NJ", "DJCand"+std::to_string(j)+"_jj_QGL_vs_NJ; DJCand"+std::to_string(j)+" jj QGL; NJ", 20, 0, 1, 3, 1.5, 4.5, TMath::Sqrt(DiJetsCand[j].daughter(0).qgLikelihood()*DiJetsCand[j].daughter(1).qgLikelihood()), selectedJets.size(), theWeight*LumiSF);
      if(selectedJets.size() > 3){
	theHistograms->fill("DJCand"+std::to_string(j)+"_j0_QGL_vs_NJ", "DJCand"+std::to_string(j)+"_j0_QGL_vs_NJ; DJCand"+std::to_string(j)+" j0 QGL; NJ", 20, 0, 1, 3, 1.5, 4.5, DiJetsCand[j].daughter(0).qgLikelihood(), 4, theWeight*LumiSF);
	theHistograms->fill("DJCand"+std::to_string(j)+"_j1_QGL_vs_NJ", "DJCand"+std::to_string(j)+"_j1_QGL_vs_NJ; DJCand"+std::to_string(j)+" j1 QGL; NJ", 20, 0, 1, 3, 1.5, 4.5, DiJetsCand[j].daughter(1).qgLikelihood(), 4, theWeight*LumiSF);
	theHistograms->fill("DJCand"+std::to_string(j)+"_jj_QGL_vs_NJ", "DJCand"+std::to_string(j)+"_jj_QGL_vs_NJ; DJCand"+std::to_string(j)+" jj QGL; NJ", 20, 0, 1, 3, 1.5, 4.5, TMath::Sqrt(DiJetsCand[j].daughter(0).qgLikelihood()*DiJetsCand[j].daughter(1).qgLikelihood()), 4, theWeight*LumiSF);
      }
      theHistograms->fill("DJCand"+std::to_string(j)+"_jj_QGL_vs_massDist", "DJCand"+std::to_string(j)+"_jj_QGL_vs_massDist; DJCand"+std::to_string(j)+" jj QGL; massDist", 20, 0, 1, 25, 0, 50, TMath::Sqrt(DiJetsCand[j].daughter(0).qgLikelihood()*DiJetsCand[j].daughter(1).qgLikelihood()), VMassDist_mDJCand, theWeight*LumiSF);
    }
  }

  if(*haveGoodRECOFJCand && !*haveGoodRECODiJetCand && doControlPlots){
    theHistograms->fill("FJCand_mass", "FJCand_mass; FJCand_mass", 35, 50, 120, V_FJCandidate->mass(), theWeight*LumiSF);
    theHistograms->fill("FJCand_PNScore_WvsQCD", "FJCand_PNScore_WvsQCD", 20, 0, 1, V_FJCandidate->particleNet().WvsQCD, theWeight*LumiSF);
    theHistograms->fill("FJCand_PNScore_ZvsQCD", "FJCand_PNScore_ZvsQCD", 20, 0, 1, V_FJCandidate->particleNet().ZvsQCD, theWeight*LumiSF);
    theHistograms->fill("FJCand_PNScore_VvsQCD", "FJCand_PNScore_VvsQCD", 20, 0, 1, (fabs(V_FJCandidate->mass()-phys::WMASS)<fabs(V_FJCandidate->mass()-phys::ZMASS) )? V_FJCandidate->particleNet().WvsQCD : V_FJCandidate->particleNet().ZvsQCD, theWeight*LumiSF);

    //    std::cout<<"FJCand_mass ="<< V_FJCandidate->mass() <<endl<<"    PNW="<<V_FJCandidate->particleNet().WvsQCD<<endl<<"    PNZ="<<V_FJCandidate->particleNet().ZvsQCD<<endl<<"    PNV="<<((fabs(V_FJCandidate->mass()-phys::WMASS)<fabs(V_FJCandidate->mass()-phys::ZMASS) )? V_FJCandidate->particleNet().WvsQCD : V_FJCandidate->particleNet().ZvsQCD )<<"   "<<endl;    

  }
  */
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
  //if(IsARunForMVAFeat)
  //hadrTopo=*haveGoodRECODiJetCand;//temporary unique topology

  //if(hadrTopo== -1) cout<<"FJ topo!!!"<<endl;
  //std::cout<<"Reconstruction returning topo"<<hadrTopo<<endl;
  return hadrTopo;
}

double VZGAnalyzer::VHadScore(phys::Boson<phys::Jet> DJCand, int algoType, double inflecPt, double smoothness){
  double dMMax=45.0;
  double dM=dMMax;
  dM=fabs(DJCand.mass()-phys::ZMASS);
  if(fabs(DJCand.mass()-phys::WMASS)<dM)    dM=fabs(DJCand.mass()-phys::WMASS);

  double dInvM=1.-(dM/dMMax);
  if(algoType>0) return dInvM; 

  double QGLV= TMath::Sqrt(DJCand.daughter(0).qgLikelihood()*DJCand.daughter(1).qgLikelihood());
  if(algoType<0) return QGLV; 

  double wgtQGLV=Sigmoid(dM/dMMax,inflecPt,smoothness);
  return wgtQGLV*QGLV+(1.-wgtQGLV)*dInvM;
    
}

int VZGAnalyzer::ReconstructAlt(phys::Boson<phys::Jet> *V_JJCandidate, phys::Jet *V_FJCandidate, bool *haveGoodRECODiJetCand, bool *haveGoodRECOFJCand, phys::Photon *gamma, bool doControlPlots, double inflecPt, double smoothness)
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
    if(DiJetsCand.size()>1){
      double BestVHadScore=0;
      double tempVHadScore;
      int bestCandIndex =0;
      for(int j = 0; j<DiJetsCand.size(); j++){
	tempVHadScore=VHadScore(DiJetsCand[j], 0, inflecPt, smoothness);
	if(tempVHadScore > BestVHadScore){
	  BestVHadScore=tempVHadScore;
	  bestCandIndex=j;
	}
      }
      *V_JJCandidate = DiJetsCand.at(bestCandIndex);
    }
  }

  *haveGoodRECODiJetCand=(V_JJCandidate->mass()>50  &&   V_JJCandidate->mass()<120);

  if(!(*haveGoodRECOFJCand || *haveGoodRECODiJetCand) ){
    hadrTopo=0;
  }else{//changed to consider the only 2j topology
    if(*haveGoodRECODiJetCand) hadrTopo=1;//    if(peakDist_mDJCand<peakDist_mFJCand) hadrTopo=1;
    else hadrTopo=0;//    else hadrTopo=-1;
  }

  hadrTopo=*haveGoodRECODiJetCand;//temporary unique topology
  
  return hadrTopo;
}


void VZGAnalyzer::PhotonSelection(std::vector<phys::Photon> *phot)
{
  /*
  bool looseGammaExists=false;
  //  bool VLGammaExists=false;

  for (auto p : *photons){
    // Pixel seed and electron veto
    // if (ph.hasPixelSeed() || !ph.passElectronVeto())
    //        continue;
    if(!looseGammaExists && p.id() == 22 && KinematicsOK(p, 20, 2.4) && !p.hasPixelSeed() && p.passElectronVeto() && p.cutBasedIDLoose()){
      looseGammaExists=true;
    }

  }
  for (auto p : *photons){
    // Pixel seed and electron veto
    // if (ph.hasPixelSeed() || !ph.passElectronVeto())
    //        continue;
    if(!looseGammaExists && p.id() == 22 && KinematicsOK(p, 20, 2.4) && !p.hasPixelSeed() && p.passElectronVeto() && p.cutBasedID(Photon::IdWp::VeryLoose)){
      phot->push_back(p);
    }

  }
  */
  std::vector<phys::Photon> gamma;
  phys::Photon tightestGamma;
  //  double tightestGammaMVAvalue = -1.1;
  std::vector<phys::Particle> fsrRecovered;

  std::bitset<2> fsrIndex = std::bitset<2>(Z->daughtersWithFSR());
  if(fsrIndex.test(0))    fsrRecovered.push_back(Z->fsrPhoton(0));
  if(fsrIndex.test(1))    fsrRecovered.push_back(Z->fsrPhoton(1));
    
  for (auto p : *photons){
    bool isFsrPh=false;
    if (p.id() == 22 && KinematicsOK(p, 20, 2.4) && !p.hasPixelSeed() && p.passElectronVeto() && p.passMVA(Photon::MVAwp::wp90)){//photon would pass signal selection
      if(fsrRecovered.size()<1) phot->push_back(p);
      else{
	for(auto fsr : fsrRecovered){
	  if( fabs(physmath::deltaR(fsr, p)) < 0.001 ) {
	    isFsrPh=true;
	    break;
	  }
	}
	if(isFsrPh) continue;
	else{
	  phot->push_back(p);
	}
      }
    }
  }

	
  /*
  if(gamma.size()<1) return;

  tightestGamma = *std::max_element(gamma.begin(), gamma.end(),
				    [](const phys::Photon PhA, const phys::Photon PhB) {
				      return PhA.MVAvalue() < PhB.MVAvalue();
				    });
  //  std::cout<<"PHOTON SELECTED, MVA ID = "<<tightestGamma.MVAvalue()<<endl;
  phot->push_back(tightestGamma);
  */
  if(phot->size()>0)    std::stable_sort(phot->begin(), phot->end(), phys::EComparator());
  
  //std::cout << "Number of selected RECO photons = " << phot->size() << std::endl;
}

void VZGAnalyzer::PhotonVLSelection(std::vector<phys::Photon> *phot, int cutBasedWP)
{//Nothing = 0 //Kin acc = 1 //at least VL = 2//at least Loose = 3//at least Medium = 4//Tight = 5
  //  std::cout<<"entering PhotonVLSelection"<<endl;
  std::vector<phys::Photon> kinGamma, VLGamma, looseGamma, mediumGamma, tightGamma;
  int idCutsPassed=0;//Nothing = 0 //Kin!VL = 1 //VL!Loose = 2//Loose!Medium = 3//Medium!Tight = 4//Tight = 5

  std::vector<phys::Particle> fsrRecovered;

  std::bitset<2> fsrIndex = std::bitset<2>(Z->daughtersWithFSR());

  if(fsrIndex.test(0))    fsrRecovered.push_back(Z->fsrPhoton(0));
  if(fsrIndex.test(1))    fsrRecovered.push_back(Z->fsrPhoton(1));
  
  foreach (auto p , *photons){
    bool isFsrPh=false;
    if (p.id() == 22 && KinematicsOK(p, 20, 2.4) && !p.hasPixelSeed() && p.passElectronVeto()){
      if(fsrRecovered.size()>0){
	for(auto fsr : fsrRecovered){
	  if( fabs(physmath::deltaR(fsr, p)) < 0.001 ) {
	    isFsrPh=true;
	    break;
	  }
	}
	if(isFsrPh) continue;
      }
      kinGamma.push_back(p);
      if (p.cutBasedID(Photon::IdWp::VeryLoose)) {
	VLGamma.push_back(p);
	if (p.cutBasedIDLoose()) {
	  looseGamma.push_back(p);
	  if (p.cutBasedIDMedium()) {
	    mediumGamma.push_back(p);
	    if (p.cutBasedIDTight()) {
	      tightGamma.push_back(p);
	    }
	  }
	}
      }
    }
  }

  if(tightGamma.size()>0){
    std::stable_sort(tightGamma.begin(), tightGamma.end(), phys::EComparator());
    phot->push_back(tightGamma.at(0));
    return;
  }
  if(mediumGamma.size()>0 && cutBasedWP<=4){
    std::stable_sort(mediumGamma.begin(), mediumGamma.end(), phys::EComparator());
    phot->push_back(mediumGamma.at(0));
    return;
  }
  if(looseGamma.size()>0 && cutBasedWP<=3){
    std::stable_sort(looseGamma.begin(), looseGamma.end(), phys::EComparator());
    phot->push_back(looseGamma.at(0));
    return;
  }
  if(VLGamma.size()>0 && cutBasedWP<=2){
    std::stable_sort(VLGamma.begin(), VLGamma.end(), phys::EComparator());
    phot->push_back(VLGamma.at(0));
    return;
  }
  if(kinGamma.size()>0 && cutBasedWP==1){
    std::stable_sort(kinGamma.begin(), kinGamma.end(), phys::EComparator());
    phot->push_back(kinGamma.at(0));
  }
  //  std::cout<<"exiting PhotonVLSelection"<<endl;
}



void VZGAnalyzer::CompatibilityTest(phys::Boson<phys::Jet> bestCandidate, phys::Boson<phys::Particle> genVB, std::string sample, std::string algorithm)
{
}

void VZGAnalyzer::printHistos(uint i, std::string histoType, phys::Boson<phys::Jet> recoV, phys::Jet recoFJ, std::vector<phys::Photon> selectedphotons, int VBTopo,  std::string region, bool isCR)
{
  if(verboseControlBlinding)
    std::cout<<"In printHistos: region: "<<region<<"          after cut"<<i<<"           histo type:  "<<histoType<<"            ph eff SF = "<<PhEffSF<<endl;
  /*
  if(i>1){
    if(i==2) std::cout<<"------------------------------------------------------"<<endl;
    std::cout<<"In printHistos: region: "<<region<<"          after cut"<<i<<"           histo type:  "<<histoType<<"            ph eff SF = "<<PhEffSF<<endl;
  }
  */
  //  if(VBTopo==-1) cout<<"FJTopo! --> entering printHistos"<<endl;
  bool isSigSample = theSampleInfo.isMC() && (theSampleInfo.fileName().find("WZG")!=std::string::npos || theSampleInfo.fileName().find("ZZG")!=std::string::npos || theSampleInfo.fileName().find("ZH")!=std::string::npos);
  bool isDYSample = theSampleInfo.isMC() && (theSampleInfo.fileName().find("DY")!=std::string::npos);
  bool isZGSample = theSampleInfo.isMC() && !isSigSample && (theSampleInfo.fileName().find("ZG")!=std::string::npos);

  if(i==0){
    //std::cout<<"In printHistos: region: "<<region<<"          after cut"<<i<<"           histo type:  "<<histoType<<endl;
    theHistograms->fill("#AAA_cut_flow_" + histoType, "Cut flow", cutsToApply, 0, cutsToApply, i, (theWeight*LumiSF));
    theHistograms->fill("#AAA_unw_cut_flow_" + histoType, "Unw. events cut flow", cutsToApply, 0, cutsToApply, i, 1.);      
    theHistograms->fill("photonID_" + histoType +"_noReq", "photonID", 4, 0, 4, 0, (theWeight*LumiSF));
    theHistograms->fill("photonID_" + histoType +"_DJtopo", "photonID", 4, 0, 4, 0, (theWeight*LumiSF));
  }
  
  if(selectedphotons.size()<1) return;
  
  //  std::cout<<"passing cut 0 in printHistos"<<endl;
  double mimicVZGMVAScore = -2.;
  double VZGMVAScore      = -2.;   

  
  if(i==1){
    //    std::cout<<"In printHistos: region: "<<region<<"          after cut"<<i<<"           histo type:  "<<histoType<<endl;
    theHistograms->fill("#AAA_cut_flow_" + histoType, "Cut flow", cutsToApply, 0, cutsToApply, i, (theWeight*PhEffSF*LumiSF));
    theHistograms->fill("#AAA_unw_cut_flow_" + histoType, "Unw. events cut flow", cutsToApply, 0, cutsToApply, i, 1.);      
    if(selectedphotons.at(0).cutBasedIDLoose())          theHistograms->fill("photonID_" + histoType +"_noReq", "photonID", 4, 0, 4, 1, (theWeight*LumiSF));
    if(selectedphotons.at(0).cutBasedIDMedium())     theHistograms->fill("photonID_" + histoType +"_noReq", "photonID", 4, 0, 4, 2, (theWeight*LumiSF));
    if(selectedphotons.at(0).cutBasedIDTight())     theHistograms->fill("photonID_" + histoType+"_noReq", "photonID", 4, 0, 4, 3, (theWeight*LumiSF));

    foreach (const phys::Jet &jet, *jets){
      if (KinematicsOK(jet,0.,4.7) && jet.passLooseJetID())
	theHistograms->fill("AUX_Jets_pt_presel" ,"AUX_Jets_pt_presel", 40, 0, 200, jet.pt(), theWeight*PhEffSF*LumiSF);
    }
  }

  if(!UNBLIND && !isCR && !theSampleInfo.isMC()) return; //SR BLINDING
  if(isCR){
    if(i>ANALYSIS_CUTs_WP) return; //to avoid producing CR plots for cut higher than 
    //i=ANALYSIS_CUTs_WP;
    if(!isRunForCR) i=ANALYSIS_CUTs_WP; // CT: otherwise no CRDY SYS plots would be produced
    histoType = histoType+"_"+region;
  }//Note: this is not active for the sys plots, bc the isCR bool is passed false as an argument when calling printHistos for the SR. Needs to be re-thought for the unblinding step  
  /*
  if(i==1 && cut(1, recoV, recoFJ, selectedphotons, VBTopo, mimicVZGMVAScore) && VBTopo==1){// && LGsolved){
    VZGMVAScore      = VZGMVAScoreEval(recoV,       recoFJ,  selectedphotons, VBTopo,       0, 0, false, 0);
  }
  */
  bool isForSys = (theSampleInfo.isMC()
		       && (
			   (isDYSample && histoType=="nonPrompt")
			   ||
			   (isZGSample && histoType=="prompt")
			   ||
			   (isSigSample && histoType=="sign")
			   ||
			   (!isSigSample && !isZGSample && !isDYSample && histoType=="all")
			   )
		   );
  if(isRunForCR) isForSys = (theSampleInfo.isMC()
			       && (
				   (isDYSample && histoType=="nonPrompt_"+REGION_FOR_FIT)
				   ||
				   (isZGSample && histoType=="prompt_"+REGION_FOR_FIT)
				   ||
				   (isSigSample && histoType=="sign_"+REGION_FOR_FIT)
				   ||
				   (!isSigSample && !isZGSample && !isDYSample && histoType=="all_"+REGION_FOR_FIT)
				   )
			       ); 
  bool isDYPro = isDYSample && histoType=="prompt";
  //isForSys=false;//turn on here for running only CR plots

  //double mVZGcand = ( selectedphotons[0].p4() + Z->p4() + recoV.p4() ).M();
  
  //_________________________________________________BLOCK_FOR_SYS_HISTOS________________________________________//
  /*
  bool LGsolved=false;
  if(i==1)   LGsolved = fabs(physmath::deltaR(Z->daughter(0), selectedphotons.at(0)))>0.5 && fabs(physmath::deltaR(Z->daughter(1), selectedphotons.at(0)) )>0.5;
  */
  if(i==1  && !isCR && !isRunForCR  && (UNBLIND && !theSampleInfo.isMC() && histoType=="all")  && cut(1, recoV, recoFJ, selectedphotons, VBTopo, mimicVZGMVAScore) && VBTopo==1){// && LGsolved){
    if(verboseControlBlinding){
      std::  cout << "----------------------------------------" << event << endl;
      std::  cout << "Run: " << run << " event: " << event << endl;
      std::cout<<"unblinding VZGMVAScore for data"<<endl;
    }
    VZGMVAScore      = VZGMVAScoreEval(recoV,       recoFJ,  selectedphotons, VBTopo,       noJERCsys, false, 0);
    if(VZGMVAScore <= 1. && VZGMVAScore >= binEdges.at(0) ){
      theHistograms->fill("SYS_BDTScore_central", "SYS_BDTScore_central" , binEdges,  VZGMVAScore, theWeight*PhEffSF*LumiSF);
    }
    
    if(VZGMVAScore >= BDT_CUT && VZGMVAScore <= 1. && Z->mass()>80){
      theHistograms->fill("SYS_mll_central", "SYS_mll_central" , mll_binEdges,  Z->mass(), theWeight*PhEffSF*LumiSF);
    }
    
  }

  if(i==1 && !isCR && !isRunForCR  && isForSys && cut(1, recoV, recoFJ, selectedphotons, VBTopo, mimicVZGMVAScore)){// && LGsolved){

    double VZGMVAScore_DFSup= -2.;
    double VZGMVAScore_DFSdn= -2.;
    
    if(VBTopo==1){
      VZGMVAScore       = VZGMVAScoreEval(recoV,       recoFJ,  selectedphotons, VBTopo, noJERCsys,  false, 0);
      VZGMVAScore_DFSup = VZGMVAScoreEval(recoV,       recoFJ,  selectedphotons, VBTopo, noJERCsys,  true, +1);
      VZGMVAScore_DFSdn = VZGMVAScoreEval(recoV,       recoFJ,  selectedphotons, VBTopo, noJERCsys,  true, -1);
    }
    
    // QCD scale
    // envelope: consider the six variations: {Do, Central, Up} x {Dn, Central, Up} - (central, central) - (Dn, Dn) - (Up, Up) and use the max and min
    float QCDscale_Up(1.), QCDscale_Dn(1.);
    std::vector<float> envelope {
      theSampleInfo.QCDscale_muR0p5F1(),
      theSampleInfo.QCDscale_muR0p5F2(),
      theSampleInfo.QCDscale_muR1F0p5(),
      theSampleInfo.QCDscale_muR1F2(),
      theSampleInfo.QCDscale_muR2F0p5(),
      theSampleInfo.QCDscale_muR2F1()
    };
    QCDscale_Up = *max_element(envelope.begin(), envelope.end());
    QCDscale_Dn = *min_element(envelope.begin(), envelope.end());

    //---miniblock for qvg sf attempt---//
    double J0qvg= recoV.daughter(0).deepFlavour().probuds/(recoV.daughter(0).deepFlavour().probuds + recoV.daughter(0).deepFlavour().probg);
    double J1qvg= recoV.daughter(1).deepFlavour().probuds/(recoV.daughter(1).deepFlavour().probuds + recoV.daughter(1).deepFlavour().probg);
    int absFlavor=0;

    double DFqJ0 = recoV.daughter(0).deepFlavour().probuds;
    double DFgJ1 = recoV.daughter(1).deepFlavour().probg;

    theHistograms->fill("J0qvg raw - btv", "J0qvg raw - btv" , 40, -1, 1,  DFqJ0 - J0qvg, theWeight*rewgt*PhEffSF*LumiSF);

    
    if(isSigSample) absFlavor=1;
    else if(isDYSample) absFlavor=21;
    else if(isZGSample) absFlavor=0;
    else absFlavor=5;
    
    QGScaleFactor J0qvgSF = getQGSF(theSampleInfo.isMC(), J0qvg, recoV.daughter(0).eta(), recoV.daughter(0).pt(), absFlavor, "M");
    QGScaleFactor J1qvgSF = getQGSF(theSampleInfo.isMC(), J1qvg, recoV.daughter(1).eta(), recoV.daughter(1).pt(), absFlavor, "M");
    
    double qgTagSF    = J0qvgSF.central*J1qvgSF.central;
    double qgTagSF_up = J0qvgSF.up     *J1qvgSF.up     ;
    double qgTagSF_dn = J0qvgSF.down   *J1qvgSF.down   ;    
    /*
    std::cout<<"Abs flavor ="        <<absFlavor <<endl;
    std::cout<<"Tot qvg SF central ="<<qgTagSF   <<endl;
    std::cout<<"Tot qvg SF up ="     <<qgTagSF_up<<endl;
    std::cout<<"Tot qvg SF down ="   <<qgTagSF_dn<<endl;

    std::cout<<"-------------------------------------"<<endl;
    */
    //---end of block for qvg sf attempt---//

    
    if(VZGMVAScore <= 1. && VZGMVAScore >= binEdges.at(0)){
      theHistograms->fill("SYS_BDTScore_central", "SYS_BDTScore_central" , binEdges,  VZGMVAScore, theWeight*rewgt*qgTagSF*PhEffSF*LumiSF);
      theHistograms->fill("AUX_BDTScore_centralUnw", "AUX_BDTScore_centralUnw" , binEdges,  VZGMVAScore, 1.);
      theHistograms->fill("AUX_weights" , 50,-5,5,  theWeight*rewgt*qgTagSF*PhEffSF, 1.);

      if(isSigSample){//CT: cap applied
	theHistograms->fill("SYS_BDTScore_alphas_Up"  , "SYS_BDTScore_alphas_Up"   , binEdges,  VZGMVAScore,
			         theWeight*rewgt*qgTagSF*PhEffSF*LumiSF*(
                 		          fabs( theSampleInfo.alphas_MZ_Up()  -1.)  < ALPHAS_CAP ?
									           theSampleInfo.alphas_MZ_Up() :
							                           ( 1.+ALPHAS_CAP*fabs(theSampleInfo.alphas_MZ_Up())  /theSampleInfo.alphas_MZ_Up())
			         )
			    );
	theHistograms->fill("SYS_BDTScore_alphas_Down", "SYS_BDTScore_alphas_Down" , binEdges,  VZGMVAScore,
			          theWeight*rewgt*qgTagSF*PhEffSF*LumiSF*(
					  fabs( theSampleInfo.alphas_MZ_Down()-1.)  < ALPHAS_CAP ?
										   2-theSampleInfo.alphas_MZ_Down()  :
					                                           2-( 1.+ALPHAS_CAP*fabs(theSampleInfo.alphas_MZ_Down())/theSampleInfo.alphas_MZ_Down())
					  )
			    );

      }else{//CT: original implementation
	theHistograms->fill("SYS_BDTScore_alphas_Up"  , "SYS_BDTScore_alphas_Up"   , binEdges,  VZGMVAScore, theSampleInfo.alphas_MZ_Up()*theWeight*rewgt*qgTagSF*PhEffSF*LumiSF);
	//	if(isSigSample || isZGSample || isDYSample) //CT: original implementation
	if(isZGSample || isDYSample)
	  theHistograms->fill("SYS_BDTScore_alphas_Down", "SYS_BDTScore_alphas_Down" , binEdges,  VZGMVAScore, (2.-theSampleInfo.alphas_MZ_Down())*theWeight*rewgt*qgTagSF*PhEffSF*LumiSF);
	else
	  theHistograms->fill("SYS_BDTScore_alphas_Down", "SYS_BDTScore_alphas_Down" , binEdges,  VZGMVAScore, theSampleInfo.alphas_MZ_Down()*theWeight*rewgt*qgTagSF*PhEffSF*LumiSF);
      }
      /*
      if( (fabs( theSampleInfo.alphas_MZ_Up()   -1.)> 1.5 ||  fabs(theSampleInfo.alphas_MZ_Down()   -1.)> 1.5) && fabs(theWeight*rewgt*qgTagSF*PhEffSF)>0.01 ){
	std::  cout << "----------------------------------------"<<endl<<"Event nb. "<< event << endl;
	std::  cout << "----------------------------------------"<<endl;
	std::cout<<"event weight: "<<theWeight*rewgt*qgTagSF*PhEffSF<<endl;
	std::cout<<"aS up       : "<<theSampleInfo.alphas_MZ_Up()<<endl;
	std::cout<<"2 - aS down : "<<2.-theSampleInfo.alphas_MZ_Down()<<endl;
	//	std::cout<<"aS up       : "<<theSampleInfo.alphas_MZ_Up()<<endl;
	std::cout<<"1/aS down   : "<<1/theSampleInfo.alphas_MZ_Down()<<endl;
      }
      *//*
      theHistograms->fill("AlphaS up variation (nominal weights)"  , "AlphaS up variation (nominal weights)"   , 22, -1.2, 3.2, fabs( theSampleInfo.alphas_MZ_Up()  -1.)  < 2. ? theSampleInfo.alphas_MZ_Up()      : 1.+2.1*fabs(theSampleInfo.alphas_MZ_Up())  /theSampleInfo.alphas_MZ_Up()  , theWeight*rewgt*qgTagSF*PhEffSF*LumiSF);
      theHistograms->fill("AlphaS 2-dn variation (nominal weights)", "AlphaS 2-dn variation (nominal weights)" , 22, -1.2, 3.2, fabs( theSampleInfo.alphas_MZ_Down()-1.)  < 2. ? 2.-theSampleInfo.alphas_MZ_Down()  : 1.-2.1*fabs(theSampleInfo.alphas_MZ_Down())/theSampleInfo.alphas_MZ_Down(),theWeight*rewgt*qgTagSF*PhEffSF*LumiSF);

      theHistograms->fill("AlphaS up variation (variation applied)"  , "AlphaS up variation (variation applied)"   , 22, -1.2, 3.2, fabs( theSampleInfo.alphas_MZ_Up()  -1.)  < 2. ? theSampleInfo.alphas_MZ_Up()     : 1.+2.1*fabs(theSampleInfo.alphas_MZ_Up())  /theSampleInfo.alphas_MZ_Up()   , theSampleInfo.alphas_MZ_Up()*theWeight*rewgt*qgTagSF*PhEffSF*LumiSF);
      theHistograms->fill("AlphaS 2-dn variation (variation applied)", "AlphaS 2-dn variation (variation applied)" , 22, -1.2, 3.2, fabs( theSampleInfo.alphas_MZ_Down()-1.)  < 2. ? 2.-theSampleInfo.alphas_MZ_Down() : 1.-2.1*fabs(theSampleInfo.alphas_MZ_Down())/theSampleInfo.alphas_MZ_Down(), (2.-theSampleInfo.alphas_MZ_Down())*theWeight*rewgt*qgTagSF*PhEffSF*LumiSF);

      theHistograms->fill("AlphaS up variation (w/cap, after variation)"  , "AlphaS up variation (w/cap, after variation)"   , 10.*ALPHAS_CAP+2., 1.-ALPHAS_CAP-0.2, 1.+ALPHAS_CAP+0.2,
			  fabs( theSampleInfo.alphas_MZ_Up()  -1.)  < ALPHAS_CAP ?
			                                                           theSampleInfo.alphas_MZ_Up()  :
			                                                           1.+(ALPHAS_CAP+.1)*fabs(theSampleInfo.alphas_MZ_Up())  /theSampleInfo.alphas_MZ_Up()  ,
			  theWeight*rewgt*qgTagSF*PhEffSF*LumiSF*( fabs( theSampleInfo.alphas_MZ_Up()  -1.)  < ALPHAS_CAP ?
										   theSampleInfo.alphas_MZ_Up() :
							                           ( 1.+(ALPHAS_CAP-.0001)*fabs(theSampleInfo.alphas_MZ_Up())  /theSampleInfo.alphas_MZ_Up())
							   )
			  );
      theHistograms->fill("AlphaS 2-dn variation (w/cap, after variation)", "AlphaS 2-dn variation (w/cap, after variation)" , 10.*ALPHAS_CAP+2., 1.-ALPHAS_CAP-0.2, 1.+ALPHAS_CAP+0.2,
                          fabs( theSampleInfo.alphas_MZ_Down()-1.)  < ALPHAS_CAP ?
			                                                           2-theSampleInfo.alphas_MZ_Down()  :
			                                                           1.-(ALPHAS_CAP+.1)*fabs(theSampleInfo.alphas_MZ_Down())/theSampleInfo.alphas_MZ_Down()  ,
			  theWeight*rewgt*qgTagSF*PhEffSF*LumiSF*( fabs( theSampleInfo.alphas_MZ_Down()-1.)  < ALPHAS_CAP ?
										   2-theSampleInfo.alphas_MZ_Down()  :
							                           2-( 1.+(ALPHAS_CAP-.0001)*fabs(theSampleInfo.alphas_MZ_Down())/theSampleInfo.alphas_MZ_Down())
							   )
			  );
      */
      theHistograms->fill("SYS_BDTScore_PDFVar_Up"  , "SYS_BDTScore_PDFVar_Up"   , binEdges,  VZGMVAScore, theSampleInfo.PDFVar_Up()*theWeight*rewgt*qgTagSF*PhEffSF*LumiSF);
      theHistograms->fill("SYS_BDTScore_PDFVar_Down", "SYS_BDTScore_PDFVar_Down" , binEdges,  VZGMVAScore, theSampleInfo.PDFVar_Down()*theWeight*rewgt*qgTagSF*PhEffSF*LumiSF);

      theHistograms->fill("SYS_BDTScore_QCDscale_Up"  , "SYS_BDTScore_QCDscale_Up"   , binEdges,  VZGMVAScore, QCDscale_Up*theWeight*rewgt*qgTagSF*PhEffSF*LumiSF);
      theHistograms->fill("SYS_BDTScore_QCDscale_Down", "SYS_BDTScore_QCDscale_Down" , binEdges,  VZGMVAScore, QCDscale_Dn*theWeight*rewgt*qgTagSF*PhEffSF*LumiSF);

      theHistograms->fill("SYS_BDTScore_DeepFlavorQGModeling_Up"  , "SYS_BDTScore_DeepFlavorQGModeling_Up"   , binEdges,  VZGMVAScore, theWeight*rewgt*qgTagSF_up*PhEffSF*LumiSF);
      theHistograms->fill("SYS_BDTScore_DeepFlavorQGModeling_Down", "SYS_BDTScore_DeepFlavorQGModeling_Down" , binEdges,  VZGMVAScore, theWeight*rewgt*qgTagSF_dn*PhEffSF*LumiSF);

      
      theHistograms->fill("SYS_BDTScore_L1Prefiring_Up"  , "SYS_BDTScore_L1Prefiring_Up"   , binEdges,  VZGMVAScore, (theSampleInfo.L1PrefiringWeightUp()/theSampleInfo.L1PrefiringWeight())*theWeight*rewgt*qgTagSF*PhEffSF*LumiSF);
      theHistograms->fill("SYS_BDTScore_L1Prefiring_Down", "SYS_BDTScore_L1Prefiring_Down" , binEdges,  VZGMVAScore, (theSampleInfo.L1PrefiringWeightDn()/theSampleInfo.L1PrefiringWeight())*theWeight*rewgt*qgTagSF*PhEffSF*LumiSF);

      double relPhEffSFUnc=(PhEffSF!=0) ? PhEffSFUnc/PhEffSF : 0.;
      theHistograms->fill("SYS_BDTScore_effPhIDMVA_Up"  , "SYS_BDTScore_effPhIDMVA_Up"   , binEdges,  VZGMVAScore, theWeight*rewgt*qgTagSF*(PhEffSF+PhEffSFUnc)*LumiSF);
      theHistograms->fill("SYS_BDTScore_effPhIDMVA_Down", "SYS_BDTScore_effPhIDMVA_Down" , binEdges,  VZGMVAScore, theWeight*rewgt*qgTagSF*(PhEffSF-PhEffSFUnc)*LumiSF);

      double eleEff_w=0., muoEff_w=0.;
      eleEff_w  = Z->eleEffSFUnc()/Z->efficiencySF();
      muoEff_w  = Z->muoEffSFUnc()/Z->efficiencySF();
      theHistograms->fill("SYS_BDTScore_electronEff_Up"  , "SYS_BDTScore_electronEff_Up"   , binEdges,  VZGMVAScore, theWeight*rewgt*qgTagSF*(1 + eleEff_w)*LumiSF);
      theHistograms->fill("SYS_BDTScore_electronEff_Down", "SYS_BDTScore_electronEff_Down" , binEdges,  VZGMVAScore, theWeight*rewgt*qgTagSF*(1 - eleEff_w)*LumiSF);
      theHistograms->fill("SYS_BDTScore_muonEff_Up"  , "SYS_BDTScore_muonEff_Up"   , binEdges,  VZGMVAScore, theWeight*rewgt*qgTagSF*(1 + muoEff_w)*LumiSF);
      theHistograms->fill("SYS_BDTScore_muonEff_Down", "SYS_BDTScore_muonEff_Down" , binEdges,  VZGMVAScore, theWeight*rewgt*qgTagSF*(1 - muoEff_w)*LumiSF);

      theHistograms->fill("SYS_BDTScore_puWeight_Up"  , "SYS_BDTScore_puWeight_Up"   , binEdges,  VZGMVAScore, (theSampleInfo.puWeightUncUp()/theSampleInfo.puWeight())*theWeight*rewgt*qgTagSF*PhEffSF*LumiSF);
      theHistograms->fill("SYS_BDTScore_puWeight_Down", "SYS_BDTScore_puWeight_Down" , binEdges,  VZGMVAScore, (theSampleInfo.puWeightUncDn()/theSampleInfo.puWeight())*theWeight*rewgt*qgTagSF*PhEffSF*LumiSF);
      
      theHistograms->fill("SYS_BDTScore_ZXreweighting_Up"  , "SYS_BDTScore_ZXreweighting_Up"   , binEdges,  VZGMVAScore, theWeight*qgTagSF*(rewgt+rewErr)*LumiSF);
      theHistograms->fill("SYS_BDTScore_ZXreweighting_Down", "SYS_BDTScore_ZXreweighting_Down" , binEdges,  VZGMVAScore, theWeight*qgTagSF*(rewgt-rewErr)*LumiSF);

      //________subBlock_for_altVar_________//
      if(VZGMVAScore >= BDT_CUT && VZGMVAScore <= 1.){
	theHistograms->fill("SYS_mll_central", "SYS_mll_central" , mll_binEdges,    Z->mass(), theWeight*rewgt*PhEffSF*LumiSF);

	theHistograms->fill("SYS_mll_alphas_Up"  , "SYS_mll_alphas_Up"   , mll_binEdges,    Z->mass(), theSampleInfo.alphas_MZ_Up()*theWeight*rewgt*PhEffSF*LumiSF);
	if(isSigSample || isZGSample || isDYSample)      theHistograms->fill("SYS_mll_alphas_Down", "SYS_mll_alphas_Down" , mll_binEdges,    Z->mass(), (2.-theSampleInfo.alphas_MZ_Down())*theWeight*rewgt*PhEffSF*LumiSF);
	else       theHistograms->fill("SYS_mll_alphas_Down", "SYS_mll_alphas_Down" , mll_binEdges,    Z->mass(), theSampleInfo.alphas_MZ_Down()*theWeight*rewgt*PhEffSF*LumiSF);
      
	theHistograms->fill("SYS_mll_PDFVar_Up"  , "SYS_mll_PDFVar_Up"   , mll_binEdges,    Z->mass(), theSampleInfo.PDFVar_Up()*theWeight*rewgt*PhEffSF*LumiSF);
	theHistograms->fill("SYS_mll_PDFVar_Down", "SYS_mll_PDFVar_Down" , mll_binEdges,    Z->mass(), theSampleInfo.PDFVar_Down()*theWeight*rewgt*PhEffSF*LumiSF);

	theHistograms->fill("SYS_mll_QCDscale_Up"  , "SYS_mll_QCDscale_Up"   , mll_binEdges,    Z->mass(), QCDscale_Up*theWeight*rewgt*PhEffSF*LumiSF);
	theHistograms->fill("SYS_mll_QCDscale_Down", "SYS_mll_QCDscale_Down" , mll_binEdges,    Z->mass(), QCDscale_Dn*theWeight*rewgt*PhEffSF*LumiSF);

	theHistograms->fill("SYS_mll_L1Prefiring_Up"  , "SYS_mll_L1Prefiring_Up"   , mll_binEdges,    Z->mass(), (theSampleInfo.L1PrefiringWeightUp()/theSampleInfo.L1PrefiringWeight())*theWeight*rewgt*PhEffSF*LumiSF);
	theHistograms->fill("SYS_mll_L1Prefiring_Down", "SYS_mll_L1Prefiring_Down" , mll_binEdges,    Z->mass(), (theSampleInfo.L1PrefiringWeightDn()/theSampleInfo.L1PrefiringWeight())*theWeight*rewgt*PhEffSF*LumiSF);

	theHistograms->fill("SYS_mll_electronEff_Up"  , "SYS_mll_electronEff_Up"   , mll_binEdges,  Z->mass(), theWeight*rewgt*(1 + eleEff_w)*LumiSF);
	theHistograms->fill("SYS_mll_electronEff_Down", "SYS_mll_electronEff_Down" , mll_binEdges,  Z->mass(), theWeight*rewgt*(1 - eleEff_w)*LumiSF);
	theHistograms->fill("SYS_mll_muonEff_Up"  , "SYS_mll_muonEff_Up"   , mll_binEdges,  Z->mass(), theWeight*rewgt*(1 + muoEff_w)*LumiSF);
	theHistograms->fill("SYS_mll_muonEff_Down", "SYS_mll_muonEff_Down" , mll_binEdges,  Z->mass(), theWeight*rewgt*(1 - muoEff_w)*LumiSF);

      }
      /*
      if(VZGMVAScore_JESup >= BDT_CUT && VZGMVAScore <= 1.) theHistograms->fill("SYS_mll_BDTcut_Up"  , "SYS_mll_BDTcut_Up"   , mll_binEdges,    Z->mass(), theWeight*rewgt*PhEffSF*LumiSF);
      if(VZGMVAScore_JESdn >= BDT_CUT && VZGMVAScore <= 1.) theHistograms->fill("SYS_mll_BDTcut_Down", "SYS_mll_BDTcut_Down" , mll_binEdges,    Z->mass(), theWeight*rewgt*PhEffSF*LumiSF);
      */ //CT momentaneously turned off to handle JES split 

    }
    //__end of sub-block for altVar__//

    for(const auto& syst: jetSysts){
      if(syst.type==SystType::Nominal || syst.direction==0) continue;
      
      phys::Boson<phys::Jet> recoV_JERC;
      bool haveGoodRECODiJetCand_JERC = false;
      bool haveGoodRECOFJCand_null    = false;

      int altVBTopo = 0;
      altVBTopo = Reconstruct(&recoV_JERC,&recoFJ,&haveGoodRECODiJetCand_JERC,&haveGoodRECOFJCand_null,&selectedphotons.at(0), false,  syst);
      

      if (altVBTopo != 1) continue;

      double jercVaried_VZGMVAScore = -2;
      jercVaried_VZGMVAScore=VZGMVAScoreEval(recoV_JERC, recoFJ,  selectedphotons, altVBTopo, syst, false, 0);

      std::string SYShistName = "SYS_BDTScore_" + systName(syst);

      if(jercVaried_VZGMVAScore <= 1. && jercVaried_VZGMVAScore >= binEdges.at(0)) theHistograms->fill(SYShistName  , SYShistName   , binEdges,  jercVaried_VZGMVAScore, theWeight*rewgt*qgTagSF*PhEffSF*LumiSF);

    }
    
  }
  //_____________________________________________END_OF_BLOCK_FOR_SYS_HISTOS_____________________________________//
  /*
  //_________________________________________________BLOCK_FOR_SYS_CRZX_HISTOS________________________________________//
  if(i==1  && isCR && region==REGION_FOR_FIT && isRunForCR  && (UNBLIND && !theSampleInfo.isMC() && histoType=="all_"+REGION_FOR_FIT)  && cut(1, recoV, recoFJ, selectedphotons, VBTopo, mimicVZGMVAScore) && VBTopo==1){// && LGsolved){
    VZGMVAScore      = VZGMVAScoreEval(recoV,       recoFJ,  selectedphotons, VBTopo,       0, 0, false, 0);
    if(VZGMVAScore <= 1. && VZGMVAScore >= binEdges.at(0) ){
      //      std::cout<<"FILLING DATA PLOTS FOR CRDY"<<endl;
      theHistograms->fill("SYS_BDTScore_central", "SYS_BDTScore_central" , binEdges,  VZGMVAScore, theWeight*PhEffSF*LumiSF);
    }
  }
  /*
  if(region=="CR2P_1VL" && isForSys){
  
    std::  cout << "----------------------------------------" << event << endl;
    std::  cout << "SYS_CRDY ready to enter the plotter" << endl;
    std::  cout << "i= " << i << endl;
    std::  cout << region<<"  " << endl;
    if(cut(1, recoV, recoFJ, selectedphotons, VBTopo, mimicVZGMVAScore)) std::  cout << "passing cut(1, recoV, recoFJ, selectedphotons, VBTopo, mimicVZGMVAScore) " << endl;
    else std::  cout << "FAILING cut(1, recoV, recoFJ, selectedphotons, VBTopo, mimicVZGMVAScore) " << endl;
    if(isForSys) std::cout<<"Passing isForSys"<<endl;
  }
  *//*
  if(i==1 && isCR && region==REGION_FOR_FIT && isRunForCR  && isForSys && cut(1, recoV, recoFJ, selectedphotons, VBTopo, mimicVZGMVAScore)){// && LGsolved){
    //std::  cout << "SYS_CRDY entering plotter" << endl;

    int VBTopo_JERup = 0;
    int VBTopo_JERdn = 0;
    int VBTopo_JESup = 0;
    int VBTopo_JESdn = 0;
    
    phys::Boson<phys::Jet> recoV_JERup, recoV_JERdn, recoV_JESup, recoV_JESdn;
    bool haveGoodRECODiJetCand_JERup=false;
    bool haveGoodRECODiJetCand_JERdn=false;
    bool haveGoodRECODiJetCand_JESup=false;
    bool haveGoodRECODiJetCand_JESdn=false;
    bool haveGoodRECOFJCand_null    =false;

    VBTopo_JERup=Reconstruct(&recoV_JERup,&recoFJ,&haveGoodRECODiJetCand_JERup,&haveGoodRECOFJCand_null,&selectedphotons.at(0), false,  1,  1);
    VBTopo_JERdn=Reconstruct(&recoV_JERdn,&recoFJ,&haveGoodRECODiJetCand_JERdn,&haveGoodRECOFJCand_null,&selectedphotons.at(0), false, -1,  1);
    VBTopo_JESup=Reconstruct(&recoV_JESup,&recoFJ,&haveGoodRECODiJetCand_JESup,&haveGoodRECOFJCand_null,&selectedphotons.at(0), false,  1, -1);
    VBTopo_JESdn=Reconstruct(&recoV_JESdn,&recoFJ,&haveGoodRECODiJetCand_JESdn,&haveGoodRECOFJCand_null,&selectedphotons.at(0), false, -1, -1);

    double VZGMVAScore_JERup= -2.;
    double VZGMVAScore_JERdn= -2.;
    double VZGMVAScore_JESup= -2.;
    double VZGMVAScore_JESdn= -2.;
    double VZGMVAScore_DFSup= -2.;
    double VZGMVAScore_DFSdn= -2.;
    
    if(VBTopo==1){
      VZGMVAScore       = VZGMVAScoreEval(recoV,       recoFJ,  selectedphotons, VBTopo,       0, 0, false, 0);
      VZGMVAScore_DFSup = VZGMVAScoreEval(recoV,       recoFJ,  selectedphotons, VBTopo,       0, 0, true, +1);
      VZGMVAScore_DFSdn = VZGMVAScoreEval(recoV,       recoFJ,  selectedphotons, VBTopo,       0, 0, true, -1);
    }
    if(VBTopo_JERup==1) VZGMVAScore_JERup= VZGMVAScoreEval(recoV_JERup, recoFJ,  selectedphotons, VBTopo_JERup, 1, 1, false, 0);
    if(VBTopo_JERdn==1) VZGMVAScore_JERdn= VZGMVAScoreEval(recoV_JERdn, recoFJ,  selectedphotons, VBTopo_JERdn,-1, 1, false, 0);
    if(VBTopo_JESup==1) VZGMVAScore_JESup= VZGMVAScoreEval(recoV_JESup, recoFJ,  selectedphotons, VBTopo_JESup, 1,-1, false, 0);
    if(VBTopo_JESdn==1) VZGMVAScore_JESdn= VZGMVAScoreEval(recoV_JESdn, recoFJ,  selectedphotons, VBTopo_JESdn,-1,-1, false, 0);

    // QCD scale
    // envelope: consider the six variations: {Do, Central, Up} x {Dn, Central, Up} - (central, central) - (Dn, Dn) - (Up, Up) and use the max and min
    float QCDscale_Up(1.), QCDscale_Dn(1.);
    std::vector<float> envelope {
      theSampleInfo.QCDscale_muR0p5F1(),
      theSampleInfo.QCDscale_muR0p5F2(),
      theSampleInfo.QCDscale_muR1F0p5(),
      theSampleInfo.QCDscale_muR1F2(),
      theSampleInfo.QCDscale_muR2F0p5(),
      theSampleInfo.QCDscale_muR2F1()
    };
    QCDscale_Up = *max_element(envelope.begin(), envelope.end());
    QCDscale_Dn = *min_element(envelope.begin(), envelope.end());

    if(VZGMVAScore <= 1. && VZGMVAScore >= binEdges.at(0)){
      theHistograms->fill("SYS_BDTScore_central", "SYS_BDTScore_central" , binEdges,  VZGMVAScore, theWeight*rewgt*qgTagSF*PhEffSF*LumiSF);

      theHistograms->fill("SYS_BDTScore_alphas_Up"  , "SYS_BDTScore_alphas_Up"   , binEdges,  VZGMVAScore, theSampleInfo.alphas_MZ_Up()*theWeight*rewgt*qgTagSF*PhEffSF*LumiSF);
      if(isSigSample || isZGSample || isDYSample)      theHistograms->fill("SYS_BDTScore_alphas_Down", "SYS_BDTScore_alphas_Down" , binEdges,  VZGMVAScore, (2.-theSampleInfo.alphas_MZ_Down())*theWeight*rewgt*qgTagSF*PhEffSF*LumiSF);
      else       theHistograms->fill("SYS_BDTScore_alphas_Down", "SYS_BDTScore_alphas_Down" , binEdges,  VZGMVAScore, theSampleInfo.alphas_MZ_Down()*theWeight*rewgt*qgTagSF*PhEffSF*LumiSF);
      
      theHistograms->fill("SYS_BDTScore_PDFVar_Up"  , "SYS_BDTScore_PDFVar_Up"   , binEdges,  VZGMVAScore, theSampleInfo.PDFVar_Up()*theWeight*rewgt*qgTagSF*PhEffSF*LumiSF);
      theHistograms->fill("SYS_BDTScore_PDFVar_Down", "SYS_BDTScore_PDFVar_Down" , binEdges,  VZGMVAScore, theSampleInfo.PDFVar_Down()*theWeight*rewgt*qgTagSF*PhEffSF*LumiSF);

      theHistograms->fill("SYS_BDTScore_QCDscale_Up"  , "SYS_BDTScore_QCDscale_Up"   , binEdges,  VZGMVAScore, QCDscale_Up*theWeight*rewgt*qgTagSF*PhEffSF*LumiSF);
      theHistograms->fill("SYS_BDTScore_QCDscale_Down", "SYS_BDTScore_QCDscale_Down" , binEdges,  VZGMVAScore, QCDscale_Dn*theWeight*rewgt*qgTagSF*PhEffSF*LumiSF);

      theHistograms->fill("SYS_BDTScore_L1Prefiring_Up"  , "SYS_BDTScore_L1Prefiring_Up"   , binEdges,  VZGMVAScore, (theSampleInfo.L1PrefiringWeightUp()/theSampleInfo.L1PrefiringWeight())*theWeight*rewgt*qgTagSF*PhEffSF*LumiSF);
      theHistograms->fill("SYS_BDTScore_L1Prefiring_Down", "SYS_BDTScore_L1Prefiring_Down" , binEdges,  VZGMVAScore, (theSampleInfo.L1PrefiringWeightDn()/theSampleInfo.L1PrefiringWeight())*theWeight*rewgt*qgTagSF*PhEffSF*LumiSF);
      
      theHistograms->fill("SYS_BDTScore_puWeight_Up"  , "SYS_BDTScore_puWeight_Up"   , binEdges,  VZGMVAScore, (theSampleInfo.puWeightUncUp()/theSampleInfo.puWeight())*theWeight*rewgt*qgTagSF*PhEffSF*LumiSF);
      theHistograms->fill("SYS_BDTScore_puWeight_Down", "SYS_BDTScore_puWeight_Down" , binEdges,  VZGMVAScore, (theSampleInfo.puWeightUncDn()/theSampleInfo.puWeight())*theWeight*rewgt*qgTagSF*PhEffSF*LumiSF);
    
      if(VZGMVAScore_JERup <= 1. && VZGMVAScore_JERup >= binEdges.at(0)) theHistograms->fill("SYS_BDTScore_JER_Up"  , "SYS_BDTScore_JER_Up"   , binEdges,  VZGMVAScore_JERup, theWeight*rewgt*qgTagSF*PhEffSF*LumiSF);
      if(VZGMVAScore_JERdn <= 1. && VZGMVAScore_JERdn >= binEdges.at(0)) theHistograms->fill("SYS_BDTScore_JER_Down", "SYS_BDTScore_JER_Down" , binEdges,  VZGMVAScore_JERdn, theWeight*rewgt*qgTagSF*PhEffSF*LumiSF);
      if(VZGMVAScore_JESup <= 1. && VZGMVAScore_JESup >= binEdges.at(0)) theHistograms->fill("SYS_BDTScore_JES_Up"  , "SYS_BDTScore_JES_Up"   , binEdges,  VZGMVAScore_JESup, theWeight*rewgt*qgTagSF*PhEffSF*LumiSF);
      if(VZGMVAScore_JESdn <= 1. && VZGMVAScore_JESdn >= binEdges.at(0)) theHistograms->fill("SYS_BDTScore_JES_Down", "SYS_BDTScore_JES_Down" , binEdges,  VZGMVAScore_JESdn, theWeight*rewgt*qgTagSF*PhEffSF*LumiSF);

      if(VZGMVAScore_DFSup <= 1. && VZGMVAScore_DFSup >= binEdges.at(0)) theHistograms->fill("SYS_BDTScore_DeepFlavorQGModeling_Down", "SYS_BDTScore_DeepFlavorQGModeling_Down" , binEdges,  VZGMVAScore_DFSdn, theWeight*rewgt*qgTagSF_up*PhEffSF*LumiSF);
      if(VZGMVAScore_DFSdn <= 1. && VZGMVAScore_DFSdn >= binEdges.at(0)) theHistograms->fill("SYS_BDTScore_DeepFlavorQGModeling_Down", "SYS_BDTScore_DeepFlavorQGModeling_Down" , binEdges,  VZGMVAScore_DFSdn, theWeight*rewgt*qgTagSF_dn*PhEffSF*LumiSF);

    }
  }
  //_____________________________________________END_OF_BLOCK_FOR_SYS_CRZX_HISTOS_____________________________________//
  */
  //if (VBTopo ==0) return; //WRONG
    //  std::cout<<"passing cut 1 in printHistos"<<endl;  

  //_____________________________________________________________BLOCK_FOR_SYS_SR2PFJ__________________________________________________//
  //cout << "Starting block for SYS in SR2PFJ" << endl;

  if(runningSR2PFJ){
    if(VBTopo==-1){
      if (i <= cutsToApply && cut(i, recoV, recoFJ, selectedphotons, VBTopo, mimicVZGMVAScore)){
	//cout<<"FJTopo! --> passing cut "<<i<<"-th "<<endl;
	//if(verbose==true) std::cout<<"cut "<<i<<" filling AAA plot "<<histoType<<endl;
	theHistograms->fill("#AAA_cut_flow_" + histoType, "Cut flow", cutsToApply, 0, cutsToApply, i, (theWeight*PhEffSF*LumiSF));
	theHistograms->fill("#AAA_unw_cut_flow_" + histoType, "Unw. events cut flow", cutsToApply, 0, cutsToApply, i, 1.);      
	if (i==2 && !isCR && isForSys && cut(2, recoV, recoFJ, selectedphotons, VBTopo, mimicVZGMVAScore)){
	  cout << "passing cut 2 in SR2PFJ" << endl;
	  double PNScore=(recoFJ.particleNet().WvsQCD + recoFJ.particleNet().ZvsQCD)/2.;
	  // QCD scale
	  // envelope: consider the six variations: {Do, Central, Up} x {Dn, Central, Up} - (central, central) - (Dn, Dn) - (Up, Up) and use the max and min
	  float QCDscale_Up(1.), QCDscale_Dn(1.);
	  std::vector<float> envelope {
	    theSampleInfo.QCDscale_muR0p5F1(),
	    theSampleInfo.QCDscale_muR0p5F2(),
	    theSampleInfo.QCDscale_muR1F0p5(),
	    theSampleInfo.QCDscale_muR1F2(),
	    theSampleInfo.QCDscale_muR2F0p5(),
	    theSampleInfo.QCDscale_muR2F1()
	  };
	  QCDscale_Up = *max_element(envelope.begin(), envelope.end());
	  QCDscale_Dn = *min_element(envelope.begin(), envelope.end());

	  theHistograms->fill("SYS_ParticleNet_central", "SYS_ParticleNet_central" , 10, 0, 1.0,  PNScore, theWeight*PhEffSF*LumiSF);

	  theHistograms->fill("SYS_ParticleNet_alphas_Up"  , "SYS_ParticleNet_alphas_Up"   , 10, 0, 1.0,  PNScore, theSampleInfo.alphas_MZ_Up()*theWeight*PhEffSF*LumiSF);
	  theHistograms->fill("SYS_ParticleNet_alphas_Down", "SYS_ParticleNet_alphas_Down" , 10, 0, 1.0,  PNScore, theSampleInfo.alphas_MZ_Down()*theWeight*PhEffSF*LumiSF);

	  theHistograms->fill("SYS_ParticleNet_PDFVar_Up"  , "SYS_ParticleNet_PDFVar_Up"   , 10, 0, 1.0,  PNScore, theSampleInfo.PDFVar_Up()*theWeight*PhEffSF*LumiSF);
	  theHistograms->fill("SYS_ParticleNet_PDFVar_Down", "SYS_ParticleNet_PDFVar_Down" , 10, 0, 1.0,  PNScore, theSampleInfo.PDFVar_Down()*theWeight*PhEffSF*LumiSF);

	  theHistograms->fill("SYS_ParticleNet_QCDscale_Up"  , "SYS_ParticleNet_QCDscale_Up"   , 10, 0, 1.0,  PNScore, QCDscale_Up*theWeight*PhEffSF*LumiSF);
	  theHistograms->fill("SYS_ParticleNet_QCDscale_Down", "SYS_ParticleNet_QCDscale_Down" , 10, 0, 1.0,  PNScore, QCDscale_Dn*theWeight*PhEffSF*LumiSF);

	  theHistograms->fill("SYS_ParticleNet_L1Prefiring_Up"  , "SYS_ParticleNet_L1Prefiring_Up"   , 10, 0, 1.0,  PNScore, (theSampleInfo.L1PrefiringWeightUp()/theSampleInfo.L1PrefiringWeight())*theWeight*PhEffSF*LumiSF);
	  theHistograms->fill("SYS_ParticleNet_L1Prefiring_Down", "SYS_ParticleNet_L1Prefiring_Down" , 10, 0, 1.0,  PNScore, (theSampleInfo.L1PrefiringWeightDn()/theSampleInfo.L1PrefiringWeight())*theWeight*PhEffSF*LumiSF);

	}
	printHistos(++i, histoType, recoV, recoFJ, selectedphotons,VBTopo, region, isCR); 
      }
    }
    return;  
  }
  //__________________________________________________________END OF BLOCK SR2PFJ__________________________________________________//
  
  std::vector<std::string> cuts = {"0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16"};
  std::vector<std::string> orders = {"0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13"};

  std::string hadTopo="no_jets";
  if(VBTopo==1)
    hadTopo="diJet_Cand";
  else if (VBTopo==-1)
    hadTopo="FJCand";

  TLorentzVector jjPh = recoV.daughter(0).p4()+recoV.daughter(1).p4()+selectedphotons.at(0).p4();
  double mjjPh=jjPh.M();
  double m2jjPh=jjPh.M2();
  
  TLorentzVector llPh = Z->daughter(0).p4()+Z->daughter(1).p4()+selectedphotons.at(0).p4();
  TLorentzVector lljjPh = Z->daughter(0).p4()+Z->daughter(1).p4()+recoV.daughter(0).p4()+recoV.daughter(1).p4()+selectedphotons.at(0).p4();
  TLorentzVector l0Ph = Z->daughter(0).p4()+selectedphotons.at(0).p4();
  TLorentzVector l1Ph = Z->daughter(1).p4()+selectedphotons.at(0).p4();

  TLorentzVector ZV = Z->daughter(0).p4()+Z->daughter(1).p4()+recoV.daughter(0).p4()+recoV.daughter(1).p4();
  
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



  phys::Photon mostEnergeticPhoton;
  //    if (selectedphotons.size() > 0)
  //  {
  std::stable_sort(selectedphotons.begin(), selectedphotons.end(), phys::EComparator());
  mostEnergeticPhoton = selectedphotons[0];
  //  }


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
  std::pair<phys::Photon, phys::Lepton> nearestChLeptToPhoton;

  if(fabs(physmath::deltaR(Z->daughter(0), mostEnergeticPhoton))<fabs(physmath::deltaR(Z->daughter(1), mostEnergeticPhoton)) )
    nearestChLeptToPhoton={mostEnergeticPhoton, Z->daughter(0)};
  else if (fabs(physmath::deltaR(Z->daughter(0), mostEnergeticPhoton))>fabs(physmath::deltaR(Z->daughter(1), mostEnergeticPhoton)) )
    nearestChLeptToPhoton={mostEnergeticPhoton, Z->daughter(1)};

  float deltaR_L0Gamma=fabs(physmath::deltaR(Z->daughter(0), selectedphotons.at(0)));
  float deltaR_L1Gamma=fabs(physmath::deltaR(Z->daughter(1), selectedphotons.at(0)));
  float deltaR_J0Gamma=fabs(physmath::deltaR(recoV.daughter(0), selectedphotons.at(0)));
  float deltaR_J1Gamma=fabs(physmath::deltaR(recoV.daughter(1), selectedphotons.at(0)));
  float dEtaJJ = fabs(recoV.daughter(0).eta() - recoV.daughter(1).eta());
  float meanEtaJ = (recoV.daughter(0).eta() + recoV.daughter(1).eta())/2.;
  float ZepCorr_G = fabs(selectedphotons.at(0).eta() -  meanEtaJ );
  float ptGamma=selectedphotons.at(0).pt();
  float mllG=  llPh.M();

  float HT =  lljjPh.Pt();

  float VZpT =  ZV.Pt();
  float ZGpT =  llPh.Pt();
  float VGpT =  jjPh.Pt();

  float VZGpTscalSum =  Z->pt()+ptGamma+recoV.pt();
  float ZGpTscalSum  =  Z->pt()+ptGamma;
  float VGpTscalSum  =  recoV.pt()+ptGamma;
  float VZpTscalSum  =  Z->pt()+recoV.pt();
  
  /*
  double VZGMVAScore      = VZGMVAScoreEval(recoV, recoFJ,  selectedphotons, VBTopo, noJERCsys);//invert the comment here to REACTIVATE it for running CRs only
  if(i==2)  VZGMVAScore = VZGMVAScoreEval(recoV, recoFJ,  selectedphotons, VBTopo, noJERCsys);/*
  double VZGMVAScore_JERup= VZGMVAScoreEval(recoV, recoFJ,  selectedphotons, VBTopo, 1, 1);
  double VZGMVAScore_JERdn= VZGMVAScoreEval(recoV, recoFJ,  selectedphotons, VBTopo,-1, 1);
  double VZGMVAScore_JESup= VZGMVAScoreEval(recoV, recoFJ,  selectedphotons, VBTopo, 1,-1);
  double VZGMVAScore_JESdn= VZGMVAScoreEval(recoV, recoFJ,  selectedphotons, VBTopo,-1,-1);
  */
  
  /*____BLOCK_CORR_PLOTS__________________________//
  if (i==1 && !isCR && isForSys && cut(1, recoV, recoFJ, selectedphotons, VBTopo, mimicVZGMVAScore)){

    theHistograms->fill("CORR BDT vs ptJ0", "CORR BDT vs ptJ0;   BDT Score ; ptJ0 [GeV] ", 40, -1, 1,  40, 0,  400, VZGMVAScore, recoV.daughter(0).pt(), theWeight*rewgt*PhEffSF*LumiSF);
    theHistograms->fill("CORR BDT vs ptJ1", "CORR BDT vs ptJ1;   BDT Score ; ptJ1 [GeV] ", 40, -1, 1,  40, 0,  400, VZGMVAScore, recoV.daughter(1).pt(), theWeight*rewgt*PhEffSF*LumiSF);
    theHistograms->fill("CORR BDT vs ptJJ", "CORR BDT vs ptJJ;   BDT Score ; ptJJ [GeV] ", 40, -1, 1,  40, 0,  400, VZGMVAScore, recoV.pt(), theWeight*rewgt*PhEffSF*LumiSF);
    theHistograms->fill("CORR BDT vs mll" , "CORR BDT vs mll;   BDT Score ; mll [GeV] "  , 40, -1, 1, 30, 60, 120, VZGMVAScore, mll, theWeight*rewgt*PhEffSF*LumiSF);
    theHistograms->fill("CORR BDT vs mjj" , "CORR BDT vs mjj;   BDT Score ; mjj [GeV] "  , 40, -1, 1, 35, 50, 120, VZGMVAScore, mjj, theWeight*rewgt*PhEffSF*LumiSF);
    theHistograms->fill("CORR BDT vs mllG", "CORR BDT vs mllG; BDT Score ; mll#gamma [GeV] " , 40, -1, 1, 32, 60, 220, VZGMVAScore, mllPh, theWeight*rewgt*PhEffSF*LumiSF);
    theHistograms->fill("CORR BDT vs mjjG", "CORR BDT vs mjjG; BDT Score ; mjj#gamma [GeV] " , 40, -1, 1, 34, 50, 220, VZGMVAScore, mjjPh, theWeight*rewgt*PhEffSF*LumiSF);
    theHistograms->fill("CORR BDT vs dRLG", "CORR BDT vs dRLG; BDT Score ; #Delta R l-#gamma", 40, -1, 1, 20, 0 , 5  , VZGMVAScore, deltaR_L0Gamma<deltaR_L1Gamma ? deltaR_L0Gamma : deltaR_L1Gamma, theWeight*rewgt*PhEffSF*LumiSF);

  }
  //______________________________*/

  if (i <= cutsToApply && cut(i, recoV, recoFJ, selectedphotons, VBTopo, mimicVZGMVAScore))
  {
    if(verbose==true) std::cout<<"cut "<<i<<" filling AAA plot "<<histoType<<endl;
    //std::cout<<"2ND CALL: In printHistos: region: "<<region<<"          after cut"<<i<<"           histo type:  "<<histoType<<endl;
    theHistograms->fill("#AAA_cut_flow_" + histoType, "Cut flow", cutsToApply, 0, cutsToApply, i, (theWeight*rewgt*PhEffSF*LumiSF));
    theHistograms->fill("#AAA_unw_cut_flow_" + histoType, "Unw. events cut flow", cutsToApply, 0, cutsToApply, i, 1.);      

    if(i==2){
      if(selectedphotons.at(0).cutBasedIDLoose())          theHistograms->fill("photonID_" + histoType +"_DJtopo", "photonID", 4, 0, 4, 1, (theWeight*rewgt*PhEffSF*LumiSF));
      if(selectedphotons.at(0).cutBasedIDMedium())     theHistograms->fill("photonID_" + histoType +"_DJtopo", "photonID", 4, 0, 4, 2, (theWeight*rewgt*PhEffSF*LumiSF));
      if(selectedphotons.at(0).cutBasedIDTight())     theHistograms->fill("photonID_" + histoType+"_DJtopo", "photonID", 4, 0, 4, 3, (theWeight*rewgt*PhEffSF*LumiSF));

     
      int nbOfLooseJets = 0;      
      foreach (const phys::Jet &jet, *jets){
	if (KinematicsOK(jet,ptcut,etacut) && jet.passLooseJetID())
	  nbOfLooseJets++;
      }
      theHistograms->fill("#LooseJetsAK4_kinAcc_" + histoType, "#LooseJetsAK4_kinAcc", 8, 0, 8, nbOfLooseJets, theWeight*rewgt*PhEffSF*LumiSF);

      //___________BLOCK_FOR_DY-REWEIGHTING____________//
      if(isCR && region.find("CR2P_1VL")!=std::string::npos){
    	if(isDYSample && histoType == "nonPrompt_CR2P_1VL" && genVBHelper_.ZtoChLep().size()>0){
	  theHistograms->fill("REWGT-DY_DY-GEN" , "Z pt", DYrewgtBinEdges, genVBHelper_.ZtoChLep()[0].pt(), theWeight*LumiSF);
	  theHistograms->fill("REWGT-DY_DY-RECO", "Z pt", DYrewgtBinEdges, Z->pt(), theWeight*LumiSF);
	}
	if(theSampleInfo.isMC() && !isDYSample && (
	                                          ( isZGSample && histoType == "prompt_CR2P_1VL")
	                                           ||
                                       	          ( isSigSample && (histoType == "sign_CR2P_1VL"|| histoType == "bckg_CR2P_1VL") )
	                                           ||
                                       	          ( !isSigSample && !isZGSample && histoType == "all_CR2P_1VL")
						     )
                                           && genVBHelper_.ZtoChLep().size()>0 ){
	  theHistograms->fill("REWGT-DY_else-GEN" , "Z pt", DYrewgtBinEdges, genVBHelper_.ZtoChLep()[0].pt(), theWeight*LumiSF);
	  theHistograms->fill("REWGT-DY_else-RECO", "Z pt", DYrewgtBinEdges, Z->pt(), theWeight*LumiSF);
	}
	if(!theSampleInfo.isMC() && histoType == "all_CR2P_1VL")//  && Z->size()>0)
	  theHistograms->fill("REWGT-DY_numDATA", "Z pt", DYrewgtBinEdges, Z->pt(), theWeight*LumiSF);
      }
      //____________________________________________//

      //___________BLOCK_FOR_ZG-REWEIGHTING____________//
      if(isCR && region.find("CRZOFF")!=std::string::npos){
    	if(isZGSample && histoType == "prompt_CRZOFF" && genVBHelper_.ZtoChLep().size()>0){
	  theHistograms->fill("REWGT-ZG_ZG-GEN" , "Z pt", ZGrewgtBinEdges, genVBHelper_.ZtoChLep()[0].pt(), theWeight*LumiSF);
	  theHistograms->fill("REWGT-ZG_ZG-RECO", "Z pt", ZGrewgtBinEdges, Z->pt(), theWeight*LumiSF);
	}
	if(theSampleInfo.isMC() && !isZGSample && (
						   ( isDYSample && histoType == "nonPrompt_CRZOFF")
	                                           ||
						   ( isSigSample && (histoType == "sign_CRZOFF"|| histoType == "bckg_CRZOFF") )
	                                           ||
						   ( !isSigSample && !isDYSample && histoType == "all_CRZOFF")
						   )
	   && genVBHelper_.ZtoChLep().size()>0 ){
	  theHistograms->fill("REWGT-ZG_else-GEN" , "Z pt", ZGrewgtBinEdges, genVBHelper_.ZtoChLep()[0].pt(), theWeight*LumiSF);
	  theHistograms->fill("REWGT-ZG_else-RECO", "Z pt", ZGrewgtBinEdges, Z->pt(), theWeight*LumiSF);
	}
	if(!theSampleInfo.isMC() && histoType == "all_CRZOFF")//  && Z->size()>0)
	  theHistograms->fill("REWGT-ZG_numDATA", "Z pt", ZGrewgtBinEdges, Z->pt(), theWeight*LumiSF);
      }
      //____________________________________________//

      
    }


    if(VBTopo==1){
      
      //      theHistograms->fill("V_vs_Z_pt_" + histoType + cuts.at(i), "V_vs_Z_pt_" + histoType + cuts.at(i) +";V p_{t} [GeV/c]; #Z p_{t} [GeV/c]", 30, 0, 300, 30, 0, 300, recoV.pt(), Z->pt(), theWeight*rewgt*PhEffSF*LumiSF);
      //      printHistos(1, "sign", recoV, recoFJ, selectedphotons, VBTopo);
      theHistograms->fill("recoVMass_" + histoType + cuts.at(i), "mass of recoV", 40, 40, 120, recoV.mass(), (theWeight*rewgt*PhEffSF*LumiSF));
      if(isCR && region.find("CR2P_1VL")!=std::string::npos) theHistograms->fill("recoVDaughter0Pt_" + histoType + cuts.at(i), "pt of recoVDaughter0", {0,50,100,150,200,250,300,400,500}, recoV.daughter(0).pt(), (theWeight*rewgt*PhEffSF*LumiSF));
      else if(isCR && region.find("CRZOFF_FSRT")!=std::string::npos) theHistograms->fill("recoVDaughter0Pt_" + histoType + cuts.at(i), "pt of recoVDaughter0", {0,50,100,150,250}, recoV.daughter(0).pt(), (theWeight*rewgt*PhEffSF*LumiSF));
      else if(!isCR) theHistograms->fill("recoVDaughter0Pt_" + histoType + cuts.at(i), "pt of recoVDaughter0", {30,60,90,130,170,210,260,320}, recoV.daughter(0).pt(), (theWeight*rewgt*PhEffSF*LumiSF));
      else theHistograms->fill("recoVDaughter0Pt_" + histoType + cuts.at(i), "pt of recoVDaughter0", 50, 0, 500, recoV.daughter(0).pt(), (theWeight*rewgt*PhEffSF*LumiSF));

      theHistograms->fill("recoVDaughter1Pt_" + histoType + cuts.at(i), "pt of recoVDaughter1", {30,50,80,130}, recoV.daughter(1).pt(), (theWeight*rewgt*PhEffSF*LumiSF));
      theHistograms->fill("recoVPt_" + histoType + cuts.at(i), "pt of recoV", {30,60,90,130,170,210,260,320}, recoV.pt(), (theWeight*rewgt*PhEffSF*LumiSF));

      theHistograms->fill("ZepCorr_" + histoType + cuts.at(i), "ZepCorr_", 50, 0, 5,       ZepCorr_G, (theWeight*rewgt*PhEffSF*LumiSF));

      theHistograms->fill("j0_p(bX)"+histoType + cuts.at(i), "j0_p(bX)"+histoType + cuts.at(i)+"; lead. jet DeepFlavour p(bX)", 20, 0, 1, recoV.daughter(0).deepFlavour().probb + recoV.daughter(0).deepFlavour().probbb + recoV.daughter(0).deepFlavour().problepb, theWeight*rewgt*PhEffSF*LumiSF);
      theHistograms->fill("j1_p(bX)"+histoType + cuts.at(i), "j1_p(bX)"+histoType + cuts.at(i)+"; sublead. jet DeepFlavour p(bX)", 20, 0, 1, recoV.daughter(1).deepFlavour().probb + recoV.daughter(1).deepFlavour().probbb + recoV.daughter(1).deepFlavour().problepb, theWeight*rewgt*PhEffSF*LumiSF);

      
      theHistograms->fill("j0_p(g)"+histoType + cuts.at(i), "j0_p(g)"+histoType + cuts.at(i)+"; lead. jet DeepFlavour p(g)", 20, 0, 1, recoV.daughter(0).deepFlavour().probg, theWeight*rewgt*PhEffSF*LumiSF);
      theHistograms->fill("j1_p(g)"+histoType + cuts.at(i), "j1_p(g)"+histoType + cuts.at(i)+"; sublead. jet DeepFlavour p(g)", 20, 0, 1, recoV.daughter(1).deepFlavour().probg, theWeight*rewgt*PhEffSF*LumiSF);
      theHistograms->fill("jj_p(g)Sum"+histoType + cuts.at(i), "jj_p(g)Sum"+histoType + cuts.at(i)+"; DiJet DeepFlavour p(g)", 20, 0, 1,  0.5*(recoV.daughter(0).deepFlavour().probg + recoV.daughter(1).deepFlavour().probg), theWeight*rewgt*PhEffSF*LumiSF);

      theHistograms->fill("j0_p(uds)"+histoType + cuts.at(i), "j0_p(uds)"+histoType + cuts.at(i)+"; lead. jet DeepFlavour p(uds)", {0,0.2,0.4,0.7,1.}, recoV.daughter(0).deepFlavour().probuds, theWeight*rewgt*PhEffSF*LumiSF);
      theHistograms->fill("j1_p(uds)"+histoType + cuts.at(i), "j1_p(uds)"+histoType + cuts.at(i)+"; sublead. jet DeepFlavour p(uds)", 20, 0, 1, recoV.daughter(1).deepFlavour().probuds, theWeight*rewgt*PhEffSF*LumiSF);
      theHistograms->fill("jj_p(uds)Sum"+histoType + cuts.at(i), "jj_p(uds)Sum"+histoType + cuts.at(i)+"; DiJet DeepFlavour p(uds)", 20, 0, 1,  0.5*(recoV.daughter(0).deepFlavour().probuds + recoV.daughter(1).deepFlavour().probuds), theWeight*rewgt*PhEffSF*LumiSF);

      theHistograms->fill("j0_Girth"+histoType + cuts.at(i), "j0_Girth"+histoType + cuts.at(i)+"; lead. jet Girth", 24, 0, 0.24, recoV.daughter(0).girth(), theWeight*rewgt*PhEffSF*LumiSF);
      theHistograms->fill("j1_Girth"+histoType + cuts.at(i), "j1_Girth"+histoType + cuts.at(i)+"; sublead. jet Girth", 24, 0, 0.24, recoV.daughter(1).girth(), theWeight*rewgt*PhEffSF*LumiSF);



      
    }else if(VBTopo==-1){      
      theHistograms->fill("recoFJMass_" + histoType + cuts.at(i), "mass of recoFJ", 40, 0, 200, recoFJ.mass(), (theWeight*rewgt*PhEffSF*LumiSF));
      theHistograms->fill("recoFJPt_" + histoType + cuts.at(i), "pt of recoFJ", 30, 0, 300, recoFJ.pt(), (theWeight*rewgt*PhEffSF*LumiSF));

    }

    theHistograms->fill("recoZDaughter0Pt_" + histoType + cuts.at(i), "pt of recoZDaughter0", 50, 0, 600, Z->daughter(0).pt(), (theWeight*rewgt*PhEffSF*LumiSF));
    theHistograms->fill("recoZDaughter1Pt_" + histoType + cuts.at(i), "pt of recoZDaughter1", 50, 0, 600, Z->daughter(1).pt(), (theWeight*rewgt*PhEffSF*LumiSF));
    

    std::vector<phys::Jet> kinRECOjets;
    std::vector<phys::Jet> kinRECOfatJets;
    foreach (const phys::Jet &jet, *jets)
      {
	if (KinematicsOK(jet,ptcut,etacut))
	  kinRECOjets.push_back(jet);
      }
    foreach (const phys::Jet &FJ, *jetsAK8)
      {
	if (KinematicsOK(FJ,ptcut,etacut))
	  kinRECOfatJets.push_back(FJ);
      }

    //PlotJets(recoV.daughter(0),recoV.daughter(1), "", theWeight*rewgt*PhEffSF*LumiSF, histoType + cuts.at(i));

    theHistograms->fill("recoZMass_" + histoType + cuts.at(i), "mass of recoZ", 30, 60, 120, Z->mass(), (theWeight*rewgt*PhEffSF*LumiSF));

    theHistograms->fill("recoZPt_" + histoType + cuts.at(i), "pt of recoZ", 50, 0, 600, Z->pt(), (theWeight*rewgt*PhEffSF*LumiSF));
    theHistograms->fill("PRE-RWGT_recoZPt_" + histoType + cuts.at(i), "pre-reweighting pt of recoZ", rewgtBinEdges, Z->pt(), (theWeight*PhEffSF*LumiSF));
    theHistograms->fill("POST-RWGT_recoZPt_" + histoType + cuts.at(i), "post-reweighting pt of recoZ", rewgtBinEdges, Z->pt(), (theWeight*rewgt*PhEffSF*LumiSF));
    
    theHistograms->fill("PRE-RWGT_ZGammaPtSum+_" + histoType + cuts.at(i),  "pre-reweighting ZGamma tot. pt ;Z#gamma p_{T} vect. sum [GeV]", 60, 0, 600,( Z->p4()+selectedphotons.at(0).p4() ).Pt(), (theWeight*PhEffSF*LumiSF));
    theHistograms->fill("POST-RWGT_ZGammaPtSum_" + histoType + cuts.at(i), "post-reweighting ZGamma tot. pt ;Z#gamma p_{T} vect. sum [GeV]", 60, 0, 600,( Z->p4()+selectedphotons.at(0).p4() ).Pt(), (theWeight*rewgt*PhEffSF*LumiSF));

    //    if(isCR && region.find("CRZOFF")!=std::string::npos){

    //PLOTS TO TEST Z PT vs Gamma PT    
    theHistograms->fill("PRE-RWGT_ptGamma_" + histoType + cuts.at(i), "pre-reweighting pt of selected photon;#gamma pt [GeV]", 36,20,200, ptGamma, theWeight*PhEffSF*LumiSF);
    theHistograms->fill("POST-RWGT_ptGamma_" + histoType + cuts.at(i), "post-reweighting pt of selected photon;#gamma pt [GeV]", 36,20,200, ptGamma, theWeight*rewgt*PhEffSF*LumiSF);

    if(ptGamma>20 && ptGamma<40) theHistograms->fill("BIN_PtG20-40_PRE-RWGT_recoZPt_" + histoType + cuts.at(i), "pre-reweighting pt of recoZ", rewgtBinEdges, Z->pt(), (theWeight*PhEffSF*LumiSF));
    else if(ptGamma>40 && ptGamma<60) theHistograms->fill("BIN_PtG40-60_PRE-RWGT_recoZPt_" + histoType + cuts.at(i), "pre-reweighting pt of recoZ", rewgtBinEdges, Z->pt(), (theWeight*PhEffSF*LumiSF));
    else if(ptGamma>60 && ptGamma<80) theHistograms->fill("BIN_PtG60-80_PRE-RWGT_recoZPt_" + histoType + cuts.at(i), "pre-reweighting pt of recoZ", rewgtBinEdges, Z->pt(), (theWeight*PhEffSF*LumiSF));
    else if(ptGamma>80 && ptGamma<100) theHistograms->fill("BIN_PtG80-100_PRE-RWGT_recoZPt_" + histoType + cuts.at(i), "pre-reweighting pt of recoZ", rewgtBinEdges, Z->pt(), (theWeight*PhEffSF*LumiSF));
    else if(ptGamma>100 && ptGamma<200) theHistograms->fill("BIN_PtG100-200_PRE-RWGT_recoZPt_" + histoType + cuts.at(i), "pre-reweighting pt of recoZ", rewgtBinEdges, Z->pt(), (theWeight*PhEffSF*LumiSF));
      //}

    //    if(isCR && region.find("CR2P_1VL")!=std::string::npos ){
    if( recoV.daughter(0).deepFlavour().probb + recoV.daughter(0).deepFlavour().probbb + recoV.daughter(0).deepFlavour().problepb < 0.2770
	&& recoV.daughter(1).deepFlavour().probb + recoV.daughter(1).deepFlavour().probbb + recoV.daughter(1).deepFlavour().problepb < 0.2770
	){//PLOTS TO TEST b-veto in CR2P_1F
      theHistograms->fill("PRE-RWGT_recoZPt_" + histoType +"_DJbVetoed"+ cuts.at(i), "pre-reweighting pt of recoZ", rewgtBinEdges, Z->pt(), (theWeight*PhEffSF*LumiSF));
      theHistograms->fill("PRE-RWGT_ptGamma_" + histoType +"_DJbVetoed"+ cuts.at(i), "pre-reweighting pt of selected photon;#gamma pt [GeV]", 36,20,200, ptGamma, theWeight*PhEffSF*LumiSF);
    }else{
      theHistograms->fill("PRE-RWGT_recoZPt_" + histoType +"_DJbContam"+ cuts.at(i), "pre-reweighting pt of recoZ", rewgtBinEdges, Z->pt(), (theWeight*PhEffSF*LumiSF));
      theHistograms->fill("PRE-RWGT_ptGamma_" + histoType +"_DJbContam"+ cuts.at(i), "pre-reweighting pt of selected photon;#gamma pt [GeV]", 36,20,200, ptGamma, theWeight*PhEffSF*LumiSF);
    }
    if( ! std::any_of(jets->begin(), jets->end(), [](const Jet& j){ auto dF = j.deepFlavour(); return dF.probb + dF.probbb + dF.problepb > 0.2770; })  ){//PLOTS TO TEST b-veto in CR2P_1F
      theHistograms->fill("PRE-RWGT_recoZPt_" + histoType +"_bVetoed"+ cuts.at(i), "pre-reweighting pt of recoZ", rewgtBinEdges, Z->pt(), (theWeight*PhEffSF*LumiSF));
      theHistograms->fill("PRE-RWGT_ptGamma_" + histoType +"_bVetoed"+ cuts.at(i), "pre-reweighting pt of selected photon;#gamma pt [GeV]", 36,20,200, ptGamma, theWeight*PhEffSF*LumiSF);
    }else{
      theHistograms->fill("PRE-RWGT_recoZPt_" + histoType +"_bContam"+ cuts.at(i), "pre-reweighting pt of recoZ", rewgtBinEdges, Z->pt(), (theWeight*PhEffSF*LumiSF));
      theHistograms->fill("PRE-RWGT_ptGamma_" + histoType +"_bContam"+ cuts.at(i), "pre-reweighting pt of selected photon;#gamma pt [GeV]", 36,20,200, ptGamma, theWeight*PhEffSF*LumiSF);      
    }

    //PLOTS TO TEST (DJ)bVetoVar(s) PT vs Gamma PT    
    theHistograms->fill("PRE-RWGT_J0bVetoVar_" + histoType + cuts.at(i), "pre-reweighting j0 b-veto var", 20,0,1, recoV.daughter(0).deepFlavour().probb + recoV.daughter(0).deepFlavour().probbb + recoV.daughter(0).deepFlavour().problepb, theWeight*PhEffSF*LumiSF);
    theHistograms->fill("PRE-RWGT_J1bVetoVar_" + histoType + cuts.at(i), "pre-reweighting j1 b-veto var", 20,0,1, recoV.daughter(1).deepFlavour().probb + recoV.daughter(1).deepFlavour().probbb + recoV.daughter(1).deepFlavour().problepb, theWeight*PhEffSF*LumiSF);
    theHistograms->fill("POST-RWGT_J0bVetoVar_" + histoType + cuts.at(i), "post-reweighting j0 b-veto var", 20,0,1, recoV.daughter(0).deepFlavour().probb + recoV.daughter(0).deepFlavour().probbb + recoV.daughter(0).deepFlavour().problepb, theWeight*rewgt*PhEffSF*LumiSF);
    theHistograms->fill("POST-RWGT_J1bVetoVar_" + histoType + cuts.at(i), "post-reweighting j1 b-veto var", 20,0,1, recoV.daughter(1).deepFlavour().probb + recoV.daughter(1).deepFlavour().probbb + recoV.daughter(1).deepFlavour().problepb, theWeight*rewgt*PhEffSF*LumiSF);

    
      //    }

    
    //! std::any_of(jets->begin(), jets->end(), [](const Jet& j){ auto dF = j.deepFlavour(); return dF.probb + dF.probbb + dF.problepb > 0.2770; })    // Note: this is the medium WP for Legacy samples (102X)

    
    if(fullPlotList)    theHistograms->fill("recoZEta_" + histoType + cuts.at(i), "eta of recoZ", 35, 0, 3.5, fabs(Z->eta()), (theWeight*rewgt*PhEffSF*LumiSF));
    //theHistograms->fill("recoZEnergy_" + histoType + cuts.at(i), "energy of  recoZ", 120, 0, 400, fabs(Z->e()), (theWeight*rewgt*PhEffSF*LumiSF));
    theHistograms->fill("recoZDeltaPhi_" + histoType + cuts.at(i), "dPhi of recoZ", 30, 0, 3.2, fabs(physmath::deltaPhi(Z->daughter(0).phi(), Z->daughter(1).phi())), (theWeight*rewgt*PhEffSF*LumiSF));

    /*
    if(VBTopo==1){
      theHistograms->fill(" HAD TOPO DJ cand mass_" + histoType + cuts.at(i), " HAD TOPO DJ cand mass cand mass_" + histoType + cuts.at(i)+" ; mjj Cand. [GeV]", 28, 50, 120, recoV.mass());
      theHistograms->fill(" HAD TOPO 2J/FJ cand mass_" + histoType + cuts.at(i), " HAD TOPO 2J/FJ cand mass_" + histoType + cuts.at(i)+" ; mVB Cand. [GeV]", 28, 50, 120, recoV.mass());
    }else if (VBTopo==-1){
      theHistograms->fill(" HAD TOPO FJ cand mass_" + histoType + cuts.at(i), " HAD TOPO FJ cand mass cand mass_" + histoType + cuts.at(i)+" ; mFj Cand. [GeV]", 28, 50, 120, recoFJ.mass());
      theHistograms->fill(" HAD TOPO 2J/FJ cand mass_" + histoType + cuts.at(i), " HAD TOPO 2J/FJ cand mass_" + histoType + cuts.at(i)+" ; mVB Cand. [GeV]", 28, 50, 120, recoFJ.mass());
    }

    theHistograms->fill(" HAD TOPO_" + histoType + cuts.at(i), " HAD TOPO_" + histoType + cuts.at(i), 3, -1.5, 1.5, VBTopo-0.5);

    */
      
    if (VBTopo == 1){
      theHistograms->fill("DR_Jets_"+histoType + cuts.at(i), "DR_Jets_"+histoType + cuts.at(i)+"; jets #DeltaR", 50, 0, 5, fabs(physmath::deltaR(recoV.daughter(0),recoV.daughter(1))), theWeight*rewgt*PhEffSF*LumiSF);
      if(fabs(physmath::deltaR(recoV.daughter(0), mostEnergeticPhoton))<fabs(physmath::deltaR(recoV.daughter(1), mostEnergeticPhoton)) )
	nearestRECOjetstoPhoton={mostEnergeticPhoton, recoV.daughter(0)};
      else if (fabs(physmath::deltaR(recoV.daughter(0), mostEnergeticPhoton))>fabs(physmath::deltaR(recoV.daughter(1), mostEnergeticPhoton)) )
	nearestRECOjetstoPhoton={mostEnergeticPhoton, recoV.daughter(1)};


      theHistograms->fill("DR_gammaClosestJet_"+histoType + cuts.at(i), "DR_gammaClosestJet_"+histoType + cuts.at(i)+"; #DeltaR", 50, 0, 5, fabs(physmath::deltaR(nearestRECOjetstoPhoton.first, nearestRECOjetstoPhoton.second)), theWeight*rewgt*PhEffSF*LumiSF);
      theHistograms->fill("dRJ0Gamma_"+histoType + cuts.at(i), "dRJ0Gamma_"+histoType + cuts.at(i)+"; #DeltaR", 36, 0.4, 4, fabs(physmath::deltaR(recoV.daughter(0),mostEnergeticPhoton)), theWeight*rewgt*PhEffSF*LumiSF);
      theHistograms->fill("dRJ1Gamma_"+histoType + cuts.at(i), "dRJ1Gamma_"+histoType + cuts.at(i)+"; #DeltaR", 36, 0.4, 4, fabs(physmath::deltaR(recoV.daughter(1),mostEnergeticPhoton)), theWeight*rewgt*PhEffSF*LumiSF);
      if(fullPlotList) theHistograms->fill("DeltaR_vs_Deltapt_gammaJet"+histoType + cuts.at(i), "DeltaR_vs_Deltapt_gammaJet"+histoType + cuts.at(i)+";#Delta pt [GeV/c] ; #DeltaR", 20, -100, 100, 50, 0, 5, nearestRECOjetstoPhoton.first.pt()-nearestRECOjetstoPhoton.second.pt(),fabs(physmath::deltaR(nearestRECOjetstoPhoton.first, nearestRECOjetstoPhoton.second)), theWeight*rewgt*PhEffSF*LumiSF);

    }
      

    theHistograms->fill("DR_Lept_"+histoType + cuts.at(i), "DR_Lept_"+histoType + cuts.at(i)+"; leptons #DeltaR", 50, 0, 5, fabs(physmath::deltaR(Z->daughter(0),Z->daughter(1))), theWeight*rewgt*PhEffSF*LumiSF);

    theHistograms->fill("DR_gammaClosestLept_"+histoType + cuts.at(i), "DR_gammaClosestLept_"+histoType + cuts.at(i)+"; #DeltaR", 50, 0, 5, fabs(physmath::deltaR(nearestChLeptToPhoton.first, nearestChLeptToPhoton.second)), theWeight*rewgt*PhEffSF*LumiSF);



    if(VBTopo==1)
      {
	theHistograms->fill("mjjG_"+histoType + cuts.at(i), 30, 50, 350, mjjPh, theWeight*rewgt*PhEffSF*LumiSF);
	if(fullPlotList){
	  theHistograms->fill("mjj_vs_mjjG_"+histoType + cuts.at(i), "mjj_vs_mjjG_"+histoType + cuts.at(i)+"; mjj [GeV] ; mjj#gamma [GeV]", 35, 50, 120, 30, 50, 350, mjj, mjjPh, theWeight*rewgt*PhEffSF*LumiSF);
	  theHistograms->fill("mll_vs_mjj_"+histoType + cuts.at(i), "mll_vs_mjj_"+histoType + cuts.at(i)+"; mll [GeV] ; mjj [GeV]", 30, 60, 120, 35, 50, 120, mll, mjj, theWeight*rewgt*PhEffSF*LumiSF);
	  theHistograms->fill("mllG_vs_mjjG_"+histoType + cuts.at(i), "mllG_vs_mjjG_"+histoType + cuts.at(i)+"; mll#gamma [GeV] ; mjj#gamma [GeV]", 80, 50, 450, 30, 50, 350, mllPh, mjjPh, theWeight*rewgt*PhEffSF*LumiSF);
	  theHistograms->fill("mllG_vs_mjj_"+histoType + cuts.at(i), "mllG_vs_mjj_"+histoType + cuts.at(i)+"; mll#gamma [GeV] ; mjj [GeV]", 80, 50, 450, 35, 50, 120, mllPh, mjj, theWeight*rewgt*PhEffSF*LumiSF);
	  theHistograms->fill("mll_vs_mjjG_"+histoType + cuts.at(i), "mll_vs_mjjG_"+histoType + cuts.at(i)+"; mll [GeV] ; mjj#gamma [GeV]", 30, 60, 120, 30, 50, 350, mllPh, mjjPh, theWeight*rewgt*PhEffSF*LumiSF);


	  theHistograms->fill("DRlGs_vs_DRjG_"+histoType + cuts.at(i), "DRlGs_vs_DRjG_"+histoType + cuts.at(i)+"; #DeltaRj#gamma [GeV] ; #DeltaRl#gamma", 50, 0, 5, 50, 0, 5, fabs(physmath::deltaR(nearestRECOjetstoPhoton.first, nearestRECOjetstoPhoton.second)),fabs(physmath::deltaR(nearestChLeptToPhoton.first, nearestChLeptToPhoton.second)), theWeight*rewgt*PhEffSF*LumiSF);
	}
	theHistograms->fill("DALITZ_PLOT_VZG_"+histoType + cuts.at(i), "DALITZ_PLOT_VZG_"+histoType + cuts.at(i)+"; m^2 ll#gamma [GeV] ; m^2  jj#gamma [GeV]", 50, 0, 100000, 50, 0, 100000, m2llPh, m2jjPh, theWeight*rewgt*PhEffSF*LumiSF);

	theHistograms->fill("DALITZ_PLOT_llG_"+histoType + cuts.at(i), "DALITZ_PLOT_llG_"+histoType + cuts.at(i)+"; m^2 l0#gamma [GeV] ; m^2  l1#gamma [GeV]", 50, 0, 50000, 50, 0, 50000, m2l0Ph, m2l1Ph, theWeight*rewgt*PhEffSF*LumiSF);

	theHistograms->fill("relativeLLGmass_vs_cosPhi_ll_"+histoType + cuts.at(i), "relativeLLGmass_vs_cosPhi_ll_"+histoType + cuts.at(i)+"; 2m^2 ll#gamma/(ptl0+ptl1) [GeV]; cos #phi ll", 100, 0, 5000, 10, -1, 1, 2*m2llPh/(Z->daughter(0).pt()+Z->daughter(1).pt()),TMath::Cos(fabs(physmath::deltaPhi(Z->daughter(0).phi(), Z->daughter(1).phi()))), theWeight*rewgt*PhEffSF*LumiSF);

	theHistograms->fill("relativeLLGmass_vs_cosPhi_l0Ph_"+histoType + cuts.at(i), "relativeLLGmass_vs_cosPhi_l0Ph_"+histoType + cuts.at(i)+"; 2m^2 ll#gamma/(ptl0+pt#gamma) [GeV]; cos #phi l0#gamma", 100, 0, 5000, 10, -1, 1, 2*m2llPh/(Z->daughter(0).pt()+selectedphotons.at(0).pt()),TMath::Cos(fabs(physmath::deltaPhi(Z->daughter(0).phi(), selectedphotons.at(0).phi()))), theWeight*rewgt*PhEffSF*LumiSF);
	theHistograms->fill("relativeLLGmass_vs_cosPhi_l1Ph_"+histoType + cuts.at(i), "relativeLLGmass_vs_cosPhi_l1Ph_"+histoType + cuts.at(i)+"; 2m^2 ll#gamma/(ptl1+pt#gamma) [GeV]; cos #phi l1#gamma", 100, 0, 5000, 10, -1, 1, 2*m2llPh/(Z->daughter(1).pt()+selectedphotons.at(0).pt()),TMath::Cos(fabs(physmath::deltaPhi(Z->daughter(1).phi(), selectedphotons.at(0).phi()))), theWeight*rewgt*PhEffSF*LumiSF);

      }    

    if(isCR && region.find("CRZOFF")!=std::string::npos) theHistograms->fill("mllG_"+histoType + cuts.at(i), {60,100,140,220}, mllPh, theWeight*rewgt*PhEffSF*LumiSF);
    else if(isCR && region.find("CR2P_1VL")!=std::string::npos) theHistograms->fill("mllG_"+histoType + cuts.at(i), {80,95,110,125,130,140,150,155,160,165,170,175,180,185,190,195,200,205,210,215,220}, mllPh, theWeight*rewgt*PhEffSF*LumiSF);
    else theHistograms->fill("mllG_"+histoType + cuts.at(i), 16,60,220, mllPh, theWeight*rewgt*PhEffSF*LumiSF);
    
    theHistograms->fill("mll_vs_mllG_"+histoType + cuts.at(i), "mll_vs_mllG_"+histoType + cuts.at(i)+"; mll [GeV] ; mll#gamma [GeV]", 30, 60, 120, 80, 50, 450, mll, mllPh, theWeight*rewgt*PhEffSF*LumiSF);
    theHistograms->fill("mllG_vs_DRlGs_"+histoType + cuts.at(i), "mllG_vs_DRlGs_"+histoType + cuts.at(i)+"; mll#gamma [GeV] ; #DeltaRl#gamma", 60, 150, 450, 50, 0, 5, mllPh, fabs(physmath::deltaR(nearestChLeptToPhoton.first, nearestChLeptToPhoton.second)), theWeight*rewgt*PhEffSF*LumiSF);
    
    theHistograms->fill("MET_"+histoType + cuts.at(i), 40, 0, 200, met->pt(), theWeight*rewgt*PhEffSF*LumiSF);



    

    int kinPhotonsCounter=0;
    int kinNoVLPhotonsCounter=0;
    int VLPhotonsCounter=0;
    int VLNoLoosePhotonsCounter=0;
    int LoosePhotonsCounter=0;
    int MediumPhotonsCounter=0;
    int TightPhotonsCounter=0;


      
    foreach (auto p , *photons)
      {
	if (p.id() == 22 && KinematicsOK(p, 20, 2.4) && !p.hasPixelSeed() && p.passElectronVeto())
	  {
	    kinPhotonsCounter++;

	    if(p.cutBasedID(Photon::IdWp::VeryLoose))
	      {
		VLPhotonsCounter++;
		if (p.cutBasedIDLoose()){
		  LoosePhotonsCounter++;
		  if (p.cutBasedIDMedium()){
		    MediumPhotonsCounter++;
		    if (p.cutBasedIDTight()) TightPhotonsCounter++;
		  }
		}else VLNoLoosePhotonsCounter++;
	      }
	    else kinNoVLPhotonsCounter++;
	  }
      }

    if(i==4){
      theHistograms->fill("#_gamma_Loose_" + histoType, "#_gamma_Loose_", 12, 0, 12, LoosePhotonsCounter, (theWeight*rewgt*PhEffSF*LumiSF));
      theHistograms->fill("#_gamma_Medium_" + histoType, "#_gamma_Medium_", 12, 0, 12, MediumPhotonsCounter, (theWeight*rewgt*PhEffSF*LumiSF));
      theHistograms->fill("#_gamma_Tight_" + histoType, "#_gamma_Tight_", 12, 0, 12, TightPhotonsCounter, (theWeight*rewgt*PhEffSF*LumiSF));
    }
    /*
    if(VBTopo==1){
      for(int l = 0; l<3; l++)
	{
	  theHistograms->fill("FWM_T"+orders.at(l)+"_jets_"+histoType+cuts.at(i), "FWM_T"+orders.at(l)+"_jets_"+histoType+cuts.at(i)+"; H_"+orders.at(l)+"^T jets + gamma", 40, -1, 3,
			      SumFWM(l, 't', jjG), theWeight*rewgt*PhEffSF*LumiSF);

	  theHistograms->fill("FWM_T"+orders.at(l)+"_fullSyst_"+histoType+cuts.at(i), "FWM_T"+orders.at(l)+"_fullSyst_"+histoType+cuts.at(i)+"; H_"+orders.at(l)+"^T jets and gamma", 40, -1, 3,
			      SumFWM(l, 't', lljjG), theWeight*rewgt*PhEffSF*LumiSF);
	}
    }
    */



    
    //    theHistograms->fill("VZGMVAScore", "VZGMVAScore", 31, -2.1, 1.0,  VZGMVAScore, theWeight*rewgt*PhEffSF*LumiSF);

    theHistograms->fill("PhotonMVAIDBinned_"+histoType + cuts.at(i), "PhotonMVAIDBinned_"+histoType + cuts.at(i) +"; Photon MVA ID", {0,0.2,0.4,0.5,0.6,0.7,0.8,0.9,1.},  selectedphotons.at(0).MVAvalue(), theWeight*rewgt*PhEffSF*LumiSF);
    theHistograms->fill("PhotonMVAID_"+histoType + cuts.at(i), "PhotonMVAID_"+histoType + cuts.at(i) +"; Photon MVA ID", 80, -1., 1.,  selectedphotons.at(0).MVAvalue(), theWeight*rewgt*PhEffSF*LumiSF);

    if(verboseControlBlinding) std::cout<<"VZGMVAScore_"<<histoType << cuts.at(i)<<"MVA Score: "<<VZGMVAScore<<endl;
    if(VZGMVAScore>-1. && VZGMVAScore<1.) theHistograms->fill("VZGMVAScore_"+histoType + cuts.at(i), "VZGMVAScore_"+histoType + cuts.at(i) +"; MVA Score", binEdges,  VZGMVAScore, theWeight*rewgt*PhEffSF*LumiSF);
    //theHistograms->fill("VZGMVAScore'shortRange_"+histoType + cuts.at(i), "VZGMVAScore_"+histoType + cuts.at(i) +"; MVA Score", 22, -0.2, 1.0,  VZGMVAScore, theWeight*rewgt*PhEffSF*LumiSF);
    
    if(isCR && region.find("CRFSRT")!=std::string::npos)     theHistograms->fill("ptGamma_"+histoType + cuts.at(i), "ptGamma_"+histoType + cuts.at(i)+";#gamma pt [GeV]", {20,35,50,65,80,100,140}, ptGamma, theWeight*rewgt*PhEffSF*LumiSF);
    else     theHistograms->fill("ptGamma_"+histoType + cuts.at(i), "ptGamma_"+histoType + cuts.at(i)+";#gamma pt [GeV]", 50, 0, 200, ptGamma, theWeight*rewgt*PhEffSF*LumiSF);

    if(isCR && region.find("CRFSRT")!=std::string::npos)     theHistograms->fill("dPhiZG_"+histoType + cuts.at(i), "dPhiZG_"+histoType + cuts.at(i)+";#delta#Phi Z-#gamma", {0,0.4,0.8,1.2,1.6,2.2,3.2}, fabs(physmath::deltaPhi(Z->phi(),selectedphotons.at(0).phi()) ), theWeight*rewgt*PhEffSF*LumiSF);
    else if(isCR && region.find("CRZOFF_DIB")!=std::string::npos)     theHistograms->fill("dPhiZG_"+histoType + cuts.at(i), "dPhiZG_"+histoType + cuts.at(i)+";#delta#Phi Z-#gamma", {0,0.8,1.4,2,2.4,2.8,3.2}, fabs(physmath::deltaPhi(Z->phi(),selectedphotons.at(0).phi()) ), theWeight*rewgt*PhEffSF*LumiSF);
    else     theHistograms->fill("dPhiZG_"+histoType + cuts.at(i), "dPhiZG_"+histoType + cuts.at(i)+";#delta#Phi Z-#gamma", 32, 0, 3.2, fabs(physmath::deltaPhi(Z->phi(),selectedphotons.at(0).phi()) ), theWeight*rewgt*PhEffSF*LumiSF);

    
    theHistograms->fill("mllG_vs_ptGamma_"+histoType + cuts.at(i), "mllG_vs_ptGamma_"+histoType + cuts.at(i)+";mll#gamma [GeV] ; #gamma pt [GeV]", 50, 80, 330, 50, 0, 200, mllG, ptGamma, theWeight*rewgt*PhEffSF*LumiSF);
    theHistograms->fill("mllG_vs_dRL0Gamma_"+histoType + cuts.at(i), "mllG_vs_dRL0Gamma_"+histoType + cuts.at(i)+";mll#gamma [GeV] ; #DeltaR l0 - #gamma", 50, 80, 330, 8, 0, 2.0, mllG, deltaR_L0Gamma, theWeight*rewgt*PhEffSF*LumiSF);
    theHistograms->fill("mllG_vs_dRL1Gamma_"+histoType + cuts.at(i), "mllG_vs_dRL1Gamma_"+histoType + cuts.at(i)+";mll#gamma [GeV] ; #DeltaR l1 - #gamma", 50, 80, 330, 8, 0, 2.0, mllG, deltaR_L1Gamma, theWeight*rewgt*PhEffSF*LumiSF);
    theHistograms->fill("dRL0Gamma_vs_dRL1Gamma_"+histoType + cuts.at(i), "dRL0Gamma_vs_dRL1Gamma_"+histoType + cuts.at(i)+"; #DeltaR l0 - #gamma; #DeltaR l1 - #gamma", 8, 0, 2.0, 8, 0, 2.0, deltaR_L0Gamma, deltaR_L1Gamma, theWeight*rewgt*PhEffSF*LumiSF);
    theHistograms->fill("dRJ0Gamma_vs_dRJ1Gamma_"+histoType + cuts.at(i), "dRJ0Gamma_vs_dRJ1Gamma_"+histoType + cuts.at(i)+"; #DeltaR J0 - #gamma; #DeltaR J1 - #gamma", 8, 0, 2.0, 8, 0, 2.0, deltaR_J0Gamma, deltaR_J1Gamma, theWeight*rewgt*PhEffSF*LumiSF);

    //    theHistograms->fill("H_T"+histoType + cuts.at(i), "H_T "+histoType + cuts.at(i)+"; H_T [GeV]", 60, 0, 600, HT, theWeight*rewgt*PhEffSF*LumiSF);
    theHistograms->fill("pre-Rwgt_VZGpT"+histoType + cuts.at(i), "pre-Rwgt_VZGpT"+ histoType + cuts.at(i)+"; VZ#gamma cand. p_{T} [GeV]", {0,30,60,90,130,170,220,300}, HT,   theWeight*PhEffSF*LumiSF);
    theHistograms->fill("VZGpT"+         histoType + cuts.at(i), "VZGpT"+          histoType + cuts.at(i)+"; VZ#gamma cand. p_{T} [GeV]", {0,30,60,90,130,170,220,300}, HT,   theWeight*rewgt*PhEffSF*LumiSF);

    theHistograms->fill("pre-Rwgt_VZpT"+histoType + cuts.at(i), "pre-Rwgt_VZpT"+   histoType + cuts.at(i)+"; VZ cand. p_{T} [GeV]",       {0,30,60,90,130,170,220,300}, VZpT, theWeight*PhEffSF*LumiSF);
    theHistograms->fill("VZpT"+         histoType + cuts.at(i), "VZpT"+            histoType + cuts.at(i)+"; VZ cand. p_{T} [GeV]",       {0,30,60,90,130,170,220,300}, VZpT, theWeight*rewgt*PhEffSF*LumiSF);

    theHistograms->fill("pre-Rwgt_ZGpT"+histoType + cuts.at(i), "pre-Rwgt_ZGpT"+   histoType + cuts.at(i)+"; Z#gamma cand. p_{T} [GeV]",  {0,30,60,90,130,170,220,300}, ZGpT, theWeight*PhEffSF*LumiSF);
    theHistograms->fill("ZGpT"+         histoType + cuts.at(i), "ZGpT"+            histoType + cuts.at(i)+"; Z#gamma cand. p_{T} [GeV]",  {0,30,60,90,130,170,220,300}, ZGpT, theWeight*rewgt*PhEffSF*LumiSF);

    theHistograms->fill("pre-Rwgt_VGpT"+histoType + cuts.at(i), "pre-Rwgt_VGpT"+   histoType + cuts.at(i)+"; V#gamma cand. p_{T} [GeV]",  {0,30,60,90,130,170,220,300}, VGpT, theWeight*PhEffSF*LumiSF);
    theHistograms->fill("VGpT"+         histoType + cuts.at(i), "VGpT"+            histoType + cuts.at(i)+"; V#gamma cand. p_{T} [GeV]",  {0,30,60,90,130,170,220,300}, VGpT, theWeight*rewgt*PhEffSF*LumiSF);

    theHistograms->fill("pre-Rwgt_VZGpTscalSum"+histoType + cuts.at(i), "pre-Rwgt_VZGpTscalSum"+ histoType + cuts.at(i)+"; p_{T}^{V}+p_{T}^{Z}+p_{T}^{#gamma} [GeV]", {0,30,60,90,130,170,220,300}, VZGpTscalSum, theWeight*PhEffSF*LumiSF);
    theHistograms->fill("VZGpTscalSum"+         histoType + cuts.at(i), "VZGpTscalSum"+          histoType + cuts.at(i)+"; p_{T}^{V}+p_{T}^{Z}+p_{T}^{#gamma} [GeV]", {0,30,60,90,130,170,220,300}, VZGpTscalSum, theWeight*rewgt*PhEffSF*LumiSF);

    theHistograms->fill("pre-Rwgt_VZpTscalSum"+histoType + cuts.at(i), "pre-Rwgt_VZpTscalSum"+   histoType + cuts.at(i)+"; p_{T}^{V}+p_{T}^{Z} [GeV]",       {0,30,60,90,130,170,220,300}, VZpTscalSum, theWeight*PhEffSF*LumiSF);
    theHistograms->fill("VZpTscalSum"+         histoType + cuts.at(i), "VZpTscalSum"+            histoType + cuts.at(i)+"; p_{T}^{V}+p_{T}^{Z} [GeV]",       {0,30,60,90,130,170,220,300}, VZpTscalSum, theWeight*rewgt*PhEffSF*LumiSF);

    theHistograms->fill("pre-Rwgt_ZGpTscalSum"+histoType + cuts.at(i), "pre-Rwgt_ZGpTscalSum"+   histoType + cuts.at(i)+"; p_{T}^{Z}+p_{T}^{#gamma} [GeV]",  {0,30,60,90,130,170,220,300}, ZGpTscalSum, theWeight*PhEffSF*LumiSF);
    theHistograms->fill("ZGpTscalSum"+         histoType + cuts.at(i), "ZGpTscalSum"+            histoType + cuts.at(i)+"; p_{T}^{Z}+p_{T}^{#gamma} [GeV]",  {0,30,60,90,130,170,220,300}, ZGpTscalSum, theWeight*rewgt*PhEffSF*LumiSF);

    theHistograms->fill("pre-Rwgt_VGpTscalSum"+histoType + cuts.at(i), "pre-Rwgt_VGpTscalSum"+   histoType + cuts.at(i)+"; p_{T}^{V}+p_{T}^{#gamma} [GeV]",  {0,30,60,90,130,170,220,300}, VGpTscalSum, theWeight*PhEffSF*LumiSF);
    theHistograms->fill("VGpTscalSum"+         histoType + cuts.at(i), "VGpTscalSum"+            histoType + cuts.at(i)+"; p_{T}^{V}+p_{T}^{#gamma} [GeV]",  {0,30,60,90,130,170,220,300}, VGpTscalSum, theWeight*rewgt*PhEffSF*LumiSF);

    
    //p.cutBasedIDLoose()    
    printHistos(++i, histoType, recoV, recoFJ, selectedphotons,VBTopo, region, isCR); 
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


