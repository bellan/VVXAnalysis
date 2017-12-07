#if !defined(__CINT__) || defined(__MAKECINT__)

#include <TFile.h>
#include <TTree.h>
#include "ResponseMatrix.h"
#include "DataToUnfold.h"
#include "PurityAndStability.h"
//#include <sys/types.h>
//#include "PersonalInfo.cxx"
#include <sys/stat.h>
#include <boost/foreach.hpp>
#define foreach BOOST_FOREACH
#endif

std::vector<std::string> Variables   = {"Mass","nJets","nIncJets","nJets_Central","Mjj","Mjj_Central","Deta","Deta_Central","PtJet1","PtJet2","EtaJet1","EtaJet2","dRZZ","PtZZ"};

//Build response matrices
void GenerateDistributions(string var, bool madgraph, bool tightregion)
{   
  struct stat st;
  if(stat((var+"_test").c_str(),&st) != 0)  system(("mkdir "+var+"_test").c_str());
  ResponseMatrix     *matrix_jesf  = new ResponseMatrix(0,madgraph,tightregion);
  DataToUnfold       *datatounfold = new DataToUnfold(); 
  PurityAndStability *pas          = new PurityAndStability(madgraph); 
  ResponseMatrix     *matrix       = new ResponseMatrix(0,madgraph,tightregion);

  for(int p=-1; p<2; p++){
    for(int q=-1; q<2; q++){    
      matrix->Build(var,"01","4e",p,q,madgraph);
      matrix->Build(var,"1","4e",p,q,madgraph);
      matrix->Build(var,"0","4e",p,q,madgraph);
      matrix->Build(var,"01","4m",p,q,madgraph);
      matrix->Build(var,"1","4m",p,q,madgraph);
      matrix->Build(var,"0","4m",p,q,madgraph); 
      matrix->Build(var,"01","2e2m",p,q,madgraph);
      matrix->Build(var,"1","2e2m",p,q,madgraph);
      matrix->Build(var,"0","2e2m",p,q,madgraph);
    } 
  }

  if(tightregion == 0){ 
    datatounfold->Build(var,"4e");
    datatounfold->Build(var,"4m");
    datatounfold->Build(var,"2e2m");
    pas->Build(var,"4e");
    pas->Build(var,"4m");
    pas->Build(var,"2e2m");
  }

  matrix_jesf->Build_Syst(var,"01","4e","SFSqDn",madgraph);
  matrix_jesf->Build_Syst(var,"01","4e","SFSqUp",madgraph);
  matrix_jesf->Build_Syst(var,"01","4m","SFSqDn",madgraph);
  matrix_jesf->Build_Syst(var,"01","4m","SFSqUp",madgraph);
  matrix_jesf->Build_Syst(var,"01","2e2m","SFSqDn",madgraph);
  matrix_jesf->Build_Syst(var,"01","2e2m","SFSqUp",madgraph);


  matrix_jesf->Build_Syst(var,"01","4e","EleSFSqDn",madgraph);
  matrix_jesf->Build_Syst(var,"01","4e","EleSFSqUp",madgraph);
  matrix_jesf->Build_Syst(var,"01","4m","EleSFSqDn",madgraph);
  matrix_jesf->Build_Syst(var,"01","4m","EleSFSqUp",madgraph);
  matrix_jesf->Build_Syst(var,"01","2e2m","EleSFSqDn",madgraph);
  matrix_jesf->Build_Syst(var,"01","2e2m","EleSFSqUp",madgraph);

  matrix_jesf->Build_Syst(var,"01","4e","MuSFSqDn",madgraph);
  matrix_jesf->Build_Syst(var,"01","4e","MuSFSqUp",madgraph);
  matrix_jesf->Build_Syst(var,"01","4m","MuSFSqDn",madgraph);
  matrix_jesf->Build_Syst(var,"01","4m","MuSFSqUp",madgraph);
  matrix_jesf->Build_Syst(var,"01","2e2m","MuSFSqDn",madgraph);
  matrix_jesf->Build_Syst(var,"01","2e2m","MuSFSqUp",madgraph);

  matrix_jesf->Build_Syst(var,"01","4e","PuDn",madgraph);
  matrix_jesf->Build_Syst(var,"01","4e","PuUp",madgraph);
  matrix_jesf->Build_Syst(var,"01","4m","PuDn",madgraph);
  matrix_jesf->Build_Syst(var,"01","4m","PuUp",madgraph);
  matrix_jesf->Build_Syst(var,"01","2e2m","PuDn",madgraph);
  matrix_jesf->Build_Syst(var,"01","2e2m","PuUp",madgraph);


  matrix_jesf->Build_Syst(var,"01","4e","PDFDn",madgraph);
  matrix_jesf->Build_Syst(var,"01","4e","PDFUp",madgraph);
  matrix_jesf->Build_Syst(var,"01","4m","PDFDn",madgraph);
  matrix_jesf->Build_Syst(var,"01","4m","PDFUp",madgraph);
  matrix_jesf->Build_Syst(var,"01","2e2m","PDFDn",madgraph);
  matrix_jesf->Build_Syst(var,"01","2e2m","PDFUp",madgraph);


  matrix_jesf->Build_Syst(var,"01","4e","AsDn",madgraph);
  matrix_jesf->Build_Syst(var,"01","4e","AsUp",madgraph);
  matrix_jesf->Build_Syst(var,"01","4m","AsDn",madgraph);
  matrix_jesf->Build_Syst(var,"01","4m","AsUp",madgraph);
  matrix_jesf->Build_Syst(var,"01","2e2m","AsDn",madgraph);
  matrix_jesf->Build_Syst(var,"01","2e2m","AsUp",madgraph);


  if(var!="Mass" &&  var!="dRZZ" && var!="PtZZ"){ 

    matrix_jesf->Build_Syst(var,"01","4e","JESDn",madgraph);
    matrix_jesf->Build_Syst(var,"01","4e","JESUp",madgraph);
    matrix_jesf->Build_Syst(var,"01","4e","JERDn",madgraph);
    matrix_jesf->Build_Syst(var,"01","4e","JERUp",madgraph);
    matrix_jesf->Build_Syst(var,"01","4m","JESDn",madgraph);
    matrix_jesf->Build_Syst(var,"01","4m","JESUp",madgraph);
    matrix_jesf->Build_Syst(var,"01","4m","JERDn",madgraph);
    matrix_jesf->Build_Syst(var,"01","4m","JERUp",madgraph);
    matrix_jesf->Build_Syst(var,"01","2e2m","JESDn",madgraph);
    matrix_jesf->Build_Syst(var,"01","2e2m","JESUp",madgraph);
    matrix_jesf->Build_Syst(var,"01","2e2m","JERDn",madgraph);
    matrix_jesf->Build_Syst(var,"01","2e2m","JERUp",madgraph); 

    // if(tightregion == 0){
    //   datatounfold->Build_JE(var,"4e");
    //   datatounfold->Build_JE(var,"4m");
    //   datatounfold->Build_JE(var,"2e2m");
    // }
  }

  matrix->CloseFiles();
  // matrix->Delete();  
  // pas->Delete();
  // matrix_jesf->Delete();
  // datatounfold->Delete();
  // delete matrix;  
  // delete pas;
  // delete matrix_jesf;
  // delete datatounfold;
 }

void Plot(string var, string date, bool madgraph, bool tightregion){

  std::string SavePage = "~/www/PlotsVV/13TeV/";
  system(("mkdir "+SavePage+date).c_str()); 
  system(("cp "+SavePage+"index.php "+SavePage+date).c_str()); 
  system(("mkdir "+SavePage+date+"/"+ var).c_str()); 
  system(("cp "+SavePage+"index.php "+SavePage+date+"/"+ var).c_str());
  system(("mkdir "+SavePage+date+"/"+ var+"/PurityAndStability").c_str()); 
  system(("cp "+SavePage+"index.php "+SavePage+date+"/"+ var+"/PurityAndStability").c_str());

  ResponseMatrix *matrix = new ResponseMatrix(0,madgraph,tightregion);
  DataToUnfold *datatounfold = new DataToUnfold();
  PurityAndStability *pas = new PurityAndStability(madgraph);
  
  matrix->Plot(var,"4e","01","st",date);
  matrix->Plot(var,"4m","01","st",date);
  matrix->Plot(var,"2e2m","01","st",date);
  
  if(tightregion == 0){
    datatounfold->Plot(var,"4e",date); 
    pas->Plot_PAS(var,"4e",date); 
    datatounfold->Plot(var,"4m",date); 
    pas->Plot_PAS(var,"4m",date); 
    datatounfold->Plot(var,"2e2m",date); 
    pas->Plot_PAS(var,"2e2m",date);
  }
  matrix->CloseFiles();
  matrix->Delete();  
  pas->Delete();
  datatounfold->Delete();
}

//Build weighted response matrices
void GenerateWeightedDistributions(string var, bool madgraph, bool tightregion)
{   
  struct stat st;
  if(stat((var+"_test").c_str(),&st) != 0)  system(("mkdir "+var+"_test").c_str());
  
  ResponseMatrix *w_matrix = new ResponseMatrix(1,madgraph,tightregion);
  w_matrix->Build(var,"01","4e",0,0,madgraph);
  w_matrix->Build(var,"1","4e",0,0,madgraph);
  w_matrix->Build(var,"0","4e",0,0,madgraph);
  w_matrix->Build(var,"01","4m",0,0,madgraph);
  w_matrix->Build(var,"1","4m",0,0,madgraph);
  w_matrix->Build(var,"0","4m",0,0,madgraph); 
  w_matrix->Build(var,"01","2e2m",0,0,madgraph);
  w_matrix->Build(var,"1","2e2m",0,0,madgraph);
  w_matrix->Build(var,"0","2e2m",0,0,madgraph);
  w_matrix->Delete();
}

void AllDistributions_var(string var){
  GenerateDistributions(var,1,0);
  GenerateDistributions(var,0,0);
  GenerateDistributions(var,1,1);
  GenerateDistributions(var,0,1);
  // GenerateWeightedDistributions(var,1,0);
  // GenerateWeightedDistributions(var,0,0);
  // GenerateWeightedDistributions(var,1,1);
  // GenerateWeightedDistributions(var,0,1);
} 

void AllPlots_var(string var, string date){
  Plot(var,date,1,0);
  Plot(var,date,0,0);
  Plot(var,date,0,1);
  Plot(var,date,1,1);
}

void GenerateGenMCUpDownDistributions(string var, bool madgraph, bool tightregion)
{   
  string mad;
  string fr;
  if(madgraph ==1) mad = "_Mad";
  else mad = "_Pow";
  if(tightregion ==1) fr = "_fr";
  else fr = "";

string folderName = "GenMCUpDownDistributions"+ fr+ mad;

struct stat st;
 if(stat((folderName).c_str(),&st) != 0)  system(("mkdir "+folderName).c_str());
  
  
  ResponseMatrix *matrix = new ResponseMatrix(0,madgraph,tightregion);
  matrix->GenMCSystDistributions(var,"01","4e", madgraph); 
  matrix->GenMCSystDistributions(var,"01","4m", madgraph);
  matrix->GenMCSystDistributions(var,"01","2e2m", madgraph);
  matrix->GenMCSystDistributions(var,"1","4e", madgraph);
  matrix->GenMCSystDistributions(var,"1","4m", madgraph);
  matrix->GenMCSystDistributions(var,"1","2e2m", madgraph);
  matrix->GenMCSystDistributions(var,"0","4e", madgraph);
  matrix->GenMCSystDistributions(var,"0","4m", madgraph);
  matrix->GenMCSystDistributions(var,"0","2e2m", madgraph);

  matrix->Delete();
}

void GenerateGenMCUpDown_Var(string var)
{   
  GenerateGenMCUpDownDistributions(var,0,0);
  GenerateGenMCUpDownDistributions(var,0,1);
  GenerateGenMCUpDownDistributions(var,1,0);
  GenerateGenMCUpDownDistributions(var,1,1);
 
}
void GenerateGenMCUpDown_All()
{   
  GenerateGenMCUpDown_Var("nJets");
  GenerateGenMCUpDown_Var("nIncJets");
  GenerateGenMCUpDown_Var("nJets_Central"); 
  GenerateGenMCUpDown_Var("Mjj_Central"); 
  GenerateGenMCUpDown_Var("Deta_Central");
  GenerateGenMCUpDown_Var("Mjj"); 
  GenerateGenMCUpDown_Var("Deta");
  GenerateGenMCUpDown_Var("Mass");
  GenerateGenMCUpDown_Var("PtJet1");
  GenerateGenMCUpDown_Var("PtJet2"); 
  GenerateGenMCUpDown_Var("EtaJet1");
  GenerateGenMCUpDown_Var("EtaJet2");
}

void GenerateGenMGatNLOUpDownDistributions(string var, bool tightregion)
{   
  string fr;
  if(tightregion ==1) fr = "_fr";
  else fr = "";

  string folderName = "GenMCUpDownDistributions"+ fr+ "_MGatNLO";

  struct stat st;
  if(stat((folderName).c_str(),&st) != 0)  system(("mkdir "+folderName).c_str());  

  ResponseMatrix *matrix = new ResponseMatrix(0,1,tightregion);
  matrix->GenMGatNLOSystDistributions(var,"01","4e"); 
  matrix->GenMGatNLOSystDistributions(var,"01","4m");
  matrix->GenMGatNLOSystDistributions(var,"01","2e2m");
  // matrix->GenMGatNLOSystDistributions(var,"1","4e");
  // matrix->GenMGatNLOSystDistributions(var,"1","4m");
  // matrix->GenMGatNLOSystDistributions(var,"1","2e2m");
  // matrix->GenMGatNLOSystDistributions(var,"0","4e");
  // matrix->GenMGatNLOSystDistributions(var,"0","4m");
  // matrix->GenMGatNLOSystDistributions(var,"0","2e2m");

}
void GenerateGenMGatNLOUpDown_All()
{   
  foreach(const std::string &var, Variables){
  GenerateGenMGatNLOUpDownDistributions(var,0);
  GenerateGenMGatNLOUpDownDistributions(var,1);
  }
}

void GenerateAll(){
  foreach(const std::string &var, Variables){
    cout<<var.c_str()<<endl;
    GenerateDistributions(var,0,0);
    GenerateDistributions(var,0,1);
    GenerateDistributions(var,1,0);
    GenerateDistributions(var,1,1);
  }
}

void plotAll(string date){
  foreach(const std::string &var, Variables){
    cout<<var.c_str()<<endl;
    Plot(var,date,0,0);
    Plot(var,date,0,1);
    Plot(var,date,1,0);
    Plot(var,date,1,1); 
  }
}

