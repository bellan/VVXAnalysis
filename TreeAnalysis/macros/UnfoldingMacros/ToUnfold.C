#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TFile.h>
#include <TTree.h>
#include "ResponseMatrix.h"
//#include "ResponseMatrix_JE.h"
//#include "ResponseMatrix_JESF.h"
#include "DataToUnfold.h"
#include "PurityAndStability.h"
//#include <sys/types.h>
#include <sys/stat.h>
#endif


//Build response matrices
void GenerateDistributions(string var, bool madgraph, bool tightregion)
{   
  struct stat st;
  if(stat((var+"_test").c_str(),&st) != 0)  system(("mkdir "+var+"_test").c_str());
  ResponseMatrix *matrix_jesf = new ResponseMatrix(0,madgraph,tightregion);
  DataToUnfold *datatounfold = new DataToUnfold(); 
  PurityAndStability *pas = new PurityAndStability(madgraph); 
  ResponseMatrix *matrix = new ResponseMatrix(0,madgraph,tightregion);
  
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

  matrix_jesf->Build_SF(var,"01","4e","SFErrSqMinus",madgraph);
  matrix_jesf->Build_SF(var,"01","4e","SFErrSqPlus",madgraph);
  matrix_jesf->Build_SF(var,"01","4m","SFErrSqMinus",madgraph);
  matrix_jesf->Build_SF(var,"01","4m","SFErrSqPlus",madgraph);
  matrix_jesf->Build_SF(var,"01","2e2m","SFErrSqMinus",madgraph);
  matrix_jesf->Build_SF(var,"01","2e2m","SFErrSqPlus",madgraph);
  if(var!="Mass" &&  var!="dRZZ"){ 
    matrix_jesf->Build_JE(var,"01","4e","JESDown",madgraph);
    matrix_jesf->Build_JE(var,"01","4e","JESUp",madgraph);
    matrix_jesf->Build_JE(var,"01","4e","JERDown",madgraph);
    matrix_jesf->Build_JE(var,"01","4e","JERUp",madgraph);
    matrix_jesf->Build_JE(var,"01","4m","JESDown",madgraph);
    matrix_jesf->Build_JE(var,"01","4m","JESUp",madgraph);
    matrix_jesf->Build_JE(var,"01","4m","JERDown",madgraph);
    matrix_jesf->Build_JE(var,"01","4m","JERUp",madgraph);
    matrix_jesf->Build_JE(var,"01","2e2m","JESDown",madgraph);
    matrix_jesf->Build_JE(var,"01","2e2m","JESUp",madgraph);
    matrix_jesf->Build_JE(var,"01","2e2m","JERDown",madgraph);
    matrix_jesf->Build_JE(var,"01","2e2m","JERUp",madgraph); 

    if(tightregion == 0){ 
      datatounfold->Build_JE(var,"4e");
      datatounfold->Build_JE(var,"4m");
      datatounfold->Build_JE(var,"2e2m");
    }
  }
 }

void Plot(string var, string date, bool madgraph, bool tightregion){

  system(("mkdir ~/www/VBS/"+date).c_str()); 
  system(("cp ~/www/VBS/index.php ~/www/VBS/"+date).c_str()); 
  system(("mkdir ~/www/VBS/"+date+"/"+ var).c_str()); 
  system(("cp ~/www/VBS/index.php ~/www/VBS/"+date+"/"+ var).c_str());
  system(("mkdir ~/www/VBS/"+date+"/"+ var+"/PurityAndStability").c_str()); 
  system(("cp ~/www/VBS/index.php ~/www/VBS/"+date+"/"+ var+"/PurityAndStability").c_str());

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
  GenerateGenMCUpDown_Var("Jets");
  GenerateGenMCUpDown_Var("Jets_Central"); 
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
  GenerateGenMGatNLOUpDownDistributions("Jets",0);
  GenerateGenMGatNLOUpDownDistributions("Jets_Central",0); 
  GenerateGenMGatNLOUpDownDistributions("Mjj_Central",0); 
  GenerateGenMGatNLOUpDownDistributions("Deta_Central",0);
  GenerateGenMGatNLOUpDownDistributions("Mjj",0); 
  GenerateGenMGatNLOUpDownDistributions("Deta",0);
  GenerateGenMGatNLOUpDownDistributions("Mass",0);
  GenerateGenMGatNLOUpDownDistributions("PtJet1",0);
  GenerateGenMGatNLOUpDownDistributions("PtJet2",0); 
  GenerateGenMGatNLOUpDownDistributions("EtaJet1",0);
  GenerateGenMGatNLOUpDownDistributions("EtaJet2",0); 
  GenerateGenMGatNLOUpDownDistributions("Jets",1);
  GenerateGenMGatNLOUpDownDistributions("Jets_Central",1); 
  GenerateGenMGatNLOUpDownDistributions("Mjj_Central",1); 
  GenerateGenMGatNLOUpDownDistributions("Deta_Central",1);
  GenerateGenMGatNLOUpDownDistributions("Mjj",1); 
  GenerateGenMGatNLOUpDownDistributions("Deta",1);
  GenerateGenMGatNLOUpDownDistributions("Mass",1);
  GenerateGenMGatNLOUpDownDistributions("PtJet1",1);
  GenerateGenMGatNLOUpDownDistributions("PtJet2",1); 
  GenerateGenMGatNLOUpDownDistributions("EtaJet1",1);
  GenerateGenMGatNLOUpDownDistributions("EtaJet2",1);
}
