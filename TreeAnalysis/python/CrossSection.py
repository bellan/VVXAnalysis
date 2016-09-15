
#! /usr/bin/env python

##################################
## G. Pinna - L. Finco (UNITO) - Jun 2015 ###
##################################

import ROOT,copy
from ROOT import gSystem, TCanvas, TH1,  TPad, gStyle, TLegend,TGraphAsymmErrors,Math,TArrayD
import CrossInfo_copy
from CrossInfo_copy import* 
import collections
from optparse import OptionParser
import sys

import math
import operator
from Colours import *

##################################################################################################################

########################################### Add inclusive systematic #############################################

##################################################################################################################

def whichGlobSystList(MCSet,optInt):

     List = [{"name":"Trig","value":0.015},{"name":"Lumi","value":0.026},{"name":"Acc","value":0.0}] 
   
     if "fr" in MCSet:  List[2]["value"]= 0.01
     else:  List[2]["value"]= 0.05

     if optInt == True:  List[1]["value"]= 0. #optInt = True if DoInclInt or DoDiffInt are True (xs from the integral of the distribution or xs as a function of #bins)

     return  copy.deepcopy(List)


def addGlobSist(hCent,hUp,hDown,MCSet,FinState,optInt):
  
    GlobSistList = [{"name":"Trig","value":0.},{"name":"Lumi","value":0.},{"name":"Acc","value":0.}] 
    GlobSistList = whichGlobSystList(MCSet,optInt)
    
    Nbins = hCent.GetNbinsX()
    for i in range(1,Nbins+1):
        valUp = (hUp.GetBinContent(i)-hCent.GetBinContent(i))**2
        valDown = (hCent.GetBinContent(i)-hDown.GetBinContent(i))**2
       
        for sist in GlobSistList:
            #for the tight region the trigger uncertainty is added when the 4l cross-section is computed
            if "fr" in MCSet and FinState == '4l' and sist["name"] == "Trig":
                valUp+= 0.
                valDown+= 0.
            else:
                valUp+=(hCent.GetBinContent(i)*sist["value"])**2
                valDown+=(hCent.GetBinContent(i)*sist["value"])**2

        hUp.SetBinContent(i,hCent.GetBinContent(i)+math.sqrt(valUp))
        hDown.SetBinContent(i,hCent.GetBinContent(i)-math.sqrt(valDown))

   
##################################################################################################################

##################### Create TGraphAsymmetricError from up and down systematic distributions #####################

##################################################################################################################


def getUncGraph(HCent,HUp,HDown,MCSet,FinState,Data,Err,onOne,optInt):#onOne is always False in ComputeCrossSection.py: it could be deleted
    
    NBins = HCent.GetNbinsX()
    
    #Add global systematic uncertainties: acceptance, luminosity and (if not tight region and not 4l) trigger  
    if Data == True:
        addGlobSist(HCent,HUp,HDown,MCSet,FinState,optInt)
  
    HCent_cl = HCent.Clone("hCent")
  
    #Systematic uncertainties on 1 and not on points (for theoretical MC uncertainties on Data/MC ratio)==> not used anymore ==> getUncGraphMC
    for j in range(1,NBins+1):
        HCent_cl.SetBinContent(j,1) 

    if onOne == True: grErr = ROOT.TGraphAsymmErrors(HCent_cl)
    else: grErr = ROOT.TGraphAsymmErrors(HCent)
    NBins = HCent.GetNbinsX()
    
    for i in range(1,NBins+1):
            
        statUnc = 0     
        sistUnc_up = 0 
        sistUnc_down = 0 
        totUnc_up = 0 
        totUnc_down = 0
        if Data == True:
            statUnc_up =  HCent.GetBinErrorUp(i)#g.GetErrorYhigh(i) - HCent.GetBinContent(i)
            statUnc_down = HCent.GetBinErrorLow(i)#HCent.GetBinContent(i) - g.GetErrorYlow(i)
        else:
            statUnc_up = HCent.GetBinError(i)
            statUnc_down = HCent.GetBinError(i)

        sistUnc_up = HUp.GetBinContent(i)-HCent.GetBinContent(i)
        sistUnc_down = HCent.GetBinContent(i)- HDown.GetBinContent(i)
        if onOne == True and (1-sistUnc_down)<0.:  sistUnc_down =1
        if Err == "syst":
            totUnc_up = sistUnc_up
            totUnc_down = sistUnc_down 
        elif Err == "statsyst":
            totUnc_up = math.sqrt(statUnc_up*statUnc_up + sistUnc_up*sistUnc_up)
            totUnc_down = math.sqrt(statUnc_down*statUnc_down + sistUnc_down*sistUnc_down)
        
        grErr.SetPointEYhigh(i-1,totUnc_up) 
        grErr.SetPointEYlow(i-1,totUnc_down)
       
    return copy.deepcopy(grErr)



#Systematic uncertainties on 1 and not on points (for theoretical MC uncertainties on Data/MC ratio)
def getUncGraphMC(HCent,HUp,HDown,MCSet,FinState):
    
    NBins = HCent.GetNbinsX()
   
    HCent_cl = HCent.Clone("hCent")
  
    for j in range(1,NBins+1):
        HCent_cl.SetBinContent(j,1) 

    grErrs = ROOT.TGraphAsymmErrors(HCent_cl)
       
    for i in range(1,NBins+1):
      
        sistUnc_up = 0 
        sistUnc_down = 0 
        sistUnc_up = HUp.GetBinContent(i)/HCent.GetBinContent(i)-1
        sistUnc_down = 1-HDown.GetBinContent(i)/HCent.GetBinContent(i)

        grErrs.SetPointEYhigh(i-1,sistUnc_up) 
        grErrs.SetPointEYlow(i-1,sistUnc_down)

    return copy.deepcopy(grErrs)

##############################################################################################################

def getCrossPlot_MC(MCSet,Type,DoNormalized):

    FInMC  = ROOT.TFile("./FinalResults_"+MCSet+"/MC.root")
      
    print Red("######################### Monte Carlo "+ MCSet+ " #######################\n")
   
    hsum2e2mu = ROOT.TH1F()
    hsum4e    = ROOT.TH1F()
    hsum4mu   = ROOT.TH1F()

    hSum = [{"state":hsum2e2mu,"name":'2e2m'},{"state":hsum4e,"name":'4e'},{"state":hsum4mu,"name":'4m'}]    
    hSum_clone = [{"state":hsum2e2mu,"name":'2e2m'},{"state":hsum4e,"name":'4e'},{"state":hsum4mu,"name":'4m'}]    

    hsum2e2muUp = ROOT.TH1F()
    hsum4eUp    = ROOT.TH1F()
    hsum4muUp   = ROOT.TH1F()

    hSumUp = [{"state":hsum2e2muUp,"name":'2e2m'},{"state":hsum4eUp,"name":'4e'},{"state":hsum4muUp,"name":'4m'}]    
    hSumUp_clone = [{"state":hsum2e2muUp,"name":'2e2m'},{"state":hsum4eUp,"name":'4e'},{"state":hsum4muUp,"name":'4m'}]    

    hsum2e2muDown = ROOT.TH1F()
    hsum4eDown    = ROOT.TH1F()
    hsum4muDown   = ROOT.TH1F()

    hSumDown = [{"state":hsum2e2muDown,"name":'2e2m'},{"state":hsum4eDown,"name":'4e'},{"state":hsum4muDown,"name":'4m'}]    
    hSumDown_clone = [{"state":hsum2e2muDown,"name":'2e2m'},{"state":hsum4eDown,"name":'4e'},{"state":hsum4muDown,"name":'4m'}]    

    #fill hSum, hSumUp and hSumDown collections and normalize 4e, 4m and 2e2m event distributions wrt the cross section
    for h,hup,hdown in zip(hSum,hSumUp,hSumDown):

        h["state"] = copy.deepcopy(FInMC.Get("ZZTo"+h["name"]+"_"+Type+"Gen_01"))
        hup["state"] = copy.deepcopy(FInMC.Get("ZZTo"+h["name"]+"_"+Type+"Gen_01_up"))
        hdown["state"] = copy.deepcopy(FInMC.Get("ZZTo"+h["name"]+"_"+Type+"Gen_01_down"))

        if h==None or hup==None or hdown==None :
            print "ERROR no data for",h["name"]
            break

        NormUp   =   hup["state"].Integral(1,-1)/h["state"].Integral(1,-1)
        NormDown = hdown["state"].Integral(1,-1)/h["state"].Integral(1,-1)
         
        print Blue("central")+(" "*(9-len("central"))),
        setCrossSectionMC(h["state"],h["name"],Type,1.,DoNormalized,MCSet) 
        print Blue("syst up")+(" "*(9-len("syst up"))),
        setCrossSectionMC(hup["state"],hup["name"],Type,NormUp,DoNormalized,MCSet)
        print Blue("syst down")+(" "*(9-len("syst down"))),
        setCrossSectionMC(hdown["state"],hdown["name"],Type,NormDown,DoNormalized,MCSet)

    #fill collection clones to keep number of events distributions (needed for normalization to one and for 4l distribution in the wide region)
    #hSum/hSumUp/hSumDown and hSum/hSumUp/hSumDown_clone are different only if DoNormalized == 0, otherwise setCrossSectionMC has no effect on them 
    for h,hup,hdown in zip(hSum_clone,hSumUp_clone,hSumDown_clone):

        h["state"] = copy.deepcopy(FInMC.Get("ZZTo"+h["name"]+"_"+Type+"Gen_01"))
        hup["state"] = copy.deepcopy(FInMC.Get("ZZTo"+h["name"]+"_"+Type+"Gen_01_up"))
        hdown["state"] = copy.deepcopy(FInMC.Get("ZZTo"+h["name"]+"_"+Type+"Gen_01_down"))

        if h==None or hup==None or hdown==None :
            print "ERROR no data for",h["name"]
            break

    #4l final state histograms
    hTOTCross = ROOT.TH1F()
    hTOTCrossUp = ROOT.TH1F()
    hTOTCrossDown = ROOT.TH1F()

    #Define 4l element as 4e final state (just an example)
    hTOTCross = copy.deepcopy(hSum[1]["state"])
    hTOTElem = {"state":hTOTCross,"name":'4l'}
       
    hTOTCrossUp = copy.deepcopy(hSumUp[1]["state"])
    hTOTElemUp = {"state":hTOTCrossUp,"name":'4l'}
        
    hTOTCrossDown = copy.deepcopy(hSumDown[1]["state"])
    hTOTElemDown = {"state":hTOTCrossDown,"name":'4l'}

    #get 4l distributions according to the fiducial region: 
    #in the tight region 4l distribution is just the sum of the three different cross-sections 4e, 4m, 2e2m
    #in the wide region 4e, 4m, 2e2m, 4l distributions are always the same IF NORMALIZED to 1 or 7.5
         # ==>Not for MGatNLO, so it is better to add distributions of the NUMBER OF EVENTS (clone) of the three different final states and then normalize to 7.5)
    
    if "fr" in MCSet or DoNormalized == 1: 
        #sum of cross-sections for fr+not norm, sum of #events for fr+norm and wide+norm
        (hTOTCross,hTOTCrossUp,hTOTCrossDown)=combineCrossTight_MC(hSum,hSumUp,hSumDown) 
       
    else: 
         # sum of #event distributions for wide+not norm (the only case I need _clone histos)
        (hTOTCross,hTOTCrossUp,hTOTCrossDown)=combineCrossTight_MC(hSum_clone,hSumUp_clone,hSumDown_clone)
        # and now compute the cross-section (normalized to 7.5)
        NormUp_TOT = hTOTCrossUp.Integral(1,-1)/hTOTCross.Integral(1,-1)
        NormDown_TOT = hTOTCrossDown.Integral(1,-1)/hTOTCross.Integral(1,-1)
        print Blue("central")+(" "*(9-len("central"))), 
        setCrossSectionMC(hTOTCross,hTOTElem["name"],Type,1.,DoNormalized,MCSet) 
        print Blue("syst up")+(" "*(9-len("syst up"))),
        setCrossSectionMC(hTOTCrossUp,hTOTElemUp["name"],Type,NormUp_TOT,DoNormalized,MCSet) 
        print Blue("syst down")+(" "*(9-len("syst down"))),
        setCrossSectionMC(hTOTCrossDown,hTOTElemDown["name"],Type,NormDown_TOT,DoNormalized,MCSet) 
    
        
    #fill 4l distribution clones to keep number of events distributions (needed for normalization to one and for the 4l distribution in the wide region)
    (hTOTCross_clone,hTOTCrossUp_clone,hTOTCrossDown_clone)=combineCrossTight_MC(hSum_clone,hSumUp_clone,hSumDown_clone)
   
    #Add 4l element to the clone lists
    for hTot,h4lTot in zip([hSum_clone,hSumUp_clone,hSumDown_clone],[hTOTCross_clone,hTOTCrossUp_clone,hTOTCrossDown_clone]):
        hTOTElem_4l = {"state":h4lTot,"name":'4l'}
        hTot.append(hTOTElem_4l)
    
    #To obtain plots in which the bin content is first divided by the width and then the distribution is normalized to one, as it was for the previous ZZ cross section measurement, uncomment the following lines (distribution integral equal to one without mulply by the bin width):    
    #for k in range(0,4):
        #hSum_clone[k]["state"].Scale(1.,"width")

  
    #Loop on central, syst up, syst down distributions
    for hTot,h4lTot,sistSt in zip([hSum,hSumUp,hSumDown],[hTOTCross,hTOTCrossUp,hTOTCrossDown],["central","syst up","syst down"]):
        hTOTElem = {"state":h4lTot,"name":'4l'}
        if "fr" in MCSet:
            print Blue(sistSt)+(" "*(9-len(sistSt))), 
            if DoNormalized: print "{0} Total Number of Events  {1} {2:.6f} \n".format("4l",(25-len("4l")*1000)*" ", hTOTElem["state"].Integral(1,-1)) 
            else: print "{0} Total Cross-Section  {1} {2:.6f} \n".format("4l",(25-len("4l")*1000)*" ", hTOTElem["state"].Integral(1,-1)) 
        
        #Add 4l element to the lists (hSum, hSumUp and hSumDown)
        hTot.append(hTOTElem)
         
        #Loop on final states
        for i in hTot: 
            if i["state"]==None:
                print i["state"]," has no enetries" 
           
            if DoNormalized: 
                if i["name"] == "2e2m": j = 0
                elif i["name"] == "4e": j = 1
                elif i["name"] == "4m": j = 2
                elif i["name"] == "4l": j = 3
                
               #Normalization to one
               #How it is for us (distribution integral equal to one if each bin is multiplied by its width) 
                i["state"].Scale(1./hSum_clone[j]["state"].Integral(0,-1),"width")
                print i["name"], "bin ", i ,"hSum_clone[j].Integral(0,-1)",hSum_clone[j]["state"].Integral(0,-1)
             
                # How it was for the previous paper (Bin content is first divided by the width and then the distribution is normalized to one) 
                #print "i[state] ciao",  i["state"].Integral(0,-1), "hSum_clone[j]",hSum_clone[j]["state"].Integral(0,-1) 
                #i["state"].Scale(1.,"width")
                ##hSum_clone[j]["state"].Scale(1.,"width")
                #print "i[state]",  i["state"].Integral(0,-1), "hSum_clone[j]",hSum_clone[j]["state"].Integral(0,-1) 
                #i["state"].Scale(1./hSum_clone[j]["state"].Integral(0,-1))#,"width")

              
            else: i["state"].Scale(1.,"width")
            #else: i["state"].Scale(1.) # Previous ZZ                                 
    return hSum,hSumUp,hSumDown


####################################################################################

def getCrossPlot_Data(MCSet,UseUnfold,Type,Sign,UseMCReco,DoNormalized,optInt):

    if UseMCReco:  print Red("########################### MC RECO ########################\n")
    else: print Red("############################ DATA  #########################\n".format(Sign))
    
    hsum2e2mu = ROOT.TH1F()
    hsum4e    = ROOT.TH1F()
    hsum4mu   = ROOT.TH1F()

    hSum = [{"state":hsum2e2mu,"name":'2e2m'},{"state":hsum4e,"name":'4e'},{"state":hsum4mu,"name":'4m'}]    
    hSum_clone = [{"state":hsum2e2mu,"name":'2e2m'},{"state":hsum4e,"name":'4e'},{"state":hsum4mu,"name":'4m'}]    

    hsum2e2muUp = ROOT.TH1F()
    hsum4eUp    = ROOT.TH1F()
    hsum4muUp   = ROOT.TH1F()

    hSumUp = [{"state":hsum2e2muUp,"name":'2e2m'},{"state":hsum4eUp,"name":'4e'},{"state":hsum4muUp,"name":'4m'}]    
    hSumUp_clone = [{"state":hsum2e2muUp,"name":'2e2m'},{"state":hsum4eUp,"name":'4e'},{"state":hsum4muUp,"name":'4m'}]    

    hsum2e2muDown = ROOT.TH1F()
    hsum4eDown    = ROOT.TH1F()
    hsum4muDown   = ROOT.TH1F()

    hSumDown = [{"state":hsum2e2muDown,"name":'2e2m'},{"state":hsum4eDown,"name":'4e'},{"state":hsum4muDown,"name":'4m'}]    
    hSumDown_clone = [{"state":hsum2e2muDown,"name":'2e2m'},{"state":hsum4eDown,"name":'4e'},{"state":hsum4muDown,"name":'4m'}]    

    if UseMCReco:
        FInCenter = ROOT.TFile("./FinalResults_"+MCSet+"/MCReco.root")
        FInUp     = ROOT.TFile("./FinalResults_"+MCSet+"/DataUp.root") 
        FInDown   = ROOT.TFile("./FinalResults_"+MCSet+"/DataDown.root")  
    else:
        if UseUnfold:
            FInCenter =  ROOT.TFile("./FinalResults_"+MCSet+"/DataUnfold.root")
            FInUp     =  ROOT.TFile( "./FinalResults_"+MCSet+"/DataUnfoldUp.root")
            FInDown   =  ROOT.TFile( "./FinalResults_"+MCSet+"/DataUnfoldDown.root")  
        else: 
            FInCenter = ROOT.TFile("./FinalResults_"+MCSet+"/Data.root") 
            FInUp     = ROOT.TFile("./FinalResults_"+MCSet+"/DataUp.root")
            FInDown   = ROOT.TFile( "./FinalResults_"+MCSet+"/DataDown.root")  
 
    #fill hSum, hSumUp and hSumDown collections and normalize 4e, 4m and 2e2m event distributions by dividing by total lumi (fr, not norm), lumi*BR (wide, not norm) - If DoNormalized == 1, distributions are not modified
    for h,hup,hdown in zip(hSum,hSumUp,hSumDown):

        h["state"] = copy.deepcopy(FInCenter.Get("ZZTo"+h["name"]+"_"+Type+"_01"))
        hup["state"] = copy.deepcopy(FInUp.Get("ZZTo"+h["name"]+"_"+Type+"_01"))
        hdown["state"] = copy.deepcopy(FInDown.Get("ZZTo"+h["name"]+"_"+Type+"_01"))
        #h["state"].Sumw2()
        #hup["state"].Sumw2()
        #hdown["state"].Sumw2()
        h["state"].SetBinErrorOption(ROOT.TH1.kPoisson)
        hup["state"].SetBinErrorOption(ROOT.TH1.kPoisson)
        hdown["state"].SetBinErrorOption(ROOT.TH1.kPoisson)
        
        if h==None:
            print "ERROR no data for",h["name"]
            break

        NormUp   =   hup["state"].Integral(1,-1)/h["state"].Integral(1,-1)
        NormDown = hdown["state"].Integral(1,-1)/h["state"].Integral(1,-1)

        print Blue("Central  "),
        setCrossSectionData(h["state"],h["name"],DoNormalized,MCSet)
        print Blue("Sist Up  "),
        setCrossSectionData(hup["state"],hup["name"],DoNormalized,MCSet)
        print Blue("Sist Down"),
        setCrossSectionData(hdown["state"],hdown["name"],DoNormalized,MCSet)
        
    hTOTCross= ROOT.TH1F() 
    hTOTCrossUp= ROOT.TH1F()
    hTOTCrossDown= ROOT.TH1F()
    
    if "fr" in MCSet: (hTOTCross,hTOTCrossUp,hTOTCrossDown)=combineCrossTight(hSum,hSumUp,hSumDown,MCSet,optInt)   
    else: (hTOTCross,hTOTCrossUp,hTOTCrossDown)=combineCross(hSum,hSumUp,hSumDown) 
   
    #fill distribution clones to keep number of events distributions (needed for normalization to one and for the 4l distribution in the wide region) 
    #At this stage _clone distributions are #events distributions only if DoNormalized == 1
    hSum_clone = copy.deepcopy(hSum)
    hSumUp_clone = copy.deepcopy(hSumUp)
    hSumDown_clone = copy.deepcopy(hSumDown)
    hTOTCross_clone = copy.deepcopy(hTOTCross)
    hTOTCrossUp_clone = copy.deepcopy(hTOTCrossUp)
    hTOTCrossDown_clone = copy.deepcopy(hTOTCrossDown)

    #Add 4l element to the clone lists
    for hTot,h4lTot in zip([hSum_clone,hSumUp_clone,hSumDown_clone],[hTOTCross_clone,hTOTCrossUp_clone,hTOTCrossDown_clone]):
        hTOTElem_4l = {"state":h4lTot,"name":'4l'}
        hTot.append(hTOTElem_4l)


    for hTot,h4lTot,sistSt in zip([hSum,hSumUp,hSumDown],[hTOTCross,hTOTCrossUp,hTOTCrossDown],["central","Sist Up","Sist Down"]):
        hTOTElem = {"state":h4lTot,"name":'4l'}
        print Blue(sistSt)+(" "*(9-len(sistSt))), 
        if DoNormalized: print "{0} Total Number of Events  {1} {2:.6f} \n".format("4l",(25-len("4l")*1000)*" ", hTOTElem["state"].Integral(1,-1)) 
        else: print "{0} Total Cross-Section  {1} {2:.6f} \n".format("4l",(25-len("4l")*1000)*" ", hTOTElem["state"].Integral(1,-1)) 
        hTot.append(hTOTElem)
       
        for i in hTot:
             if i["state"]==None:
                 print i["state"]," has no enetries" 

             if DoNormalized: 
                 if i["name"] == "2e2m": j = 0
                 elif i["name"] == "4e": j = 1
                 elif i["name"] == "4m": j = 2
                 elif i["name"] == "4l": j = 3
                
                 #Normalization to one 
                 err_to_norm = ROOT.Double(0)
                 int_to_norm = 0
                
                 #How it is for us (distribution integral equal to one if each bin is multiplied by its width)
                 #hSum_clone[j]["state"].GetBinContent(1) and i["state"].GetBinContent(1) (before the normalization) have the same integral
                 #i["state"].GetBinError(1)/int_to_norm before the normalization is the same as i["state"].GetBinError(1) after the normalization
                 #Even if they have the same integral value, the normalization is obtain using hSum_clone[j]["state"].Integral(0,-1) and not hSum_clone[j]["state"].IntegralAndError(0,-1,err_to_norm,"width") because the latter sometimes gave some strange results (I think because of "width") 
                 int_to_norm = hSum_clone[j]["state"].IntegralAndError(0,-1,err_to_norm,"width")
                 i["state"].Scale(1./hSum_clone[j]["state"].Integral(0,-1),"width")
                 #i["state"].Scale(1./int_to_norm) 

                 
                 # How it was for the previous paper (Bin content is first divided by the width and then the distribution is normalized to one) 
                 #int_to_norm = hSum_clone[j]["state"].IntegralAndError(0,-1,err_to_norm)
                 #i["state"].Scale(1.,"width")
                 #i["state"].Scale(1./i["state"].Integral(0,-1))#,"width")FIXME
 
                
                 #Add the uncertainty due to the normalization (the error on the normalization value) to the statistical uncertainty
                 Nbins = i["state"].GetNbinsX()  
                 comb_err_stat = 0
                 for k in range(1,Nbins+1):
                     comb_err_stat = 0
                     #linear sum
                     #comb_err_stat = i["state"].GetBinError(k) + err_to_norm*hSum_clone[j]["state"].GetBinContent(k)/(int_to_norm*int_to_norm)
                     #squared sum (i["state"].GetBinError(k) is already divided by the normalization) 
                     comb_err_stat = math.sqrt(i["state"].GetBinError(k)*i["state"].GetBinError(k) + (err_to_norm*hSum_clone[j]["state"].GetBinContent(k)/(int_to_norm*int_to_norm))*(err_to_norm*hSum_clone[j]["state"].GetBinContent(k)/(int_to_norm*int_to_norm)))
                     
                     #print "ERRORS********", i["state"].GetBinError(k),comb_err_stat
                     i["state"].SetBinError(k,comb_err_stat)

              
             else: i["state"].Scale(1.,"width")
             #else: i["state"].Scale(1.) #Previous ZZ

    return hSum,hSumUp,hSumDown




#########################################################
################# Set Cross Section ######################

def setCrossSectionData(h1,FinState,DoNormalized,MCSet):
    
    if FinState=='4e':     BR=BRele*BRele
    elif FinState=='4m':   BR=BRmu*BRmu
    elif FinState=='2e2m': BR=2*BRmu*BRele
       
    #h1.Sumw2()
    if DoNormalized: h1.Scale(1.)#do nothing
    else:  
        if "fr" in MCSet: h1.Scale(1./(Lumi))#normalization to the total integrated lumi
        else: h1.Scale(1./(Lumi*BR))#normalization to lumi and BR
        
    if DoNormalized: print "{0} Total Number of Events {1} {2:.6f} \n".format(FinState,(25-len(FinState)*1000)*" ", (h1.Integral(1,-1))) # Check total cross section without normalization
    else:  print "{0} Total Cross-Section {1} {2:.6f} \n".format(FinState,(25-len(FinState)*1000)*" ", (h1.Integral(1,-1))) 
  


#####################################################
def setCrossSectionMC(h1,FinState,Type,NormSyst,DoNormalized,MCSet):
     
    if DoNormalized: h1.Scale(1.)#do nothing
    else: #no branching ratio because the normalization is fixed and taken from 1507.06257v1 at NNLO (for tight) and from MCFM 6.6 (for wide)
        if "fr" in MCSet: 
            if FinState=='4e':     norm = 0.00516 #pb
            elif FinState=='4m':   norm = 0.00490 #pb
            elif FinState=='2e2m': norm = 0.01015 #pb
            else: norm = 0.02021 #pb
        else: norm = 7.5 #pb
        totnorm = norm*NormSyst

        Integral = h1.Integral(0,-1) 
        if Type == "Mass" or Type == "Jets" or Type == "CentralJets":   
            h1.Scale(totnorm/Integral)#,"width") 
        else:
            print "For now, only the normalization to 1 is available for variables requiring more than 1 jet."
            sys.exit()
    
    if DoNormalized: print "{0} Total Number of Events {1} {2:.6f} \n".format(FinState,(25-len(FinState)*1000)*" ", (h1.Integral(1,-1))) # Check total cross section without normalization
    else:  print "{0} Total Cross-Section {1} {2:.6f} \n".format(FinState,(25-len(FinState)*1000)*" ", (h1.Integral(1,-1))) 
    

##################################################################################################################
##################### Combine different final state cross section distribution for final 4l ######################
##################################################################################################################

#Compute the 4l cross section measurement in the full phase space (60 <m_Z < 120 GeV), as the weighted average of the three final state distributions, with both statistic ad sytematic errors
def combineCross(HList,HListUp,HListDown):

    HCross=ROOT.TH1F()
    HCrossUp=ROOT.TH1F()
    HCrossDown=ROOT.TH1F()

    #just equal to one of them
    HCross=copy.deepcopy(HList[1]["state"])
    HCrossUp=copy.deepcopy(HListUp[1]["state"])
    HCrossDown=copy.deepcopy(HListDown[1]["state"])

    Nbins= HList[1]["state"].GetNbinsX()

    for i in range(1,Nbins+1):

        Hlist  = zip(HList,HListUp,HListDown)
        #Sort List by entries magnitude, from higher to lower to skip 0 entries bins. Sort wrt the bin content of each bin of the central histo 
        SortedHlist = sorted(Hlist,key=lambda value: value[0]["state"].GetBinContent(i),reverse = True)        

        Cross = 0
        ErrSistUp = 0
        ErrSistDown = 0
        ErrStat = 0

        WeightStat = 0.
        WeightTot = 0.
        WeightSistUp = 0.
        WeightSistDown = 0.

        for elem in SortedHlist: #loop over final states 
            Entries = elem[0]["state"].GetBinContent(i) #bin content of each bin of the central histo 
            
          
            if Entries == 0.:  break   # Because of sorting also other final states will be 0, so use break instead of continue
            
            #print "central, up, down element: ",  elem[0]["state"].GetBinContent(i),elem[1]["state"].GetBinContent(i),elem[2]["state"].GetBinContent(i)

            errStatSq = (elem[0]["state"].GetBinError(i))**2
            errSistUpSq = (elem[1]["state"].GetBinContent(i)-Entries)**2
            errSistDownSq = (elem[2]["state"].GetBinContent(i)-Entries)**2
            errSistSq =  ((elem[1]["state"].GetBinContent(i)-elem[2]["state"].GetBinContent(i))/2.)**2 #Use the average of sistematic up and down

            weightStat =  1./errStatSq

            if errSistUpSq==0:  weightSistUp = 0.
            else:  weightSistUp =  1./errSistUpSq

            if errSistDownSq==0:  weightSistDown = 0.
            else:  weightSistDown =  1./errSistDownSq

            weightTot = 1./(errStatSq+errSistSq)
           
            WeightStat += weightStat
            WeightTot += weightTot
            WeightSistUp += weightSistUp
            WeightSistDown += weightSistDown 
            Cross +=  weightTot*Entries

        Cross = Cross/WeightTot

        ErrStat = math.sqrt(1./WeightStat)        
        ErrSistUp =math.sqrt(1./WeightSistUp)
        ErrSistDown =math.sqrt(1./WeightSistDown)

        HCrossUp.SetBinContent(i,Cross+ErrSistUp)
        HCrossDown.SetBinContent(i,Cross-ErrSistDown)
        HCross.SetBinContent(i,Cross)
        HCross.SetBinError(i,ErrStat)      
 
    return HCross,HCrossUp,HCrossDown

#Compute the 4l cross section measurement in the tighter fiducial region, as the sum of the three final state distributions
def combineCrossTight(HList,HListUp,HListDown,MCSet,optInt):
    
    GlobSistList_ = whichGlobSystList(MCSet,optInt)

    HCross=ROOT.TH1F()
    HCrossUp=ROOT.TH1F()
    HCrossDown=ROOT.TH1F()
    
    trig_unc = GlobSistList_[0]["value"]

    #just equal to one of them
    HCross=copy.deepcopy(HList[1]["state"])
    HCrossUp=copy.deepcopy(HListUp[1]["state"])
    HCrossDown=copy.deepcopy(HListDown[1]["state"])

    Nbins= HList[1]["state"].GetNbinsX()

    for i in range(1,Nbins+1):

        Hlist  = zip(HList, HListUp,HListDown)
        #Sort List by entries magnitude, from higher to lower to skip 0 entries bins (sort final states wrt the bin content of the central distribution and combine only those final states with no 0 entries-> the sum does not care about the order of the different final states)
        SortedHlist = sorted(Hlist,key=lambda value: value[0]["state"].GetBinContent(i),reverse = True)        

        Cross = 0
        ErrSistUpSq = 0
        ErrSistDownSq = 0
        ErrStatSq = 0

        for elem in SortedHlist:#loop over the different final states until a bin of the central distribution has 0 entries
            Entries = elem[0]["state"].GetBinContent(i)
            if Entries == 0.:  break   # Because of sorting also other final states will be 0, so use break instead of continue
            
            #print elem[0]["state"].GetBinContent(i),elem[1]["state"].GetBinContent(i),elem[2]["state"].GetBinContent(i)

            errStatSq = (elem[0]["state"].GetBinError(i))**2
            errSistUpSq = (elem[1]["state"].GetBinContent(i)-Entries)**2
            errSistDownSq = (elem[2]["state"].GetBinContent(i)-Entries)**2
            errSistSq =  ((elem[1]["state"].GetBinContent(i)+elem[2]["state"].GetBinContent(i))/2.)**2 #Use the average of sistematic up and down

            Cross +=  Entries 
            ErrStatSq += errStatSq 
            ErrSistUpSq += errSistUpSq
            ErrSistDownSq += errSistDownSq
       
        sigma_4e = Hlist[0][1]["state"].GetBinContent(i)
        sigma_4m = Hlist[0][2]["state"].GetBinContent(i)
        sigma_2e2m = Hlist[0][0]["state"].GetBinContent(i)
        
        sigma_4l_corrunc = math.sqrt(sigma_4e*sigma_4e +sigma_4m*sigma_4m) +sigma_2e2m 
        corrunc_4l = (trig_unc*sigma_4l_corrunc)**2 #trigger uncertainty 
        
        ErrStat = math.sqrt(ErrStatSq)  
        ErrSistUp =math.sqrt(ErrSistUpSq + corrunc_4l) 
        ErrSistDown =math.sqrt(ErrSistDownSq + corrunc_4l)

        #ErrSistUpSq =math.sqrt(ErrSistUp+((0.078+0.029+0.044)*Cross)**2)#uncertainties coming from the unfolding summed in quadrature + lumi 7.8% (2.6*3), Triggher 2,9% (sqrt(2*(trig_unc)^2)+ trig_unc), SF 4.4% (sqrt(2*(0.03)^2)+ 0.03)
        #ErrSistDownSq =math.sqrt(ErrSistDown+((0.078+0.029+0.044)*Cross)**2)

        HCrossUp.SetBinContent(i,Cross+ErrSistUp)
        HCrossDown.SetBinContent(i,Cross-ErrSistDown)
        HCross.SetBinContent(i,Cross)
        HCross.SetBinError(i,ErrStat)      
 
        #print  "cross bin", i, HCross.GetBinContent(i), HCrossUp.GetBinContent(i),  HCrossDown.GetBinContent(i)

    return HCross,HCrossUp,HCrossDown

def combineCrossTight_MC(HList,HListUp,HListDown):#this function sums the three different final state distributions for the central, up and down histograms

    HCross=ROOT.TH1F()
    HCrossUp=ROOT.TH1F()
    HCrossDown=ROOT.TH1F()

    #now equal to 2e2m distributions, then they will be filled with the sum of 2e2m, 4e and 4m
    HCross=copy.deepcopy(HList[1]["state"])
    HCrossUp=copy.deepcopy(HListUp[1]["state"])
    HCrossDown=copy.deepcopy(HListDown[1]["state"])

    #print "HCross, HCrossUp, HCrossDown ", HCross, HCrossUp,  HCrossDown

    Nbins= HList[1]["state"].GetNbinsX()

    for i in range(1,Nbins+1):

        Hlist  = zip(HList, HListUp,HListDown)
        #Sort List by entries magnitude, from higher to lower to skip 0 entries bins (sort final states wrt the bin content of the central distribution and combine only those final states with no 0 entries-> the sum does not care about the order of the different final states)
        SortedHlist = sorted(Hlist,key=lambda value: value[0]["state"].GetBinContent(i),reverse = True)        
        Cross = 0
        CrossUp = 0
        CrossDown = 0
        Entries = 0
        EntriesUp =0
        EntriesDown =0
       
        for elem in SortedHlist:#loop on the different final states until a bin of the central distribution has 0 entries
            Entries = elem[0]["state"].GetBinContent(i)
            EntriesUp = elem[1]["state"].GetBinContent(i)
            EntriesDown = elem[2]["state"].GetBinContent(i)
            if Entries == 0.:  break   # Because of sorting also other final states will be 0 so use break and not continue
           
            Cross +=  Entries #Sum of the bin content of the different final states
            CrossUp +=  EntriesUp 
            CrossDown +=  EntriesDown
       
        # sigma_4e = HList[1]["state"].GetBinContent(i)
        # sigma_4m = HList[2]["state"].GetBinContent(i)
        # sigma_2e2m = HList[0]["state"].GetBinContent(i)
        # sigma_4e_up = HListUp[1]["state"].GetBinContent(i)
        # sigma_4m_up = HListUp[2]["state"].GetBinContent(i)
        # sigma_2e2m_up = HListUp[0]["state"].GetBinContent(i)
        # sigma_4e_down = HListDown[1]["state"].GetBinContent(i)
        # sigma_4m_down = HListDown[2]["state"].GetBinContent(i)
        # sigma_2e2m_down = HListDown[0]["state"].GetBinContent(i) 
        
        HCrossUp.SetBinContent(i,CrossUp)
        HCrossDown.SetBinContent(i,CrossDown)
        HCross.SetBinContent(i,Cross)
        
#        print "central ", sigma_4e, sigma_4m, sigma_2e2m, Cross, HCross.GetBinContent(i)
 #       print "up ", sigma_4e_up, sigma_4m_up, sigma_2e2m_up, CrossUp, HCrossUp.GetBinContent(i)
  #      print "down ", sigma_4e_down, sigma_4m_down, sigma_2e2m_down, CrossDown, HCrossDown.GetBinContent(i)
       
    return HCross,HCrossUp,HCrossDown

def getHisto(Sign,MCSet,Type):
    
    fIn = ROOT.TFile() 
  
    if Sign==0: fIn = ROOT.TFile("./FinalResults_"+MCSet+"/Data.root") 
    elif Sign==1: fIn= ROOT.TFile("./FinalResults_"+MCSet+"/DataUp.root") 
    elif Sign==-1: fIn= ROOT.TFile( "./FinalResults_"+MCSet+"/DataDown.root")

    #use the following distributions if you want an INCLUSIVE measurement from unfolded data distributions!
    # if Sign==0: fIn = ROOT.TFile("./FinalResults_"+MCSet+"/DataUnfold.root") 
    # elif Sign==1: fIn= ROOT.TFile("./FinalResults_"+MCSet+"/DataUnfoldUp.root") 
    # elif Sign==-1: fIn= ROOT.TFile( "./FinalResults_"+MCSet+"/DataUnfoldDown.root")

    hsum2e2mu = ROOT.TH1F()
    hsum4e    = ROOT.TH1F()
    hsum4mu   = ROOT.TH1F()
   
    hSum = [{"state":hsum2e2mu,"name":'2e2m'},{"state":hsum4e,"name":'4e'},{"state":hsum4mu,"name":'4m'}]    
    for h in hSum:
        if (Type == "Mass") or (Type == "Jets"):
            h1 = fIn.Get("ZZTo"+h["name"]+"_"+Type+"_01")
            h["state"] = copy.deepcopy(h1)
        else:
            print Red("\n WRONG TYPE: Mass or Jets only\n") 
       
    return hSum


def TotalCross(MCSet,Type):

    hData = getHisto(0,MCSet,Type)
    hDataUp = getHisto(1,MCSet,Type)
    hDataDown = getHisto(-1,MCSet,Type)
   
    list2e2m = [0,0,0,0,0]
    list4e = [0,0,0,0,0]
    list4m = [0,0,0,0,0]
 
    CrossDic = collections.OrderedDict()
    CrossDic["2e2m"]=list2e2m
    CrossDic["4e"]=list4e
    CrossDic["4m"]=list4m
    

    for i,j,k in zip(hData,hDataUp,hDataDown):

        if "fr" in MCSet: BR = 1
        else:
            if i["name"]=='4e':     BR=BRele*BRele
            elif i["name"]=='4m':   BR=BRmu*BRmu
            elif i["name"]=='2e2m': BR=2*BRmu*BRele
            if i["name"]=='4l': break
        
     
        ErrStat=ROOT.Double(0.)
        #print ErrStat
        CrossDic[i["name"]][0] = i["state"].IntegralAndError(0,-1,ErrStat)/(BR*Lumi) 
        CrossDic[i["name"]][1] = ErrStat/(BR*Lumi)
        CrossDic[i["name"]][2] = j["state"].Integral(0,-1)/(BR*Lumi) - CrossDic[i["name"]][0]     
        CrossDic[i["name"]][3] = CrossDic[i["name"]][0] - k["state"].Integral(0,-1)/(BR*Lumi)
        
        MidSist =  (CrossDic[i["name"]][2]+ CrossDic[i["name"]][3])/2.
        CrossDic[i["name"]][4] = math.sqrt(CrossDic[i["name"]][1]**2+MidSist**2)
    
    if "fr" in MCSet: list4l = combineInclusiveCrossTight(CrossDic)
    else: list4l = combineInclusiveCross(CrossDic)
    dic4l = {"4l":list4l}
    CrossDic.update(dic4l)

    lumi_unc = GlobSistList_[1]["value"]
    acc_unc = GlobSistList_[2]["value"]

    print Red("\n##################### RESULTS FOR INCLUSIVE CROSS SECTION ####################\n") 
    for fin,value in CrossDic.iteritems():
        value[2]=value[2]**2
        value[3]=value[3]**2 
        if "fr" in MCSet and fin == "4l":#only lumi and acceptance are added (trigger unc is added in combineInclusiveCrossTight)

            value[2]= value[2]+(value[0]*lumi_unc)**2+(value[0]*acc_unc)**2 #added also acceptance uncertainty - FIXME
            value[3]= value[3]+(value[0]*lumi_unc)**2+(value[0]*acc_unc)**2
        else:            
            for sist in GlobSistList_: #Add inclusive systematic errors
                value[2]+=(value[0]*sist["value"])**2
                value[3]+=(value[0]*sist["value"])**2 #CHECK
        value[4]= math.sqrt(value[2]+value[3])
        value[2]=math.sqrt(value[2])
        value[3]=math.sqrt(value[3])

        print Blue("{0}".format(fin)),
        if "fr" in MCSet: print  " {0:.2f} +- {1:.2f} (stat) + {2:.2f} (sist) - {3:.2f} (sist) +- {4:.2f} (Total) [fb]\n".format(value[0]*1000,value[1]*1000,value[2]*1000,value[3]*1000,value[4]*1000)
        else: print " {0:.2f} +- {1:.2f} (stat) + {2:.2f} (sist) - {3:.2f} (sist) +- {4:.2f} (Total) [pb] \n".format(value[0],value[1],value[2],value[3],value[4])



def combineInclusiveCross(Dic):
  
    TotStat= math.sqrt(1./(1./(Dic["2e2m"][1]*Dic["2e2m"][1])+ 1./(Dic["4e"][1]*Dic["4e"][1])+ 1./(Dic["4m"][1]*Dic["4m"][1])))

    TotSistUp= math.sqrt(1/(1/(Dic["2e2m"][2]*Dic["2e2m"][2])+ 1/(Dic["4e"][2]*Dic["4e"][2])+ 1/(Dic["4m"][2]*Dic["4m"][2])))

    TotSistDown= math.sqrt(1/(1/(Dic["2e2m"][3]*Dic["2e2m"][3])+ 1/(Dic["4e"][3]*Dic["4e"][3])+ 1/(Dic["4m"][3]*Dic["4m"][3])))
         
    TotErr= math.sqrt(1/(1/(Dic["2e2m"][4]*Dic["2e2m"][4])+ 1/(Dic["4e"][4]*Dic["4e"][4])+ 1/(Dic["4m"][4]*Dic["4m"][4])))

    WhTot= 1/(1/(Dic["2e2m"][4]*Dic["2e2m"][4])+1/(Dic["4e"][4]*Dic["4e"][4])+1/(Dic["4m"][4]*Dic["4m"][4]))

    TotCross = (Dic["2e2m"][0]/(Dic["2e2m"][4]*Dic["2e2m"][4])+Dic["4e"][0]/(Dic["4e"][4]*Dic["4e"][4])+Dic["4m"][0]/(Dic["4m"][4]*Dic["4m"][4]))*WhTot

    return [TotCross,TotStat,TotSistUp,TotSistDown,TotErr]


def combineInclusiveCrossTight(Dic):
    print Red("Combine inclusive cross sections!!!!  ")
    
    GlobSistList1 = GlobSistList_tight
    
    trig_unc = GlobSistList1[0]["value"]
    
    TotStat= math.sqrt((Dic["2e2m"][1]*Dic["2e2m"][1])+(Dic["4e"][1]*Dic["4e"][1])+(Dic["4m"][1]*Dic["4m"][1]))

    TriggerUnc =  (trig_unc*(math.sqrt(Dic["4e"][0]*Dic["4e"][0]+Dic["4m"][0]*+Dic["4m"][0]) + Dic["2e2m"][0]))**2

    TotSistUp= math.sqrt((Dic["2e2m"][2])**2+(Dic["4e"][2])**2+(Dic["4m"][2])**2 +TriggerUnc)

    TotSistDown= math.sqrt((Dic["2e2m"][3])**2+(Dic["4e"][3])**2+(Dic["4m"][3])**2 +TriggerUnc)

    TotSist = (TotSistUp + TotSistDown)/2          

    TotErr= math.sqrt(TotStat*TotStat+TotSist*TotSist)

    TotCross= Dic["2e2m"][0]+Dic["4e"][0]+Dic["4m"][0]

    print "cross sections: ", Dic["2e2m"][0], Dic["4e"][0],Dic["4m"][0]

    return [TotCross,TotStsat,TotSistUp,TotSistDown,TotErr]


   

           
   
