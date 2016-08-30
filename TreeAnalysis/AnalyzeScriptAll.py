#! /usr/bin/env python
import os
import sys
from os import walk
from python.CrossInfo import*
What= sys.argv[1]
An = sys.argv[2]


if An == "fake":
    AnalyzerList = [{"analyzer":"FakeRateAnalyzer","region":"MC"}]

elif An == "MC":
    AnalyzerList = [{"analyzer":"ZZMCAnalyzer","region":"MC"}]

elif An == "ZZj":
    AnalyzerList = [{"analyzer":"ZZjAnalyzer","region":"SR"}]#,{"analyzer":"ZZjAnalyzer","region":"CR"}]

elif An == "Reco":
    AnalyzerList = [{"analyzer":"ZZRecoAnalyzer","region":"SR"},{"analyzer":"ZZRecoAnalyzer","region":"CR"}]

elif An == "Data":
    AnalyzerList = [{"analyzer":"ZZRecoAnalyzer","region":"SR"},{"analyzer":"ZZRecoAnalyzer","region":"CR"},{"analyzer":"ZZjAnalyzer","region":"SR"},{"analyzer":"ZZjAnalyzer","region":"CR"}]



datalist    = ["DoubleEG2016B","DoubleMu2016B","MuonEG2016B","SingleEle2016B","DoubleEG2016C","DoubleMu2016C","MuonEG2016C","SingleEle2016C","DoubleEG2016D","DoubleMu2016D","MuonEG2016D","SingleEle2016D","DoubleEG2016E","DoubleMu2016E","MuonEG2016E","SingleEle2016E","DoubleEG2016F","DoubleMu2016F","MuonEG2016F","SingleEle2016F"]

MClist_bkg  = ["DYJetsToLL_M50","TTTo2L2Nu","WZ","WWZ","TTZToLL"]
MClist_sig = ["ZZTo4l","ggZZ2e2mu","ggZZ4mu","ggZZ4e","ZZTo4lamcatnlo","ZZTo2e2muJJ","ZZTo4muJJ","ZZTo4eJJ","ggH126","ZZTo4lamcatnlo"]


MClist = MClist_sig+MClist_bkg

Data_Directory = "samples"
MC_Directory   = "samples"

Csv = "../Producers/python/samples_13TeV_2016.csv"
bashCommand="./python/run.py"

for analyzer in AnalyzerList:
    if analyzer["region"] is "CR":    
        
        AddbashCommand1 = "hadd -f results/"+analyzer["analyzer"]+"_"+analyzer["region"]+"3P1F/data.root"
        AddbashCommand2 = "hadd -f results/"+analyzer["analyzer"]+"_"+analyzer["region"]+"2P2F/data.root"
        AddbashCommand3 = "hadd -f results/"+analyzer["analyzer"]+"_"+analyzer["region"]+"/data.root"
        
        if What=="All" or What=="Data":
            if os.path.exists("results/"+analyzer["analyzer"]+"_"+analyzer["region"]+"3P1F/data.root"): os.system("rm results/"+analyzer["analyzer"]+"_"+analyzer["region"]+"3P1F/data.root")
            if os.path.exists("results/"+analyzer["analyzer"]+"_"+analyzer["region"]+"2P2F/data.root"): os.system("rm results/"+analyzer["analyzer"]+"_"+analyzer["region"]+"2P2F/data.root")
            if os.path.exists("results/"+analyzer["analyzer"]+"_"+analyzer["region"]+"/data.root"): os.system("rm results/"+analyzer["analyzer"]+"_"+analyzer["region"]+"/data.root")      
        
    elif analyzer["region"] is "CR_HZZ":

        AddbashCommand1 = "hadd -f results/"+analyzer["analyzer"]+"_CR3P1F_HZZ/data.root"
        AddbashCommand2 = "hadd -f results/"+analyzer["analyzer"]+"_CR2P2F_HZZ/data.root"
        AddbashCommand3 = "hadd -f results/"+analyzer["analyzer"]+"_CR_HZZ/data.root"
        if What=="All" or What=="Data":            
            if os.path.exists("results/"+analyzer["analyzer"]+"_CR3P1F_HZZ/data.root"): os.system("rm results/"+analyzer["analyzer"]+"_CR3P1F_HZZ/data.root")
            if os.path.exists("results/"+analyzer["analyzer"]+"_CR2P2F_HZZ/data.root"): os.system("rm results/"+analyzer["analyzer"]+"_CR2P2F_HZZ/data.root")
            if os.path.exists("results/"+analyzer["analyzer"]+"_CR_HZZ/data.root"): os.system("rm results/"+analyzer["analyzer"]+"_CR_HZZ/data.root")

    else:
        print "results/"+analyzer["analyzer"]+"_"+analyzer["region"]+"/data.root"
        AddbashCommand1 = "hadd -f results/"+analyzer["analyzer"]+"_"+analyzer["region"]+"/data.root"
        AddbashCommand2 = ""
        AddbashCommand3 = ""

        #if os.path.exists("results/"+analyzer["analyzer"]+"_"+analyzer["region"]+"/data.root"): os.system("rm results/"+analyzer["analyzer"]+"_"+analyzer["region"]+"/data.root")

    if(What=="Data" or What=="All"):

        for dataset in datalist:
            if analyzer['analyzer']!="FakeRateAnalyzer" and analyzer["region"]=="MC": continue 
            print bashCommand+" "+analyzer["analyzer"]+" "+dataset+" -r "+analyzer["region"]+" -d "+Data_Directory+" -c"+Csv
            os.system(bashCommand+" "+analyzer["analyzer"]+" "+dataset+" -r "+analyzer["region"]+" -d "+Data_Directory+" -c"+Csv)
            
            if analyzer["region"] is "CR":

                if os.path.exists("results/"+analyzer["analyzer"]+"_"+analyzer["region"]+"3P1F/"+dataset+".root"):
                    AddbashCommand1 = AddbashCommand1+" "+"results/"+analyzer["analyzer"]+"_"+analyzer["region"]+"3P1F/"+dataset+".root"
                    print "results/"+analyzer["analyzer"]+"_"+analyzer["region"]+"3P1F/"+dataset+".root"
                  
                if os.path.isfile("results/"+analyzer["analyzer"]+"_"+analyzer["region"]+"2P2F/"+dataset+".root"): 
                    print "results/"+analyzer["analyzer"]+"_"+analyzer["region"]+"2P2F/"+dataset+".root"                   
                    AddbashCommand2 = AddbashCommand2+" "+"results/"+analyzer["analyzer"]+"_"+analyzer["region"]+"2P2F/"+dataset+".root"

                if os.path.isfile("results/"+analyzer["analyzer"]+"_"+analyzer["region"]+"/reducible_background_from_"+dataset+".root"): 
                    print "results/"+analyzer["analyzer"]+"_"+analyzer["region"]+"/reducible_background_from_"+dataset+".root"
                    AddbashCommand3 = AddbashCommand3+" "+"results/"+analyzer["analyzer"]+"_"+analyzer["region"]+"/reducible_background_from_"+dataset+".root"
                    print AddbashCommand3

            elif analyzer["region"] is "CR_HZZ":
     
                if os.path.exists("results/"+analyzer["analyzer"]+"_CR3P1F_HZZ/data.root"): os.system("rm results/"+analyzer["analyzer"]+"_CR3P1F_HZZ/data.root")

                if os.path.exists("results/"+analyzer["analyzer"]+"_CR3P1F_HZZ/"+dataset+".root"): AddbashCommand1 = AddbashCommand1+" "+"results/"+analyzer["analyzer"]+"_CR3P1F/"+dataset+".root"
                
                print "results/"+analyzer["analyzer"]+"_CR3P1F_HZZ/"+dataset+".root"
                  
                if os.path.isfile("results/"+analyzer["analyzer"]+"_CR2P2F_HZZ/"+dataset+".root"): 
                    print "results/"+analyzer["analyzer"]+"_CR2P2F_HZZ/"+dataset+".root"                   
                    AddbashCommand2 = AddbashCommand2+" "+"results/"+analyzer["analyzer"]+"_CR2P2F_HZZ/"+dataset+".root"

                if os.path.isfile("results/"+analyzer["analyzer"]+"_CR/reducible_background_from_"+dataset+".root"): 
                    print "results/"+analyzer["analyzer"]+"_CR/reducible_background_from_"+dataset+".root"
                    AddbashCommand3 = AddbashCommand3+" "+"results/"+analyzer["analyzer"]+"_CR_HZZ/reducible_background_from_"+dataset+".root"
                    print AddbashCommand3

            else:
                if os.path.exists("results/"+analyzer["analyzer"]+"_"+analyzer["region"]+"/"+dataset+".root"):
                    AddbashCommand1 = AddbashCommand1+" "+"results/"+analyzer["analyzer"]+"_"+analyzer["region"]+"/"+dataset+".root"
                    AddbashCommand2 = ""
                    AddbashCommand3 = ""
                
        print AddbashCommand1
        os.system(AddbashCommand1)

        print AddbashCommand2
        os.system(AddbashCommand2)

        print AddbashCommand3
        os.system(AddbashCommand3)


    if(What=="MC" or What=="All"):

        for MCset in MClist:
            
            print bashCommand+" "+analyzer["analyzer"]+" "+MCset+" -r "+analyzer["region"]+" -d "+MC_Directory+" -l "+str(Lumi)+" -c"+Csv                
            os.system(bashCommand+" "+analyzer["analyzer"]+" "+MCset+" -r "+analyzer["region"]+" -d "+MC_Directory+" -l "+str(Lumi)+" -c"+Csv)



