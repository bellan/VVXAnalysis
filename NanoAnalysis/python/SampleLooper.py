from __future__ import print_function
import pkgutil
import importlib
from  VVXAnalysis.NanoAnalysis import analyses


import math
import ROOT
ROOT.PyConfig.IgnoreCommandLineOptions = True
from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection
from ZZAnalysis.NanoAnalysis.tools import getLeptons, get_genEventSumw

from VVXAnalysis.NanoAnalysis.EventAnalyzer import EventAnalyzer


pathMC = "/eos/cms/store/group/phys_higgs/cmshzz4l/cjlst/RunIII_byZ1Z2/240820/2022EE/"
pathDATA = "/eos/cms/store/group/phys_higgs/cmshzz4l/cjlst/RunIII_byZ1Z2/240820/2022EE/"


maxEntriesPerSample = None # Use only up to this number of events in each MC sample, for quick tests; use None for no scaling


class SampleLooper:
    
    def __init__(self,dataType='MC'):
        self.samples = [
            dict(name = "VBS",filename = "/home/bellan/Workspace/NanoAOD/ZZTo4l_2Jets_EW.root")
        ]
        self.dataType = dataType
        self.isMC = (self.dataType == 'MC')
        
        self.outFile = ROOT.TFile.Open("VVX_"+ dataType +".root","recreate")

        self.load_analyses()

    def load_analyses(self):           
        for _, module_name, _ in pkgutil.iter_modules(
                analyses.__path__
        ):           
            importlib.import_module(
                f"analyses.{module_name}"
            )


        
        
    def loop(self):

        ## Loop over the samples
        for sample in self.samples:
            sampleName = sample["name"]
            filename   = sample["filename"]
            inputFile = ROOT.TFile.Open(filename)
            
            event = inputFile["Events"]
            genEventSumw = 1.
            if self.isMC:
                genEventSumw = get_genEventSumw(inputFile, maxEntriesPerSample)
          
            nEntries = event.GetEntries()
            iEntry=0
            printEntries=max(5000,nEntries/10)

            ######### Analyse the events in a sample! #############
            eventAnalyzer = EventAnalyzer.registry["VVXAnalyzer"](event, sampleName, self.isMC, genEventSumw) #FIXME
            eventAnalyzer.begin()
            
            while iEntry<nEntries and event.GetEntry(iEntry):
                iEntry+=1
                if iEntry%printEntries == 0 : print("Processing", iEntry)
                eventAnalyzer.init()
                eventAnalyzer.analyze()

            eventAnalyzer.end(self.outFile)
            #######################################################

            
    def end(self):
        self.outFile.Close()


         
