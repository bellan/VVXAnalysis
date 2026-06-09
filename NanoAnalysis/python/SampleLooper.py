from __future__ import print_function
import pkgutil
import importlib
from  VVXAnalysis.NanoAnalysis import analyses

import json
import math
import ROOT
ROOT.PyConfig.IgnoreCommandLineOptions = True
from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection
from ZZAnalysis.NanoAnalysis.tools import getLeptons, get_genEventSumw

from VVXAnalysis.NanoAnalysis.EventAnalyzer import EventAnalyzer
from VVXAnalysis.NanoAnalysis.Sample import Sample
from VVXAnalysis.NanoAnalysis.SampleLoader import SampleLoader


maxEntriesPerSample = None # Use only up to this number of events in each MC sample, for quick tests; use None for no scaling


class SampleLooper:
    
    def __init__(self,cfg, samples):
        self.analyzer = cfg.analyzer
        self.regions  = cfg.regions
        self.samples  = samples
        
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
            #FIXME: add check that file exists
            print(sample.path())
            inputFile = ROOT.TFile.Open(sample.path())

            '''
            Set the Gen Event Sum Weight, needed to properly normalize the sample weight
            for data it is always 1, for MC we need to extract it from the counters.
            We then need to pass this information to the Analyzer
            '''
            genEventSumw = get_genEventSumw(inputFile, maxEntriesPerSample) if sample.isMC() else 1.
           
            events = inputFile["Events"]
            nEntries = events.GetEntries()
            iEntry=0
            printEntries=max(5000,nEntries/10)

            ######### Analyse the events in a sample! #############
            eventAnalyzer = EventAnalyzer.registry[self.analyzer]()#(base_configuration)
            eventAnalyzer.init(events, genEventSumw, sample.isMC())
            eventAnalyzer.begin()
            
            while iEntry<nEntries and events.GetEntry(iEntry):
                iEntry+=1
                if iEntry%printEntries == 0 : print("Processing", iEntry)
                eventAnalyzer.getCollections()
                eventAnalyzer.analyze()

            eventAnalyzer.end(sample)
            #######################################################

            
    def end(self):
        self.outFile.Close()


         
