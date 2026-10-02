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

import os, sys, traceback
from concurrent.futures import ProcessPoolExecutor, as_completed

maxEntriesPerSample = 100 # Use only up to this number of events in each MC sample, for quick tests; use None for no scaling


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


    def analyzeSample(self, sample):       
        #FIXME: add check that file exists
        print("analyzing", sample.explainYourself())
        print(sample.path())
        inputFile = ROOT.TFile.Open(sample.path())
        
        
        '''
        Set the Gen Event Sum Weight, needed to properly normalize the sample weight
        for data it is always 1, for MC we need to extract it from the counters.
        We then need to pass this information to the Analyzer
        '''
        genEventSumw = get_genEventSumw(inputFile, maxEntriesPerSample) if sample.isMC() else 1.
        
        events = inputFile.Get("Events")
        nEntries = events.GetEntries()
        iEntry=0
        printEntries=max(5000,nEntries/10)
        
        ######### Analyse the events in a sample! #############
        eventAnalyzer = EventAnalyzer.registry[self.analyzer](self.regions)#(base_configuration)
        eventAnalyzer.init(events, genEventSumw, sample.luminosity, sample.isMC())
        eventAnalyzer.begin()
        
        while iEntry<nEntries and events.GetEntry(iEntry):
            iEntry+=1
            if iEntry%printEntries == 0 : print("Processing", iEntry)
            eventAnalyzer.eventSetup()
            eventAnalyzer.analyze()
            
        eventAnalyzer.end(sample)
        #######################################################




    def _runSampleLogged(analysis, sample, logDir):
        """Eseguita nel worker: manda stdout/stderr (Python e C++) su un log dedicato."""
        os.makedirs(logDir, exist_ok=True)
        logPath = os.path.join(logDir, f"{sample}.log")

        sys.stdout.flush(); sys.stderr.flush()
        saved = os.dup(1), os.dup(2)
        try:
            with open(logPath, "w", buffering=1) as log:
                os.dup2(log.fileno(), 1)
                os.dup2(log.fileno(), 2)
                try:
                    analysis.analyzeSample(sample)
                except Exception:
                    traceback.print_exc()   # il traceback completo finisce nel log
                    raise
                finally:
                    sys.stdout.flush(); sys.stderr.flush()
        finally:
            os.dup2(saved[0], 1); os.dup2(saved[1], 2)
            os.close(saved[0]); os.close(saved[1])
        return logPath

    def loop(self, nWorkers=None):  
        # ## Loop over the samples
        # for sample in self.samples:
        #     self.analyzeSample(sample)
        print(f"Number of samples to be processed: {len(self.samples)}")
        print(f"Number of workers set by user: {nWorkers}")        
        ## Loop over the samples, analyzing up to nWorkers of them in parallel
        if nWorkers is None:
            print(f"workers not set, checking how many CPU are present in the system. Number of CPU detected: {os.cpu_count()}")
            nWorkers = os.cpu_count()
        nWorkers = max(1, min(nWorkers, len(self.samples)))               
        print(f"Using {nWorkers} workers to process the samples")
        
        # Serial fallback: comodo per il debug e per i traceback leggibili
        if nWorkers == 1:
            for sample in self.samples:
                self.analyzeSample(sample)
            return

        failed = []
        with ProcessPoolExecutor(max_workers=nWorkers) as executor:
            futures = {executor.submit(self.analyzeSample, sample): sample
                       for sample in self.samples}
            for future in as_completed(futures):
                sample = futures[future]
                try:
                    future.result()
                    print(f"[done]   {sample}")
                except Exception as e:
                    print(f"[FAILED] {sample}: {e}")
                    failed.append(sample)

        if failed:
            raise RuntimeError(f"Analysis failed for {len(failed)} sample(s): {failed}")






            
            
    def end(self):
        self.outFile.Close()


         
