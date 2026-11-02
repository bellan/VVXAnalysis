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

import os, sys, traceback, shutil
from concurrent.futures import ProcessPoolExecutor, as_completed

import subprocess

import ctypes
_libc = ctypes.CDLL(None)
    
maxEntriesPerSample = None # Use only up to this number of events in each MC sample, for quick tests; use None for no scaling


def computeTotalSumw(sample):
    """Sum of genEventSumw over ALL the chunks of a sample (only the Runs tree is read)"""
    chunks = sample.chunkFiles()
    values = []
    for i, (n, path) in enumerate(chunks, 1):
        f = ROOT.TFile.Open(path)
        if not f or f.IsZombie():
            raise IOError(f"Cannot open {path}")
        values.append(get_genEventSumw(f))
        f.Close()
        print(f"  [{i}/{len(chunks)}] Chunk{n}: sumw = {values[-1]}")
    return math.fsum(values), [n for n, _ in chunks]


class SampleLooper:
    
    def __init__(self, cfg, samples, loader=None, keepChunks=False):

        self.analyzer   = cfg.analyzer
        self.regions    = cfg.regions
        self.samples    = samples
        self.loader     = loader       # needed only to write genEventSumw in the DB
        self.keepChunks = keepChunks   # debug: keep the chunk outputs after the merge

        self.load_analyses()
        

    def load_analyses(self):           
        for _, module_name, _ in pkgutil.iter_modules(
                analyses.__path__
        ):           
            importlib.import_module(
                f"analyses.{module_name}"
            )


    # ==================================================================
    # Tasks: one per sample (single file) or one per chunk
    # ==================================================================
    def buildTasks(self, chunks=None):
        tasks = []
        for s in self.samples:
            try:
                tasks.extend(s.tasks(chunks))
            except FileNotFoundError as e:
                print(f"[SKIPPED] {s}: {e}")
        return tasks


    def chunkedSamples(self, tasks):
        return list(dict.fromkeys(t.sample for t in tasks if t.isChunk()))


    def logPath(self, task):
        return f"logs/{self.analyzer}/{task.sample.year}/{task.logStem()}.log"


    # ==================================================================
    # Total genEventSumw of the chunked MC samples (cached in the DB)
    # ==================================================================
    def prepareSumw(self, tasks, allowCompute=True, force=False):
        for s in self.chunkedSamples(tasks):
            if not s.isMC():
                continue
            if not force and s.validSumw() is not None:
                continue
            if not allowCompute:
                raise RuntimeError(
                    f"genEventSumw of {s} is missing in the DB or computed on a different set of chunks. "
                    "Run once without --chunk (or with --update-sumw) to compute it.")
            if self.loader is None:
                raise RuntimeError("A SampleLoader is needed to store genEventSumw in the DB")
            print(f"Computing the total genEventSumw of {s}")
            value, chunkIds = computeTotalSumw(s)
            self.loader.store_sumw(s, value, len(chunkIds))
            print(f"  total genEventSumw = {value} ({len(chunkIds)} chunks), stored in the DB")


    def updateSumw(self):
        """--update-sumw: recompute and store genEventSumw of the chunked MC samples, then stop"""
        self.prepareSumw(self.buildTasks(), force=True)


    # ==================================================================
    # Manifest: what the merge has to find (analyzer, regions, chunks)
    # ==================================================================
    def writeManifests(self, tasks):
        for s in self.chunkedSamples(tasks):
            staging = s.stagingDir()
            if os.path.isdir(staging):
                print(f"Removing old chunk outputs in {staging}")
                shutil.rmtree(staging)
            os.makedirs(staging)
            manifest = {"analyzer": self.analyzer,
                        "regions" : list(self.regions),
                        "year"    : s.year,
                        "origin"  : s.origin,
                        "name"    : s.name,
                        "chunks"  : sorted(t.chunk for t in tasks if t.sample is s)}
            with open(os.path.join(staging, "manifest.json"), "w") as f:
                json.dump(manifest, f, indent=1)


    # ==================================================================
    # The analysis of one task
    # ==================================================================
    def analyzeSample(self, task):
        sample = task.sample
        print("analyzing", task.label(), "-", sample.explainYourself())
        print(task.inputPath)
        inputFile = ROOT.TFile.Open(task.inputPath)
        if not inputFile or inputFile.IsZombie():
            raise IOError(f"Cannot open {task.inputPath}")
        
        '''
        Set the Gen Event Sum Weight, needed to properly normalize the sample weight
        for data it is always 1. For MC in a single file we extract it from the counters.
        For a chunk it must be the sum over ALL the chunks of the sample (stored in the DB),
        otherwise each chunk would be normalized to its own sum of weights.
        We then need to pass this information to the Analyzer
        '''
        if not sample.isMC():
            genEventSumw = 1.
        elif task.isChunk():
            genEventSumw = sample.validSumw()
            if genEventSumw is None:
                raise RuntimeError(f"genEventSumw of {sample} is missing or stale in the DB, "
                                   "run with --update-sumw")
        else:
            genEventSumw = get_genEventSumw(inputFile, maxEntriesPerSample)
        print(f"genEventSumw used for the normalization: {genEventSumw}")
        
        events = inputFile.Get("Events")
        nEntries = events.GetEntries()
        iEntry=0
        printEntries=max(5000,nEntries//10)
        
        ######### Analyse the events in a sample! #############
        eventAnalyzer = EventAnalyzer.registry[self.analyzer](self.regions)#(base_configuration)
        eventAnalyzer.init(events, genEventSumw, sample.luminosity, sample.isMC())
        eventAnalyzer.begin()
        
        while iEntry<nEntries and events.GetEntry(iEntry):
            iEntry+=1
            if iEntry%printEntries == 0 : print("Processing", iEntry)
            eventAnalyzer.eventSetup()
            eventAnalyzer.analyze()
        
        eventAnalyzer.end(sample, task.outputBase())
        eventAnalyzer.endOfAnalysis()
        #######################################################


    def _flushAll(self):
        sys.stdout.flush(); sys.stderr.flush()
        _libc.fflush(None)          #flush all stdio C (ROOT included) buffers


    def _runSampleLogged(self, task):
        """stdout/stderr (Python e C++) in a dedicate log."""
        logPath = self.logPath(task)
        os.makedirs(os.path.dirname(logPath), exist_ok=True)
        print(f"log is in {logPath}")
        self._flushAll()
        saved = os.dup(1), os.dup(2)
        try:
            with open(logPath, "w", buffering=1) as log:
                os.dup2(log.fileno(), 1)
                os.dup2(log.fileno(), 2)
                try:
                    self.analyzeSample(task)
                except Exception:
                    traceback.print_exc()   # traceback in the log
                    raise
                finally:
                    self._flushAll()
        finally:
            os.dup2(saved[0], 1); os.dup2(saved[1], 2)
            os.close(saved[0]); os.close(saved[1])
        return logPath


    # ==================================================================
    # Local backend
    # ==================================================================
    def loop(self, nWorkers=None, chunks=None, merge=True):
        tasks = self.buildTasks(chunks)
        if not tasks:
            print("No task to process")
            return

        # Total sumw in the DB (never computed when only some chunks are requested)
        self.prepareSumw(tasks, allowCompute=(chunks is None))
        # A new full run starts from a clean staging area
        if chunks is None:
            self.writeManifests(tasks)

        print(f"Number of tasks to be processed (samples or chunks): {len(tasks)}")
        print(f"Number of workers set by user: {nWorkers}")        
        ## Loop over the tasks, analyzing up to nWorkers of them in parallel
        if nWorkers is None:
            print(f"workers not set, checking how many CPU are present in the system. Number of CPU detected: {os.cpu_count()}")
            nWorkers = os.cpu_count()
        nWorkers = max(1, min(nWorkers, len(tasks)))               
        print(f"Using {nWorkers} workers to process the tasks")
        
        failed = []
        # Serial fallback: comodo per il debug e per i traceback leggibili
        if nWorkers == 1:
            for task in tasks:
                self.analyzeSample(task)
        else:
            total, nDone = len(tasks), 0
            with ProcessPoolExecutor(max_workers=nWorkers) as executor:
                futures = {executor.submit(self._runSampleLogged, task): task
                           for task in tasks}
                for future in as_completed(futures):
                    task = futures[future]
                    nDone +=1
                    try:
                        future.result()
                        print(f"[{nDone}/{total}] done    {task.label()}")
                    except Exception as e:
                        print(f"[{nDone}/{total}] FAILED  {task.label()}: {e}  -> {self.logPath(task)}")
                        failed.append(task)

        # Merge the chunks of the samples that were completely processed in this run
        # (with --chunk we are re-running a few chunks: merge later with --merge)
        if merge and chunks is None:
            failedSamples = {t.sample for t in failed}
            for s in self.chunkedSamples(tasks):
                if s in failedSamples:
                    print(f"{s}: some chunks failed, merge skipped (fix them and use --merge)")
                    continue
                self.merge(s)

        if failed:
            raise RuntimeError(f"Analysis failed for {len(failed)} task(s): {[t.label() for t in failed]}")
         

    # ==================================================================
    # Merge of the chunk outputs
    # ==================================================================
    def merge(self, sample):
        """hadd the chunk outputs of a sample, region by region, and clean the staging area.
        Returns True if the merge was completed."""
        staging = sample.stagingDir()
        manifestFile = os.path.join(staging, "manifest.json")
        if not os.path.isfile(manifestFile):
            print(f"{sample}: no {manifestFile}, nothing to merge")
            return False
        with open(manifestFile) as f:
            m = json.load(f)

        def chunkOutput(n, region):
            return os.path.join(staging, f"Chunk{n}", f"{m['analyzer']}_{region}", f"{m['name']}.root")

        # Never merge an incomplete set of chunks
        missing = [(n, r) for n in m["chunks"] for r in m["regions"]
                   if not os.path.isfile(chunkOutput(n, r))]
        if missing:
            chunksMissing = sorted({n for n, _ in missing})
            print(f"{sample}: merge NOT done, outputs missing for chunk(s) {chunksMissing}")
            return False

        for region in m["regions"]:
            outDir = f"results/{m['year']}/{m['analyzer']}_{region}"
            os.makedirs(outDir, exist_ok=True)
            final = f"{outDir}/{m['name']}.root"
            tmp   = f"{outDir}/.{m['name']}.tmp.root"      # renamed only if hadd succeeds
            inputs = [chunkOutput(n, region) for n in m["chunks"]]
            print(f"{sample}: merging {len(inputs)} chunks -> {final}")
            try:
                subprocess.run(["hadd", "-f", "-n", "100", tmp] + inputs, check=True)
            except Exception:
                if os.path.exists(tmp):
                    os.remove(tmp)
                raise
            os.replace(tmp, final)

        if self.keepChunks:
            print(f"{sample}: chunk outputs kept in {staging} (--keep-chunks)")
        else:
            shutil.rmtree(staging)
        return True


    def mergeAll(self):
        """--merge: merge all the selected chunked samples that have outputs in the staging area"""
        nFound, notDone = 0, []
        for s in self.samples:
            if not os.path.isfile(os.path.join(s.stagingDir(), "manifest.json")):
                continue
            nFound += 1
            if not self.merge(s):
                notDone.append(s)
        if nFound == 0:
            print("No chunk outputs to merge for the selected samples")
        if notDone:
            raise SystemExit(f"Merge not completed for: {notDone}")
         
            
    def end(self):
        self.outFile.Close()


    # ==================================================================
    # Condor backend: one job per sample (single file) or per chunk
    # ==================================================================
    def submitCondor(self, runScript, args, flavour="longlunch", chunks=None, dryRun=False):
        logDir=f"logs/{args[0]}"
        tasks = self.buildTasks(chunks)
        if not tasks:
            print("No task to submit")
            return

        if dryRun:
            for s in self.chunkedSamples(tasks):
                if s.isMC() and s.validSumw() is None:
                    print(f"(dry run) the total genEventSumw of {s} would be computed and stored in the DB")
        else:
            # Done here, on the submission machine: the jobs only read the value from the DB
            self.prepareSumw(tasks, allowCompute=(chunks is None))
            if chunks is None:
                self.writeManifests(tasks)

        for t in tasks:
            os.makedirs(os.path.dirname(f"{logDir}/{t.sample.year}/{t.logStem()}"), exist_ok=True)

        cwd = os.getcwd()

        # Wrapper: ricrea l'ambiente CMSSW sul nodo ed esegue run.py per un sample (o un chunk)
        wrapper = os.path.join(cwd, "condor_wrapper.sh")
        with open(wrapper, "w") as f:
            f.write(f"""#!/bin/bash
set -e
source /cvmfs/cms.cern.ch/cmsset_default.sh
cd {os.environ['CMSSW_BASE']}/src
eval `scramv1 runtime -sh`
cd {cwd}
CHUNK=""
if [ "$3" != "all" ]; then CHUNK="--chunk $3"; fi
python3 {runScript} {' '.join(args)} --sample "$1" "$2" $CHUNK
""")
        os.chmod(wrapper, 0o755)
            
        # year, name, chunk ("all" = whole file), log path (without extension, relative to logs/<analysis>/<year>/)
        with open("samples.txt", "w") as f:
            for t in tasks:
                chunk = t.chunk if t.isChunk() else "all"
                f.write(f"{t.sample.year} {t.sample.name} {chunk} {t.logStem()}\n")
                
        with open("analysis.sub", "w") as f:
            f.write(f"""executable            = {wrapper}
arguments             = $(sampleYear) $(sampleName) $(sampleChunk)
output                = {logDir}/$(sampleYear)/$(logStem).out
error                 = {logDir}/$(sampleYear)/$(logStem).err
log                   = {logDir}/condor.log
should_transfer_files = NO
+JobFlavour           = "{flavour}"
queue sampleYear, sampleName, sampleChunk, logStem from samples.txt
""")
                
        if dryRun:
            print(f"condor_wrapper.sh, samples.txt, analysis.sub (dry run) created: {len(tasks)} job(s)")
            return
        subprocess.run(["condor_submit", "analysis.sub"], check=True)
        print(f"{len(tasks)} job(s) submitted. When they are done: run.py {args[0]} --merge")
