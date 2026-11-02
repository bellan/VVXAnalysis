#!/bin/env python3

from optparse import OptionParser

from VVXAnalysis.NanoAnalysis.SampleLooper import SampleLooper

from VVXAnalysis.NanoAnalysis.AnalysisConfig import AnalysisConfig
from VVXAnalysis.NanoAnalysis.Colours import * 

import yaml
import json
import os,sys

from VVXAnalysis.NanoAnalysis.SampleLoader import SampleLoader

if __name__ == "__main__" :

    parser = OptionParser(usage="usage: %prog <analysis> [options]")
   
    parser.add_option("-j", "--jobs", dest="nJobs",
                      type='int',
                      default=None,
                      help=f"Set number of jobs. Default is None, that means it will use all CPU in the systems")

    parser.add_option("-c", "--condor", dest="condor",
                  action="store_true",
                  default=False,
                  help="submit the jobs through Condor")

    parser.add_option("-e", "--eos", dest="eos",
                      action="store_true",
                      default=False,
                      help="use the location of the samples written in the DB")

    parser.add_option("-s", "--sample", dest="selectedSample", nargs=2,
                      metavar=("YEAR", "NAME"), default=None, action="append",
                      help="Analyze just this sample (can be repeated and used by Condor jobs")

    parser.add_option("--chunk", dest="chunks", type="int", action="append",
                      metavar="N", default=None,
                      help="Process only this chunk of a sample split in chunks (can be repeated). "
                           "Requires exactly one -s YEAR NAME. Used by the Condor jobs and to re-run failed chunks")

    parser.add_option("--merge", dest="merge",
                      action="store_true",
                      default=False,
                      help="hadd the chunk outputs of the selected samples (after the jobs are finished)")

    parser.add_option("--keep-chunks", dest="keepChunks",
                      action="store_true",
                      default=False,
                      help="debug: do not remove the outputs of the single chunks after the merge")

    parser.add_option("--update-sumw", dest="updateSumw",
                      action="store_true",
                      default=False,
                      help="recompute the total genEventSumw of the chunked MC samples, store it in the DB and exit")

    flavours = ["espresso", "microcentury", "longlunch", "workday", "tomorrow", "testmatch", "nextweek"]
    parser.add_option("-f","--flavour", dest="flavour", type="choice",
                      choices= flavours,
                      metavar="{" + ",".join(flavours) + "}",
                      default="longlunch",
                      help="JobFlavour Condor. Default is (%default)")

    parser.add_option("--dry-run", dest="dryRun",
                      action="store_true",
                      help="write Condor files without submitting them")

    (options, args) = parser.parse_args()
    if not args:
        parser.error("the analysis name is missing")
    analysis       = args[0]

    if options.chunks and (not options.selectedSample or len(options.selectedSample) != 1):
        parser.error("--chunk can be used only together with exactly one -s YEAR NAME")
    
    ## Read the configuration in Pydantic mode
    with open(f'configurations/{analysis}.yaml') as f:
        cfg = AnalysisConfig(**yaml.safe_load(f))

            
    sampleLoader = SampleLoader(cfg.samples) # --> check against data/samples_DB.json
    samples = sampleLoader.load(not options.eos)

    for s in samples:
        print(s)
    
    if options.selectedSample:
        wanted = {(int(y), n) for y, n in options.selectedSample}
        samples = [s for s in samples if (s.year, s.name) in wanted]

        found = {(s.year, s.name) for s in samples}
        missing = wanted - found

        if missing:
            if not found:
                raise SystemExit(f"None of the sample in the list is present in the YAML selection: {sorted(missing)}")
            else:
                print(f"Some of the samples in the list are not resent in the YAML selection: {sorted(missing)}. Processing only the one that are aslo in the YALM.")
        
                      
    ## Banner ##
    print(f"""

    \t\t\t{White('*** UNITO Framework ***')}
    cfg: {cfg}
    Analyzer: {cfg.analysis.analyzer}
    Chosen selection for samples: {cfg.samples}
    Actual samples selection:  {samples}
    """)
    
    sampleLooper = SampleLooper(cfg.analysis, samples,
                                loader=sampleLoader,
                                keepChunks=options.keepChunks)

    if options.updateSumw:
        sampleLooper.updateSumw()
    elif options.merge:
        sampleLooper.mergeAll()
    elif options.condor:
        ## To be fixed
        args = [analysis]
        ## With condor, run on eos samples by default
        #if options.eos:
        args.append("-e")
        sampleLooper.submitCondor(runScript=os.path.abspath(sys.argv[0]),
                                  args=args,
                                  flavour=options.flavour,
                                  chunks=options.chunks,
                                  dryRun=options.dryRun)
    else:
        sampleLooper.loop(options.nJobs, chunks=options.chunks)
