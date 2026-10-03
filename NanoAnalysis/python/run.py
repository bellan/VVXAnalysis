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


    
    parser.add_option("-s", "--sample",
                      dest="selectedSample",
                      action="append", default=None,
                      help="Analyze just this sample (can be ripeted and used by Condor jobs")
    
    parser.add_option("-f","--flavour", dest="flavour",
                      default="longlunch",
                      help="JobFlavour Condor")

    parser.add_option("--dry-run", dest="dryRun",
                      action="store_true",
                      help="write Condor files without submitting them")

    (options, args) = parser.parse_args()
    analysis       = args[0]
    
    ## Read the configuration in Pydantic mode
    with open(f'configurations/{analysis}.yaml') as f:
        cfg = AnalysisConfig(**yaml.safe_load(f))

            
    sampleLoader = SampleLoader(cfg.samples) # --> check against data/samples_DB.json
    samples = sampleLoader.load(not options.eos)

    for s in samples:
        print(s)
    
    if options.selectedSample:
        samples = [s for s in samples if s.name in options.selectedSample]

    
    ## Banner ##
    print(f"""

    \t\t\t{White('*** UNITO Framework ***')}
    cfg: {cfg}
    Analyzer: {cfg.analysis.analyzer}
    Chosen selection for samples: {cfg.samples}
    Actual samples selection:  {samples}
    """)
    
    sampleLooper = SampleLooper(cfg.analysis, samples)
    #sampleLooper.loop(options.nJobs)
    ##sampleLooper.end()


    if options.condor:
        ## To be fixed
        args = [analysis]
        #if options.eos:
        args.append("-e")
        
        sampleLooper.submitCondor(runScript=os.path.abspath(sys.argv[0]),
                                  args=args,
                                  flavour=options.flavour,
                                  dryRun=options.dryRun)
    else:
        sampleLooper.loop(options.nJobs)
