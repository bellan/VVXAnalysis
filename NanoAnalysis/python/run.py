#!/bin/env python3

from optparse import OptionParser

from VVXAnalysis.NanoAnalysis.SampleLooper import SampleLooper

from VVXAnalysis.NanoAnalysis.AnalysisConfig import AnalysisConfig
from VVXAnalysis.NanoAnalysis.Colours import * 

import yaml
import json

from VVXAnalysis.NanoAnalysis.SampleLoader import SampleLoader

if __name__ == "__main__" :

    parser = OptionParser(usage="usage: %prog <analysis> [options]")
   
    parser.add_option("-j", "--jobs", dest="nJobs",
                      type='int',
                      default=None,
                      help=f"Set number of jobs. Default is None, that means it will use all CPU in the systems")
    

    (options, args) = parser.parse_args()
    analysis       = args[0]
    
    ## Read the configuration in Pydantic mode
    with open(f'configurations/{analysis}.yaml') as f:
        cfg = AnalysisConfig(**yaml.safe_load(f))

    sampleLoader = SampleLoader(cfg.samples) # --> check against data/samples_DB.json
    samples = sampleLoader.load()
    
    ## Banner ##
    print(f"""

    \t\t\t{White('*** UNITO Framework ***')}
    cfg: {cfg}
    Analyzer: {cfg.analysis.analyzer}
    Chosen selection for samples: {cfg.samples}
    Actual samples selection:  {samples}
    """)
    
    sampleLooper = SampleLooper(cfg.analysis, samples)
    sampleLooper.loop(options.nJobs)
    #sampleLooper.end()
