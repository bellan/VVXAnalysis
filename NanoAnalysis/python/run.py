#!/bin/env python3

from VVXAnalysis.NanoAnalysis.SampleLooper import SampleLooper

from VVXAnalysis.NanoAnalysis.AnalysisConfig import AnalysisConfig
from VVXAnalysis.NanoAnalysis.Colours import * 

import yaml
import json

from VVXAnalysis.NanoAnalysis.SampleLoader import SampleLoader

if __name__ == "__main__" :


    ## Read the configuration in Pydantic mode
    with open('configurations/VVXAnalyzer.yaml') as f:
        cfg = AnalysisConfig(**yaml.safe_load(f))

    sampleLoader = SampleLoader(cfg.samples) # --> check against data/samples_DB.json
    samples = sampleLoader.load()
    
    ## Banner ##
    print(f"""

    \t\t\t{White('*** UNITO Framework ***')}
    cfg: {cfg}
    Analyzer: {cfg.analysis.analyzer}
    Chosen selection for samples: {cfg.samples}
    Actual sample selection:  {samples}
    """)
    
    for s in samples:
        print(s.isMC())

    print(2)
    sampleLooper = SampleLooper(cfg.analysis, samples)
    sampleLooper.loop()
    #sampleLooper.end()
    print(3)
