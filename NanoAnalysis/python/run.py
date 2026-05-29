#!/bin/env python3

from VVXAnalysis.NanoAnalysis.SampleLooper import SampleLooper 


if __name__ == "__main__" :

    ## To be fixed
    dataType = 'MC'
    
    sampleLooper = SampleLooper(dataType)
    sampleLooper.loop()
    sampleLooper.end()
