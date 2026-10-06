#! /usr/bin/env python

import sys, os, subprocess

print("Merging ZZ samples")

outputdir = sys.argv[1]
inputdir  = outputdir+'ZZunmerged'

print(inputdir)

failure, output = subprocess.getstatusoutput('ls {0:s}/*.root | grep -v ext'.format(inputdir))

samples = output.split()
for sample in samples:
    failure, name = subprocess.getstatusoutput('basename {0:s} .root'.format(sample))
    command =  'hadd {0:s}/{1:s}.root {2:s}/{1:s}.root {2:s}/{1:s}_ext.root'.format(outputdir, name, inputdir)
    print(command)
    failure, output = subprocess.getstatusoutput(command)
