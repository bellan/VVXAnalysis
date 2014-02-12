#! /usr/bin/env python


##################################
## R. Bellan (UNITO) - Feb 2014 ##
##################################


import sys, os, commands, math

from readSampleInfo import *

print 'ciao'

typeofsample = sys.argv[1]
cregion = 'none' #sys.argv[2]
typeofsamples = ['test','mudata', 'edata', 'W', 'Z', 'ttbar','QCDPT', 'diboson','ttZ']

baseinputdir = '/afs/cern.ch/work/b/bellan/public/samples/'

def run(typeofsample, cregion):
    inputdir  = '../samples/newvars20130711/'
    outputdir = 'output'
    
    executable = 'bin/eventAnalyzer'
    lumi = 19029.853

    runperiods = []
    sample = ''

    outputdir = outputdir+"_"+cregion+"/"
    if not os.path.exists(outputdir):
        os.popen('mkdir "%s"' %outputdir)

    if typeofsample == 'test':
        runperiods = ['1']
        inputdir   = 'samples'
        sample     = 'test_'

    if typeofsample == 'WZZ':
        runperiods = ['']
        inputdir   = baseinputdir+"WZZJets"
        sample     = 'ZZWAnalysis'

    if typeofsample == 'mudata' or typeofsample == 'edata': 
        runperiods = ['2012A-13Jul', '2012A-06Aug', '2012B-13Jul', '2012C-24Aug', '2012C-11Dec', '2012C-PromptReco', '2012D-PromptReco']
        lumi = -1
        if typeofsample == 'mudata':
            inputdir  = 'samples/newvars20130722/'
            sample = 'MuData-'
        if typeofsample == 'edata':
            sample = 'EData-'
            inputdir  = 'samples/newvars20130722/'

    if typeofsample == 'w' or typeofsample == 'W':
        runperiods = ['Wlnu1J-madgraph', 'Wlnu2J-madgraph', 'Wlnu3J-madgraph', 'Wlnu4J-madgraph']

    if typeofsample == 'tt' or typeofsample == 'ttbar':
        runperiods = ['ttbar-madgraph', 'tbarW-powheg', 'tW-powheg'] # run on powheg too
        #runperiods = ['ttbar-powheg']

    if typeofsample == 'z' or typeofsample == 'Z':
        runperiods = ['Zll1J-madgraph', 'Zll2J-madgraph', 'Zll3J-madgraph', 'Zll4J-madgraph']

    if typeofsample == 'diboson':
        runperiods = ['WW-pythia', 'WZ-pythia', 'ZZ-pythia']

    if typeofsample == 'ttZ':
        runperiods = ['ttZ-madgraph']

    if typeofsample == 'qcdht' or typeofsample == 'QCDHT':
        runperiods = ['QCD-HT100to250-madgraph', 'QCD-HT250to500-madgraph', 'QCD-HT500to1000-madgraph', 'QCD-HT1000toInf-madgraph']
    
    if typeofsample == 'qcdpt' or typeofsample == 'QCDPT':
        runperiods = ["QCD-PT0to5-pythia6",        "QCD-PT15to30-pythia6",     "QCD-PT30to50-pythia6",    "QCD-PT600to800-pythia6",
                      "QCD-PT1000to1400-pythia6",  "QCD-PT170to300-pythia6",   "QCD-PT470to600-pythia6",  "QCD-PT800to1000-pythia6",
                      "QCD-PT120to170-pythia6",    "QCD-PT1800toInf-pythia6",  "QCD-PT50to80-pythia6",    "QCD-PT80to120-pythia6",
                      "QCD-PT1400to1800-pythia6",  "QCD-PT300to470-pythia6",   "QCD-PT5to15-pythia6"]

    # ----- Run over the run periods -----
    hadd = 'hadd ' + outputdir+typeofsample + '.root'
    for period in runperiods:
        basefile = sample+period
        if os.path.exists(outputdir+basefile+'.root'):
            os.popen('rm "%s".root' %(outputdir+basefile))

        externalXsec = -1
        if not typeofsample == 'mudata' and not typeofsample == 'edata':
            externalXsec = crossSection(period)
            print period, " ---> External cross section: ", externalXsec
        command = "./{0:s} {1:s}/{3:s}.root {2:s}/{3:s}.root {4:.0f} {5:.3f}".format(executable,inputdir,outputdir, basefile, lumi, externalXsec)
        print command
        failure, output = commands.getstatusoutput(command)
        print output
        hadd = hadd + " " + outputdir+basefile+'.root'

    if os.path.exists( outputdir+typeofsample + '.root'):
        os.popen('rm "%s".root' %(outputdir+typeofsample))
    print hadd
    failure, output = commands.getstatusoutput(hadd)


if typeofsample == 'all':
    for sample in typeofsamples:
        if cregion == 'all':
            for cr in range(0,4):
                run(sample, cr)    # runs over all samples in all control reagions
        else:
            run(sample, cregion)   # runs over all samples in a specific control reagions
else:
    if cregion == 'all':
        for cr in range(0,4):     
            run(typeofsample, cr)  # runs over a specific sample in all control regions
    else:
        run(typeofsample, cregion) # runs over a specific sample in a specific region
