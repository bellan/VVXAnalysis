#! /usr/bin/env python


##################################
## R. Bellan (UNITO) - Feb 2014 ##
##################################


import sys, os, commands, math

from readSampleInfo import *

print 'ciao'

typeofsample = sys.argv[1]
cregion = 'baseline' #sys.argv[2]

getExternalCrossSectionFromFile = False
if len(sys.argv) > 2:
    getExternalCrossSectionFromFile = sys.argv[2] 

typeofsamples = ['test','mudata', 'edata', 'W', 'Z', 'ttbar','QCDPT', 'diboson','ttZ']

baseinputdir = '/afs/cern.ch/work/b/bellan/public/samples/'

def run(typeofsample, cregion):
    inputdir  = 'samples/'
    outputdir = 'output'
    
    executable = 'bin/eventAnalyzer'
    lumi = 19029.853

    runperiods = []
    sample = ''

    outputdir = outputdir+"_"+cregion+"/"
    if not os.path.exists(outputdir):
        os.popen('mkdir "%s"' %outputdir)
        
    if typeofsample == 'test':
        runperiods = ['test']
        
    #################################################################################

    ### SIGNAL ###
    if typeofsample == 'WZZ':
        runperiods = ['WZZJets']

    #################################################################################

    ### BACKGROUNDS ###
    if  typeofsample == 'ZZ':
        runperiods = ['ZZ2e2mu', 'ZZ2mu2tau', 'ZZ4mu', 'ZZ2e2tau', 'ZZ4e', 'ZZ4tau', 'ZZ2e2mu_ext', 'ZZ2mu2tau_ext', 'ZZ4mu_ext', 'ZZ2e2tau_ext', 'ZZ4e_ext', 'ZZ4tau_ext','ggZZ4l','ggZZ2l2l']
        
    if typeofsample == 'ZZJetsTo4L':
        runperiods = ['ZZJetsTo4L']

    if typeofsample == 'tt' or typeofsample == 'ttbar':
        runperiods = ['TTTo2L2Nu2B','TTZJets','TTWJets','TTWWJets']

    if typeofsample == 'H':
        runperiods = ['powheg15H126','VBFH126','ttH126','WH126','ZH126']
        
    if typeofsample == 'triboson':
        runperiods = ['ZZZJets','WWWJets','WWZJets']

    if typeofsample == 'WZ':
        runperiods = ['WZ']

    if typeofsample == 'Z' or typeofsample == 'z':
        runperiods = ['DYJetsToLLTuneZ2M10','DYJetsToLLTuneZ2M50']

    #missing: powheg15jhuGenV3H126, minloH126, WZZ_aMCatNLO

    #################################################################################

    ### INDIVIDUAL SAMPLES, FOR CONVENIENCE ###

    if typeofsample == 'ZZZ':
        runperiods = ['ZZZJets']

    if typeofsample == 'ZZ2e2mu':
        runperiods = ['ZZ2e2mu']

    if typeofsample == 'ZZ2m2tau':
        runperiods = ['ZZ2mu2tau']

    if typeofsample == 'ZZ2e2tau':
        runperiods = ['ZZ2e2tau']

    if typeofsample == 'ZZ4mu':
        runperiods = ['ZZ4mu']

    if typeofsample == 'ZZ4e':
        runperiods = ['ZZ4e']

    if typeofsample == 'ZZ4tau':
        runperiods = ['ZZ4tau']

    if typeofsample == 'TTZ':
        runperiods = ['TTZJets']

    #################################################################################

    ### DATA ###

    #################################################################################

    if typeofsample == 'mudata' or typeofsample == 'edata': 
        runperiods = ['2012A-13Jul', '2012A-06Aug', '2012B-13Jul', '2012C-24Aug', '2012C-11Dec', '2012C-PromptReco', '2012D-PromptReco']
        lumi = -1
        if typeofsample == 'mudata':
            inputdir  = 'samples/'
            sample = 'MuData-'
        if typeofsample == 'edata':
            sample = 'EData-'
            inputdir  = 'samples/'



    # ----- Run over the run periods -----
    hadd = 'hadd ' + outputdir+typeofsample + '.root'
    for period in runperiods:
        basefile = sample+period
        if os.path.exists(outputdir+basefile+'.root'):
            os.popen('rm "%s".root' %(outputdir+basefile))

        externalXsec = -1
        if not typeofsample == 'mudata' and not typeofsample == 'edata' and getExternalCrossSectionFromFile:
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
