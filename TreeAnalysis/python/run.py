#! /usr/bin/env python


##################################
## R. Bellan (UNITO) - Feb 2014 ##
##################################


print 'Ciao'



import sys, os, commands, math
from readSampleInfo import *



############################## Configuration ##############################
typeofsamples = typeOfSamples('../Producers/python/samples_8TeV.csv')
typeofsamples.append('test')

DB = readSampleDB('../Producers/python/samples_8TeV.csv')

typeofsample = sys.argv[1]

baseinputdir = '/afs/cern.ch/work/b/bellan/public/samples/'

cregion = 'baseline' #sys.argv[2]

executable = 'bin/eventAnalyzer'

getExternalCrossSectionFromFile = False
if len(sys.argv) > 2:
    getExternalCrossSectionFromFile = sys.argv[2] 


############################################################################




def run(executable, typeofsample, cregion):
    inputdir  = 'samples/'
    outputdir = 'output'
    
    lumi = 19029.853

    datasets = getSamplesBy('process',typeofsample,'../Producers/python/samples_8TeV.csv')

    sampleprefix = ''

    outputdir = outputdir+"_"+cregion+"/"
    if not os.path.exists(outputdir):
        os.popen('mkdir "%s"' %outputdir)

    ### Override configuration, for test only ###
    if typeofsample == 'test':
        datasets = ['test']
        

    #################################################################################

    ### DATA ###

    #################################################################################

    if typeofsample == 'mudata' or typeofsample == 'edata': 
        datasets = ['2012A-13Jul', '2012A-06Aug', '2012B-13Jul', '2012C-24Aug', '2012C-11Dec', '2012C-PromptReco', '2012D-PromptReco']
        lumi = -1
        if typeofsample == 'mudata':
            inputdir  = 'samples/'
            sampleprefix = 'MuData-'
        if typeofsample == 'edata':
            sampleprefix = 'EData-'
            inputdir  = 'samples/'


    # ----- Run over the run periods -----
    hadd = 'hadd ' + outputdir+typeofsample + '.root'
    for period in datasets:
        basefile = sampleprefix+period
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
                run(executable, sample, cr)    # runs over all samples in all control reagions
        else:
            run(executable, sample, cregion)   # runs over all samples in a specific control reagions
else:
    if cregion == 'all':
        for cr in range(0,4):     
            run(executable, typeofsample, cr)  # runs over a specific sample in all control regions
    else:
        run(executable, typeofsample, cregion) # runs over a specific sample in a specific region
