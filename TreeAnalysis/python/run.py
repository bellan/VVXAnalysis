#! /usr/bin/env python


##################################
## R. Bellan (UNITO) - Feb 2014 ##
##################################


print '*** UNITO Framework ***'


import sys, os, commands, math
from readSampleInfo import *




############################################################################
############################## User's inputs ###############################
############################################################################

executable = sys.argv[1]
typeofsample = sys.argv[2]

getExternalCrossSectionFromFile = False
if len(sys.argv) > 3: getExternalCrossSectionFromFile = sys.argv[3] 

cregion = 'baseline' # in case, make it a further external argument

luminosity = 19029.853
baseinputdir = 'samples/'



############################## Configuration ##############################
################ Change this if you know what you are doing ###############
###########################################################################


csvfile = '../Producers/python/samples_8TeV.csv'
DB = readSampleDB(csvfile)

typeofsamples = typeOfSamples(csvfile)
typeofsamples.append('test')

failure, output = commands.getstatusoutput('ls ./bin/ | grep -v .cpp | grep -v .xml | grep -v .md')
availableExecutable = output.split()
if executable in availableExecutable: 
    executable = 'bin/'+executable
else:
    print "ERROR! Unknown executable. Availble executables are:",availableExecutable
    sys.exit(1)


############################################################################
############################ Print the configuration #######################
############################################################################

print "Configuration:"
print "Executable:", executable
print "Sample/type of samples:", typeofsample
print "Get (again) cross section from csv file:", getExternalCrossSectionFromFile
print "CSV file:", csvfile
print "Control region type:", cregion
print "Integrated luminosity:", luminosity



############################################################################


def run(executable, typeofsample, cregion, luminosity):
    inputdir  = 'samples/'
    outputdir = 'results/'
    if not os.path.exists(outputdir): os.popen('mkdir "%s"' %outputdir)

    outputdir = outputdir+executable[4:]+"_"+cregion+"/"
    if not os.path.exists(outputdir): os.popen('mkdir "%s"' %outputdir)


    datasets = getSamplesBy('process',typeofsample,csvfile)
    if len(datasets) == 0:
        datasets = getSamplesBy('identifier',typeofsample,csvfile)

    sampleprefix = ''

    

    ### Override configuration, for test only ###
    if typeofsample == 'test':
        datasets = ['test']
        

    #################################################################################

    ### Special treatment for DATA ###

    #################################################################################

    if typeofsample == 'mudata' or typeofsample == 'edata': 
        datasets = ['2012A-13Jul', '2012A-06Aug', '2012B-13Jul', '2012C-24Aug', '2012C-11Dec', '2012C-PromptReco', '2012D-PromptReco']
        luminosity = -1
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
            externalXsec = crossSection(period, csvfile)
            print period, " ---> External cross section: ", externalXsec
        command = "./{0:s} {1:s}/{3:s}.root {2:s}/{3:s}.root {4:.0f} {5:.3f}".format(executable,inputdir,outputdir, basefile, luminosity, externalXsec)
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
                run(executable, sample, cr, luminosity)    # runs over all samples in all control reagions
        else:
            run(executable, sample, cregion, luminosity)   # runs over all samples in a specific control reagions
else:
    if cregion == 'all':
        for cr in range(0,4):     
            run(executable, typeofsample, cr, luminosity)  # runs over a specific sample in all control regions
    else:
        run(executable, typeofsample, cregion, luminosity) # runs over a specific sample in a specific region
