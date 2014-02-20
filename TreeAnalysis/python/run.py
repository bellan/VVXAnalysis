#! /usr/bin/env python


##################################
## R. Bellan (UNITO) - Feb 2014 ##
##################################


import sys, os, commands, math
from readSampleInfo import *
from Colours import *


print "\n\n"
print "\t\t\t",White('*** UNITO Framework ***')
print "\n\n"


############################################################################
############################## User's inputs ###############################
############################################################################

executable = sys.argv[1]
typeofsample = sys.argv[2]

getExternalCrossSectionFromFile = False
if len(sys.argv) > 3: getExternalCrossSectionFromFile = sys.argv[3] 

cregion = 'baseline' # in case, make it a further external argument

luminosity = 19029.853
baseinputdir = 'samples'



############################## Configuration ##############################
################ Change this if you know what you are doing ###############
###########################################################################


csvfile = '../Producers/python/samples_8TeV.csv'

typeofsamples = typeOfSamples(csvfile)
typeofsamples.append('test')

failure, output = commands.getstatusoutput('ls ./bin/ | grep -v .cpp | grep -v .xml | grep -v .md')
availableExecutable = output.split()
if executable in availableExecutable: 
    executable = 'bin/'+executable
else:
    print Important("ERROR! Unknown executable."),"Availble executables are:",availableExecutable
    sys.exit(1)


############################################################################
############################ Print the configuration #######################
############################################################################

print Blue("----------------------------------------------------------------------")
print Red("Configuration:")
print "Executable:", Blue(executable)
print "Sample/type of samples:", Blue(typeofsample)
print "CSV file: ", Blue(csvfile)
print "Get (again) cross sections from csv file: ", Blue(getExternalCrossSectionFromFile)
print "Control region type: ", Blue(cregion)
print "Integrated luminosity: ", Blue(luminosity)
print Blue("----------------------------------------------------------------------")
print "\n"


############################################################################


def run(executable, typeofsample, cregion, luminosity):
    inputdir  = baseinputdir
    outputdir = 'results'
    if not os.path.exists(outputdir): os.popen('mkdir "%s"' %outputdir)

    outputdir = outputdir+"/"+executable[4:]+"_"+cregion
    if not os.path.exists(outputdir): os.popen('mkdir "%s"' %outputdir)


    datasets = getSamplesBy('process',typeofsample,csvfile)
    if len(datasets) == 0:
        datasets = getSamplesBy('identifier',typeofsample,csvfile)
    if len(datasets) == 0:
        print Important('Error! This sample is not available!'), typeofsample

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
            sampleprefix = 'MuData-'
        if typeofsample == 'edata':
            sampleprefix = 'EData-'


    # ----- Run over the run periods -----
    hadd = 'hadd {0:s}/{1:s}.root'.format(outputdir,typeofsample)
    for period in datasets:
        basefile = sampleprefix+period
        if os.path.exists('{0:s}/{1:s}.root'.format(outputdir,basefile)):
            os.popen('rm {0:s}/{1:s}.root'.format(outputdir,basefile))

        externalXsec = -1
        if not typeofsample == 'mudata' and not typeofsample == 'edata' and getExternalCrossSectionFromFile:
            externalXsec = crossSection(period, csvfile)
            print "For {0:s} {1:s} {2:.6f}".format(period, Warning("Using external cross section:"), externalXsec)

        print Red('\n------------------------------ {0:s} -------------------------------\n'.format(basefile))
        command = "./{0:s} {1:s}/{3:s}.root {2:s}/{3:s}.root {4:.0f} {5:.10f}".format(executable,inputdir,outputdir, basefile, luminosity, externalXsec)
        print "Command going to be executed:", Violet(command)
        failure, output = commands.getstatusoutput(command)
        print "\n",output
        hadd = '{0:s} {1:s}/{2:s}.root'.format(hadd, outputdir, basefile)

    print Red('----------------------------------------------------------------------\n')
    if len(datasets) > 1:
        if os.path.exists('{0:s}/{1:s}.root'.format(outputdir,typeofsample)):
            os.popen('rm {0:s}/{1:s}.root'.format(outputdir,typeofsample))
        print "Command going to be executed:", Violet(hadd)
        failure, output = commands.getstatusoutput(hadd)

    print "The output is in", Green('{0:s}/{1:s}.root'.format(outputdir,typeofsample))  


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

print "\nJob status: ", OK("DONE"),"\n"
