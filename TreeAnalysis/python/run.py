#! /usr/bin/env python


##################################
## R. Bellan (UNITO) - Feb 2014 ##
##################################


import sys, os, commands, math
from optparse import OptionParser
from readSampleInfo import *
from Colours import *


############################################################################
############################## User's inputs ###############################
############################################################################
regions = ['SR','CR2P2F','CR3P1F','CR3P1FXOR']

parser = OptionParser(usage="usage: %prog <analysis> <sample> [options]")
parser.add_option("-r", "--region", dest="region",
                  default="SR",
                  help="Region type are {0:s}. Default is SR.".format(', '.join(regions)))

parser.add_option("-e", "--external-cross-section", dest="getExternalCrossSectionFromFile",
                  action="store_true",
                  help="Use this option if you want to read the cross-section from the csv file")

parser.add_option("-l", "--luminosity", dest="luminosity",
                  type='int',
                  default=19712,
                  help="Set luminosity scenario from command line. Default is 19712/pb.")

parser.add_option("-d", "--directory", dest="directory",
                  default="samples",
                  help="Sample location, default is ./samples")



(options, args) = parser.parse_args()

analysis   = args[0]
typeofsample = args[1]
region = options.region

if region not in regions:
    print region, "is an unknown region. Run {0:s} -h for more details.".format(sys.argv[0])
    sys.exit(1)

getExternalCrossSectionFromFile = False if options.getExternalCrossSectionFromFile is None else options.getExternalCrossSectionFromFile

luminosity = options.luminosity
baseinputdir = options.directory
if not region == 'SR': baseinputdir = baseinputdir+"/"+region




###########################################################################
###########################################################################
print "\n\n"
print "\t\t\t",White('*** UNITO Framework ***')
print "\n\n"
###########################################################################
###########################################################################



############################## Configuration ##############################
################ Change this if you know what you are doing ###############
###########################################################################

executable = "eventAnalyzer" 
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

registry = './src/AnalysisFactory.cc'
checkIfAnalysisIsRegistered = 'grep -r {0:s} {1:s} | grep create | wc -l'.format(analysis, registry)
failure, output = commands.getstatusoutput(checkIfAnalysisIsRegistered)
if (not int(output) ==1):
    print Important("ERROR! The analysis {0:s} is not registered in {1:2}.".format(analysis,registry))
    howToRegister = Yellow('Register("{0:s}", &{0:s}::create);'.format(analysis))
    print "If you have not mispelled the name of your analysis, then please register it adding {0:s} in the constructor of {1:s} and recompile the code.".format(howToRegister, registry)
    sys.exit(1)
    

############################################################################
############################ Print the configuration #######################
############################################################################

print Blue("----------------------------------------------------------------------")
print Red("Configuration:")
print "Executable: {0:s} and analysis: {1:s}".format(Blue(executable), Blue(analysis)) 
print "Sample/type of samples:", Blue(typeofsample)
print "CSV file: ", Blue(csvfile)
print "Get (again) cross sections from csv file: ", Blue(getExternalCrossSectionFromFile)
print "Region type: ", Blue(region)
print "Integrated luminosity: ", Blue(luminosity)
print Blue("----------------------------------------------------------------------")
print "\n"


############################################################################


def run(executable, analysis, typeofsample, region, luminosity):
    inputdir  = baseinputdir
    outputdir = 'results'
    if not os.path.exists(outputdir): os.popen('mkdir "%s"' %outputdir)

    outputdir = outputdir+"/"+analysis+"_"+region
    if not os.path.exists(outputdir): os.popen('mkdir "%s"' %outputdir)


    datasets = getSamplesBy('process',typeofsample,csvfile)
    if len(datasets) == 0:
        datasets = getSamplesBy('identifier',typeofsample,csvfile)
    if len(datasets) == 0 and not typeofsample == 'test':
        print Important('Error! This sample is not available!'), typeofsample

    sampleprefix = ''

    

    ### Override configuration, for test only ###
    if typeofsample == 'test':
        datasets = ['test']
        

    #################################################################################

    ### Special treatment for DATA ###

    #################################################################################
    isData = False
    if typeofsample[0:8] == 'DoubleMu' or typeofsample[0:9] == 'DoubleEle' or typeofsample[0:4] == 'MuEG':
        luminosity = -1
        isData = True
        

    # ----- Run over the run periods -----
    hadd = 'hadd {0:s}/{1:s}.root'.format(outputdir,typeofsample)
    for period in datasets:
        basefile = sampleprefix+period
        if os.path.exists('{0:s}/{1:s}.root'.format(outputdir,basefile)):
            os.popen('rm {0:s}/{1:s}.root'.format(outputdir,basefile))

        externalXsec = -1
        if not isData and getExternalCrossSectionFromFile:
            externalXsec = crossSection(period, csvfile)
            print "For {0:s} {1:s} {2:.6f}".format(period, Warning("Using external cross section:"), externalXsec)

        print Red('\n------------------------------ {0:s} -------------------------------\n'.format(basefile))
        command = "./{0:s} {1:s} {2:s} {3:s}/{5:s}.root {4:s}/{5:s}.root {6:.0f} {7:.10f}".format(executable,analysis,region,inputdir,outputdir, basefile, luminosity, externalXsec)
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
    elif len(datasets) == 1 and not datasets[0] == typeofsample:
        print "One sample in the dataset, just copying it."
        os.popen('cp {0:s}/{1:s}.root {0:s}/{2:s}.root'.format(outputdir,datasets[0],typeofsample))

    print "The output is in", Green('{0:s}/{1:s}.root'.format(outputdir,typeofsample))  


if typeofsample == 'all' or typeofsample == 'data':
    for sample in typeofsamples:
        if typeofsample == 'all' or sample[0:8] == 'DoubleMu' or sample[0:9] == 'DoubleEle' or sample[0:5] == 'MuEG':
            if region == 'all':
                for cr in range(0,4):
                    run(executable, analysis, sample, cr, luminosity)    # runs over all samples in all control reagions
            else:
                run(executable, analysis, sample, region, luminosity)   # runs over all samples in a specific control reagions
else:
    if region == 'all':
        for cr in range(0,4):     
            run(executable, analysis, typeofsample, cr, luminosity)  # runs over a specific sample in all control regions
    else:
        run(executable, analysis, typeofsample, region, luminosity) # runs over a specific sample in a specific region

print "\nJob status: ", OK("DONE"),"\n"
