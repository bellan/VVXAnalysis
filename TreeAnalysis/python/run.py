#! /usr/bin/env python2.7

##################################
## R. Bellan (UNITO) - Feb 2014 ##
##################################

import sys, os, commands, math, subprocess
from optparse import OptionParser
from readSampleInfo import *
from Colours import *


############################################################################
############################## User's inputs ###############################
############################################################################
regions = ['SR','CR','CR2P2F','CR3P1F','SR_HZZ','CR_HZZ','CR2P2F_HZZ','CR3P1F_HZZ', 'MC','MC_HZZ']
years   = [2016,2017,2018]

parser = OptionParser(usage="usage: %prog <analysis> <sample> [options]")
parser.add_option("-r", "--region", dest="region",
                  default="SR",
                  help="Region type are {0:s}. Default is SR.".format(', '.join(regions)))

parser.add_option("-e", "--external-cross-section", dest="getExternalCrossSectionFromFile",
                  action="store_true",
                  help="Use this option if you want to force to read again the cross-section from the csv file")

parser.add_option("-y", "--year", dest="year",
                  type='int',
                  default= 2016,
                  help="Set year scenario from command line. Default is 2016.")

parser.add_option("-l", "--luminosity", dest="luminosity",
                  type='int',
                  default= None,
                  help="Set luminosity scenario from command line. The default is the one taken from the target year. Default year is 2016 and L = 35900/pb. Full lumi is 137100/pb.")

parser.add_option("-d", "--directory", dest="directory",
                  default="samples",
                  help="Sample location, default is ./samples/2016")

parser.add_option("-c", "--csv", dest="csvfile",
                  #default="../Producers/python/samples_2016_MC.csv",
                  default=None,
                  help="csv path, default is ../Producers/python/samples_<year>_MC.csv")

parser.add_option("-s", "--scalefactor", dest="doSF",
                  action="store_true",
                  default=False,
                  help="take SF from tree-analyzer instead of entuple internal values")

parser.add_option("-n", "--nevents", dest="maxNumEvents",
                  type='int',
                  default= -1,
                  help="Set max number of events to run over. Default is -1, meaning all events in the tree")



(options, args) = parser.parse_args()

analysis     = args[0]
typeofsample = args[1]
region       = options.region
doSF         = options.doSF
if doSF:
    print "Option temporarily disabled!"
    sys.exit(1)
maxNumEvents = options.maxNumEvents
year         = options.year
luminosity   = options.luminosity

if region not in regions:
    print region, "is an unknown region. Run {0:s} -h for more details.".format(sys.argv[0])
    sys.exit(1)

getExternalCrossSectionFromFile = False if options.getExternalCrossSectionFromFile is None else options.getExternalCrossSectionFromFile

csvfile = options.csvfile

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
print "Get (again) cross sections from csv file: ", Blue(getExternalCrossSectionFromFile)
print "Region type: ", Blue(region)
print "Use internal scale factor",Blue(doSF)
#print "Integrated luminosity: ", Blue(luminosity)



############################################################################


def run(executable, analysis, typeofsample, region, year, luminosity, maxNumEvents, doSF):

    # if luminosity is specified thorugh -l option, overwrite the year <-> luminosity decision
    if luminosity is None:
        if year == 2016:
            luminosity =  35900
        elif year == 2017: 
            luminosity =  41500
        elif year == 2018: 
            luminosity =  59700
        else :
            print"{0:s}: Unknown year, please specify a luminosity with -l option".format(year)
            sys.exit(1)

    #################################################################################
    ### Special treatment for DATA ###
    #################################################################################
    isData = False

    if typeofsample[0:8] == 'DoubleMu' or typeofsample[0:9] == 'DoubleEle' or typeofsample[0:4] == 'MuEG' or typeofsample[0:6] == 'Single' or typeofsample[0:4] == 'test' or  typeofsample[0:6] == 'MuonEG' or  typeofsample[0:6] == 'MuonEG' or  typeofsample[0:8] == 'DoubleEG':
        luminosity = -1
        isData = True

    print "Year: ", Blue(year)
    print "Integrated luminosity: ", Blue(luminosity)
    print Blue("----------------------------------------------------------------------")
    print "\n"

    
    inputdir = baseinputdir
    
    outputdir = 'results'
    if not os.path.exists(outputdir): os.popen('mkdir "%s"' %outputdir)

    outputdir = outputdir+"/"+str(year)
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
        command = "./{0:s} {1:s} {2:s} {3:s}/{5:s}.root {4:s}/{5:s}.root {6:.0f} {7:.0f} {8:.5f} {9:.0f} {10:b}".format(executable,analysis,region,inputdir,outputdir, basefile, year, luminosity, externalXsec, maxNumEvents, doSF)
        print "Command going to be executed:", Violet(command)

        output = subprocess.call(command,shell=True)
        print "\n",output
        hadd = '{0:s} {1:s}/{2:s}.root'.format(hadd, outputdir, basefile)

    print Red('----------------------------------------------------------------------\n')
    if len(datasets) > 1:
        if os.path.exists('{0:s}/{1:s}.root'.format(outputdir,typeofsample)):
            os.popen('rm {0:s}/{1:s}.root'.format(outputdir,typeofsample))
        print "Command going to be executed:", Violet(hadd)
        output = subprocess.call(hadd,shell=True)
    elif len(datasets) == 1 and not datasets[0] == typeofsample:
        print "One sample in the dataset, just copying it."
        os.popen('cp {0:s}/{1:s}.root {0:s}/{2:s}.root'.format(outputdir,datasets[0],typeofsample))

    output = '{0:s}/{1:s}.root'.format(outputdir,typeofsample)
    print "The output is in", Green(output)
    return output



def mergeDataSamples(outputLocations):
    print outputLocations
    failure, basename = commands.getstatusoutput('basename {0:s}'.format(outputLocations[0]))
    outputdir = outputLocations[0].replace(basename,'')
    hadd = 'hadd {0:s}/data.root '.format(outputdir)
    for result in outputLocations:
        hadd = '{0:s} {1:s}'.format(hadd, result)
    if os.path.exists('{0:s}/data.root'.format(outputdir)):
        os.popen('rm {0:s}/data.root'.format(outputdir))
    print "Command going to be executed:", Violet(hadd)
    output = subprocess.call(hadd.split(),shell=True)

def runOverCRs(executable, analysis, sample, year, luminosity, maxNumEvents, doSF, postfix = '', outputLocations = []):
    outputCR2P2F = run(executable, analysis, sample, 'CR2P2F'+postfix, year, luminosity, maxNumEvents, doSF)    # runs over all samples in the CR2P2F control reagion
    outputCR3P1F = run(executable, analysis, sample, 'CR3P1F'+postfix, year, luminosity, maxNumEvents, doSF)    # runs over all samples in the CR3P1F control reagion

    if not os.path.exists('results/{0:s}_CR{1:s}'.format(analysis,postfix)): os.popen('mkdir results/{0:s}_CR{1:s}'.format(analysis,postfix))
    outputRedBkg = 'results/{0:s}_CR{1:s}/reducible_background_from_{2:s}.root'.format(analysis, postfix, sample)
    hadd = 'hadd {0:s} {1:s} {2:s}'.format(outputRedBkg, outputCR2P2F, outputCR3P1F)
    if os.path.exists('{0:s}'.format(outputRedBkg)):
        os.popen('rm {0:s}'.format(outputRedBkg))
    print "Command going to be executed:", Violet(hadd)
    output = subprocess.call(hadd.split(),shell=True)
    outputLocations.append(outputRedBkg)

    
def runOverSamples(executable, analysis, typeofsample, region, year, luminosity, maxNumEvents, doSF):

    if typeofsample == 'all' or typeofsample == 'data':
        outputLocations = []
        for sample in typeofsamples:
            if typeofsample == 'all' or sample[0:8] == 'DoubleMu' or sample[0:9] == 'DoubleEle' or sample[0:4] == 'MuEG' or sample[0:9]== "SingleEle" or sample[0:8]== "SingleMu":

                if region == 'all':
                    for cr in regions:
                        run(executable, analysis, sample, cr, year, luminosity, maxNumEvents, doSF)    # runs over all samples in all control reagions
                elif region == 'CR':
                    runOverCRs(executable, analysis, sample, year, luminosity, maxNumEvents, doSF, "",outputLocations)
                elif region == 'CR_HZZ': 
                    runOverCRs(executable, analysis, sample, year, luminosity, maxNumEvents, doSF, '_HZZ',outputLocations)
                else:
                    outputLocations.append(run(executable, analysis, sample, region, year, luminosity, maxNumEvents, doSF))   # runs over all samples in a specific control reagions
        if typeofsample == 'data':
            mergeDataSamples(outputLocations)
            
    else:
        if region == 'all':
            for cr in range(0,4):     
                run(executable, analysis, typeofsample, cr, year, luminosity, maxNumEvents, doSF)  # runs over a specific sample in all control regions

        elif region == 'CR':
            runOverCRs(executable, analysis, typeofsample, year, luminosity, maxNumEvents, doSF)
        elif region == 'CR_HZZ':
            runOverCRs(executable, analysis, typeofsample, year, luminosity, maxNumEvents, doSF, postfix='_HZZ')
        else:
            run(executable, analysis, typeofsample, region, year, luminosity, maxNumEvents, doSF) # runs over a specific sample in a specific region


if year == 1618:
    for year in years:
        if options.csvfile is None:
            csvfile = "../Producers/python/samples_"+str(year)+"_MC.csv"

        print "CSV file: ", Blue(csvfile)
        typeofsamples = typeOfSamples(csvfile)
        typeofsamples.append('test')
        baseinputdir = options.directory+'/'+str(year)
        
        runOverSamples(executable, analysis, typeofsample, region, year, luminosity, maxNumEvents, doSF)
elif year in years:
    if options.csvfile is None:
        csvfile = "../Producers/python/samples_"+str(year)+"_MC.csv"
    print "CSV file: ", Blue(csvfile)
    typeofsamples = typeOfSamples(csvfile)
    typeofsamples.append('test')
    baseinputdir = options.directory+'/'+str(year)
    
    runOverSamples(executable, analysis, typeofsample, region, year, luminosity, maxNumEvents, doSF)
else:
    print "Unknown year"
    sys.exit(1)    

        
print "\nJob status: ", OK("DONE"),"\n"
