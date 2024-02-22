#!/usr/bin/env python

##################################
## R. Bellan (UNITO) - Feb 2014 ##
##################################

from __future__ import print_function
import sys, os, math, subprocess
from optparse import OptionParser
from readSampleInfo import *
from Colours import *
from utils23 import lumi_dict

if(sys.version_info.major >= 3):
    popen_extra_args = {'encoding': 'utf-8'}
else:
    popen_extra_args = {}

############################################################################
############################## User's inputs ###############################
############################################################################
# FIXME: change name
allowed_regions = ['SR4P', 'CR3P1F' , 'CR2P2F' , 'SR4P_1L', 'SR4P_1P', 'CR4P_1F', 'CR4L',
                   'SR3P', 'CR110'  , 'CR101'  , 'CR011'  , 'CR100'  , 'CR001'  , 'CR010', 'CR000', 'SR3P_1L', 'SR3P_1P', 'CR3P_1F', 'CRLFR', 'CR3L',
                   'SR2P', 'SR2P_1L', 'SR2P_1P', 'CR2P_1F',
                   'SR_HZZ', 'CR2P2F_HZZ', 'CR3P1F_HZZ', 'CR_HZZ', 'MC_HZZ',
                   'MC']

_all_regions = ['SR4P', 'CR3P1F', 'CR2P2F', 'SR3P', 'CR110', 'CR101', 'CR011', 'CR100', 'CR001', 'CR010', 'CR000', 'SR2P']  # Regions tu run on when "all" are requested

CR4L_regions = ['CR3P1F', 'CR2P2F']

CR3L_regions = ['CR110', 'CR101', 'CR011', 'CR100', 'CR001', 'CR010', 'CR000']

CR_HZZ_regions = ['CR3P1F_HZZ', 'CR2P2F_HZZ']

years   = ['2016preVFP', '2016postVFP', '2016', '2017', '2018']

parser = OptionParser(usage="usage: %prog <analysis> <sample> [options]")

parser.add_option("-r", "--regions", dest="regions",
                  default="SR4P",
                  help="Region type are {0:s}. Default is %default.".format(', '.join(allowed_regions)))


parser.add_option("-e", "--external-cross-section", dest="getExternalCrossSectionFromFile",
                  action="store_true",
                  default = False,
                  help="Use this option if you want to force to read again the cross-section from the csv file")

parser.add_option("-i", "--internal-cross-section", dest="useInternalCrossSection",
                  action="store_true",
                  default = False,
                  help="Use this option if you want to force to use the internal cross section of the original sample")

parser.add_option("-y", "--year", dest="year",
                  type=str,
                  default='2018',
                  help="Set year scenario from command line (default: %default)")

parser.add_option("-l", "--luminosity", dest="luminosity",
                  type='int',
                  default= None,
                  help="Set luminosity scenario from command line. The default is the one taken from the target year. Full lumi for Run2 is 137100/pb.")

parser.add_option("-d", "--directory", dest="directory",
                  default="samples",
                  help="Sample location (default: %default)")

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

parser.add_option("-o", "--outputdir",
                  default="results",
                  help="Directory in which the output will be stored (default: %default)")

parser.add_option("-u", "--unblind", dest="unblind",
                  action="store_true",
                  default=False,
                  help="unblind the search regions. This option is active only on data.")

parser.add_option("--nofr", dest="nofr",
                  action="store_true",
                  default=False,
                  help="do not apply the fake rate scale factors.")

parser.add_option("--fpw","--force-pos-weight", dest="forcePosWeight",
                  action="store_true",
                  default=False,
                  help="Force the weight to be positive.")


(options, args) = parser.parse_args()

analysis       = args[0]
typeofsample   = args[1]
regions        = options.regions
doSF           = options.doSF
unblind        = options.unblind
nofr           = options.nofr
forcePosWeight = options.forcePosWeight

#samples = list(set(typeofsample.strip(';,').split( (';' if ';' in typeofsample else ',') )))

if doSF:
    print("Option temporarily disabled!")
    sys.exit(1)
maxNumEvents = options.maxNumEvents

year         = options.year
if year not in years:
    print("Unknown year", year)
    sys.exit(1)

luminosity   = options.luminosity

regions = list(set(regions.strip(';,').split( (';' if ';' in regions else ',') )))
    
for region in regions:  
    if region not in allowed_regions and not region == 'all':
        print(region, "is an unknown region. Run {0:s} -h for more details.".format(sys.argv[0]))
        print("Multiple regions must be separated by a ';'")
        sys.exit(1)

if len(regions) == 1:
    region = regions[0]

# Some more logic for regions
tmp = list(regions)  # temporary copy
if 'all' in regions:
    regions = _all_regions
elif 'MC' in regions:
    regions = ['MC']
elif 'CR4L' in regions:
    regions = CR4L_regions
elif 'CR3L' in regions:
    regions = CR3L_regions
else:
    if 'SR4P' in regions:
        regions = [ r for r in regions if not r.startswith(('SR4P_', 'CR4P_')) ]
    if 'SR3P' in regions:
        regions = [ r for r in regions if not r.startswith(('SR3P_', 'CR3P_')) ]
    if 'SR2P' in regions:
        regions = [ r for r in regions if not r.startswith(('SR2P_', 'CR2P_')) ]
        
if(len(regions) != len(tmp)):
    print("WARN: some regions have been dropped to avoid overlap:", [r for r in tmp if r not in regions])
    print("Remaining regions:", regions)

isData = False
if typeofsample.startswith( ('DoubleMu', 'DoubleEle', 'EGamma', 'MuEG', 'Single', 'MuonEG', 'DoubleEG', str(year), 'test', 'data') ):  # giving a tuple of prefixes returns true if at least one matches
    luminosity = -1
    isData = True
    
csvfile = options.csvfile

###########################################################################
###########################################################################
print("\n\n")
print("\t\t\t",White('*** UNITO Framework ***'))
print("\n\n")
###########################################################################
###########################################################################



############################## Configuration ##############################
################ Change this if you know what you are doing ###############
###########################################################################

executable = "eventAnalyzer" 


output = subprocess.check_output('ls ./bin/ | grep -v .cpp | grep -v .xml | grep -v .md', shell=True, **popen_extra_args).rstrip('\n')
availableExecutable = output.split()
if executable in availableExecutable: 
    executable = 'bin/'+executable
else:
    print(Important("ERROR! Unknown executable."),"Availble executables are:",availableExecutable)
    sys.exit(1)

registry = './src/AnalysisFactory.cc'
checkIfAnalysisIsRegistered = 'grep -r {0:s} {1:s} | grep create | wc -l'.format(analysis, registry)
output = subprocess.check_output(checkIfAnalysisIsRegistered, shell=True, **popen_extra_args).rstrip('\n')
if (not int(output) ==1):
    print(Important("ERROR! The analysis {0:s} is not registered in {1:2}.".format(analysis,registry)))
    howToRegister = Yellow('Register("{0:s}", &{0:s}::create);'.format(analysis))
    print("If you have not mispelled the name of your analysis, then please register it adding {0:s} in the constructor of {1:s} and recompile the code.".format(howToRegister, registry))
    sys.exit(1)
    


    
############################################################################
############################ Print the configuration #######################
############################################################################

print(Blue("----------------------------------------------------------------------"))
print(Red("Configuration:"))
print("Executable: {0:s} and analysis: {1:s}".format(Blue(executable), Blue(analysis)))
if isData:
    print("Running on", Blue("data"))
else:
    print("Running on", Blue("MC"))
print("Sample/type of samples:", Blue(typeofsample))
print("Get (again) cross sections from csv file: ", Blue(options.getExternalCrossSectionFromFile))
print("Use internal cross section from sample: ", Blue(options.useInternalCrossSection), '[TODO]')
print("Region type: ", Blue(regions))
print("Use internal scale factor: ",Blue(doSF))

print(region)

############################################################################


def run(executable, analysis, typeofsample, regions, year, luminosity, maxNumEvents, doSF, unblind, nofr, forcePosWeight):

    # if luminosity is specified thorugh -l option, overwrite the year <-> luminosity decision
    if luminosity is None:
        luminosity = lumi_dict[year]['value']

    #################################################################################
    ### Blind DATA by default ###
    #################################################################################

    if isData is False:
        unblind = True
        
    print("Year: ", Blue(year))
    print("Integrated luminosity: ", Blue(luminosity))
    if isData: 
        if unblind:
            print("Running mode:", Evidence("unblinded!"))
        else:
            print("Running mode:", Blue("blinded"))

    print(Blue("----------------------------------------------------------------------"))
    print("\n")

    
    inputdir = options.directory+'/'+str(year)
    
    outputdir_format = os.path.join(options.outputdir, str(year), analysis+'_%s')
    outputdirs = {}
    for region in regions:
        odir = outputdir_format %region
        outputdirs[region] = odir
        os.popen('mkdir -p "%s"' %odir)


    datasets = getSamplesBy('process',typeofsample,csvfile)
    if len(datasets) == 0:
        datasets = getSamplesBy('identifier',typeofsample,csvfile)
    if len(datasets) == 0 and not typeofsample == 'test':
        raise ValueError('The sample {0:s} is not available in the CSV!'.format(typeofsample))
    
    sampleprefix = ''
  

    ### Override configuration, for test only ###
    if typeofsample == 'test':
        datasets = ['test']


    
    # ----- Run over the run periods -----
    hadd_cmds = {}
    for region, odir in outputdirs.items():
      hadd_cmds[region] = 'hadd -k -f {0:s}/{1:s}.root'.format(odir,typeofsample)

    for period in datasets:
        basefile = sampleprefix+period

        externalXsec = -1
        if not isData and options.getExternalCrossSectionFromFile:
            externalXsec = crossSection(period, csvfile)
            print("For {0:s} {1:s} {2:.6f}".format(period, Warn("Using external cross section:"), externalXsec))
        elif not isData and options.useInternalCrossSection:
            externalXsec = -10  # This is a hack: values < -2 are used to signal internal cross section
            print("For {0:s} {1:s}".format(period, Warn("Using internal cross section")))

        if not os.path.exists('{0:s}/{1}.root'.format(inputdir,basefile)):
            print(Warn("The ROOT file for the sample {0:s} does not exist in {1:s}".format(basefile,inputdir)))
            return 
        print(Red('\n------------------------------ {0:s} -------------------------------\n'.format(basefile)))
          
        # outputdir_format is something like "results/2016/VVXAnalyzer_%s". EventAnalyzer will use Form() to replace %s with the various regions to obtain the filenames
        keywords = locals()
        keywords.update({'regions':';'.join(regions)})
        command = "./{executable} {analysis} '{regions}' {inputdir}/{basefile}.root {outputdir_format}/{basefile}.root {year} {luminosity:.0f} {externalXsec:.5f} {maxNumEvents:d} {doSF:b} {unblind:b} {nofr:b} {forcePosWeight:b}".format(**keywords)
        print("Command going to be executed (run::command):", Violet(command))
        subprocess.check_call(command, shell=True)

        for r,outputdir in outputdirs.items():
            filepath = '{0:s}/{1:s}.root'.format(outputdir,basefile)
            if os.path.exists(filepath):
                hadd_cmds[r] += ' {0:s}'.format(filepath)
                

    print(Red('----------------------------------------------------------------------\n'))
    for region,outputdir in outputdirs.items():
        if len(datasets) > 1:
            if os.path.exists('{0:s}/{1:s}.root'.format(outputdir,typeofsample)):
                os.popen('rm {0:s}/{1:s}.root'.format(outputdir,typeofsample))
            print("Command going to be executed (run::hadd):", Violet(hadd_cmds[region]))
            subprocess.check_call(hadd_cmds[region], shell=True)
        elif len(datasets) == 1 and not datasets[0] == typeofsample:
            if os.path.exists('{0:s}/{1:s}.root {0:s}/{2:s}.root'.format(outputdir,datasets[0],typeofsample)):
                print("One sample in the dataset, just copying it.")
                os.popen('cp {0:s}/{1:s}.root {0:s}/{2:s}.root'.format(outputdir,datasets[0],typeofsample))

    output = {reg:'{0:s}/{1:s}.root'.format(odir,typeofsample) for reg,odir in outputdirs.items()}
    print("The output is in", Green([v for k,v in output.items()]))
    return output

###------------------------------------------------------------------------------------###

def mergeDataSamples(outputLocationsDict):

    if(len(outputLocationsDict) == 0):
        print(Red("Error") + ": outputLocations is empty!")
        exit(1)
        
    for region,outputLocations in outputLocationsDict.items():
        basename = subprocess.check_output(['basename', outputLocations[0]], **popen_extra_args).rstrip('\n')
        outputdir = outputLocations[0].replace(basename,'').rstrip('/')
        hadd = 'hadd -k -v 0 {0:s}/data.root {1:s}'.format(outputdir, ' '.join(outputLocations))
        
        if os.path.exists('{0:s}/data.root'.format(outputdir)):
            os.popen('rm {0:s}/data.root'.format(outputdir))
            
        print("Command going to be executed (mergeDataSamples::hadd):", Violet(hadd))
        subprocess.check_call(hadd, shell=True, **popen_extra_args)

###------------------------------------------------------------------------------------###
        
def mergeCRs(analysis, year, inputLocs, antype):
    
    inputLocations = {}

    if antype == 'CR4L':
        inputLocations = {k: inputLocs[k] for k in CR4L_regions}
    if antype == 'CR3L':
        inputLocations = {k: inputLocs[k] for k in CR3L_regions}
    if antype == 'CR_HZZ':
        inputLocations = {k: inputLocs[k] for k in CR_HZZ_regions}

        
    outdir = '{0:s}/{1:s}/{2:s}_{3:s}'.format(options.outputdir, str(year), analysis, antype)
    if not os.path.exists(outdir): os.popen('mkdir {0:s}'.format(outdir))
    if isData:
        outputRedBkg = '{0:s}/reducible_background.root'.format(outdir)
    else:
        outputRedBkg = '{0:s}/reducible_background_MC.root'.format(outdir)
    hadd = 'hadd -k -v 0 {0:s}'.format(outputRedBkg)
    for key, values in inputLocations.items():
        for value in values:
            hadd += ' {0:s}'.format(value)
    
    if os.path.exists('{0:s}'.format(outputRedBkg)):
        os.popen('rm {0:s}'.format(outputRedBkg))
        
    print("Command going to be executed (runOverCRs::hadd):", Violet(hadd))
    subprocess.check_call(hadd, shell=True)
    
    if isData:
        subprocess.check_call('ln -sf reducible_background.root {0:s}/data.root'.format(outdir),shell=True)

###------------------------------------------------------------------------------------###        
    
def runOverSamples(executable, analysis, typeofsample, regions, year, luminosity, maxNumEvents, knownProcesses, doSF, unblind, nofr, forcePosWeight):

    outputLocations = {}
    
    if typeofsample == 'all' or typeofsample == 'data':
        for sample in knownProcesses:
            if typeofsample == 'all' or sample[0:4] == str(year):
                if region == 'all':                
                    outputLocs = run(executable, analysis, sample, _all_regions, year, luminosity, maxNumEvents, doSF, unblind, nofr, forcePosWeight) # runs over all samples in all regions
                    if(outputLocs is None):
                        print("No output produced")
                        return 1
                    for r,loc in outputLocs.items():
                        outputLocations.setdefault(r, []).append(loc)

                else:
                    outputLocs = run(executable, analysis, sample, regions, year, luminosity, maxNumEvents, doSF, unblind, nofr, forcePosWeight)      # runs over all samples in specific regions
                    if(outputLocs is None):
                        print("No output produced")
                        return 1
                    for r,loc in outputLocs.items():
                        outputLocations.setdefault(r, []).append(loc)

        if typeofsample == 'data':
            mergeDataSamples(outputLocations)

        if region == 'CR4L' or region == 'CR_HZZ' or region == 'CR3L':
            mergeCRs(analysis, year, outputLocations, region)

        if region == 'all':
            mergeCRs(analysis, year, outputLocations, 'CR4L')
            mergeCRs(analysis, year, outputLocations, 'CR3L')

    else:
        if region == 'all':
            outputLocs = run(executable, analysis, typeofsample, _all_regions, year, luminosity, maxNumEvents, doSF, unblind, nofr, forcePosWeight)   # runs over a specific sample in all signal/control regions
            if(outputLocs is None):
                print("No output produced")
                return 1
            for r,loc in outputLocs.items():
                outputLocations.setdefault(r, []).append(loc)

            mergeCRs(analysis, year, outputLocations, 'CR4L')
            mergeCRs(analysis, year, outputLocations, 'CR3L')
           
        else:
            outputLocs = run(executable, analysis, typeofsample, regions, year, luminosity, maxNumEvents, doSF, unblind, nofr, forcePosWeight)        # runs over a specific sample in specific regions
            if(outputLocs is None):
                print("No output produced")
                return 1
            for r,loc in outputLocs.items():
                outputLocations.setdefault(r, []).append(loc)

            if region == 'CR4L' or region == 'CR_HZZ' or region == 'CR3L':
                mergeCRs(analysis, year, outputLocations, region)
    return 0

###------------------------------------------------------------------------------------###        
                
###################################
### Actual steering of the code ###
###################################

runStatus = 0
if year == 1618:
    for year in years:
        if options.csvfile is None:
            if isData:
                csvfile = "../Producers/python/samples_"+str(year)+"UL_Data.csv"
            else:
                csvfile = "../Producers/python/samples_"+str(year)+"UL_MC.csv"
                
        print("CSV file: ", Blue(csvfile))

        knownProcesses = typeOfSamples(csvfile)
        # knownProcesses.append('test')
        
        runStatus = runOverSamples(executable, analysis, typeofsample, regions, year, luminosity, maxNumEvents, knownProcesses, doSF, unblind, nofr, forcePosWeight)

elif year in years:
    if options.csvfile is None:
        if isData:
            csvfile = "../Producers/python/samples_"+str(year)+"UL_Data.csv"
        else:
            csvfile = "../Producers/python/samples_"+str(year)+"UL_MC.csv"
        
    print("CSV file: ", Blue(csvfile))

    knownProcesses = typeOfSamples(csvfile)
    # knownProcesses.append('test')

    runStatus = runOverSamples(executable, analysis, typeofsample, regions, year, luminosity, maxNumEvents, knownProcesses, doSF, unblind, nofr, forcePosWeight)

else:
    print("Unknown year", year)
    sys.exit(1)    
       
print("\nJob status: ", OK("DONE") if runStatus == 0 else Red("FAILED"),"\n")
