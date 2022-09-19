#!/usr/bin/env python2

#############################################################################################
# Get histograms from an analyzer's results (file: sample, inside has histograms: variable) #
# and write it in a form usable by Combine (file: variable, inside histograms: sample       #
# Wuthor A. Mecca                                                                           #
#############################################################################################

from __future__ import print_function
import os, sys
import ROOT
from plotUtils import TFileContext, makedirs_ok
from subprocess import call

# The configuration
blinded = False #True
inputDir = 'results'
outputDir = 'histogramsForCombine'
analyzer = 'VVGammaAnalyzer'
regions = ['SR4P', 'CR3P1F' , 'CR2P2F' , 'SR4P_1L', 'SR4P_1P', 'CR4P_1F', 'CR4L',    
           'SR3P', 'CR110'  , 'CR101'  , 'CR011'  , 'CR100'  , 'CR001'  , 'CR010', 'CR000', 'SR3P_1L', 'SR3P_1P', 'CR3P_1F', 'CRLFR', 'CR3L',
           'SR2P', 'SR2P_1L', 'SR2P_1P', 'CR2P_1F'
           # 'SR_HZZ', 'CR2P2F_HZZ', 'CR3P1F_HZZ', 'CR_HZZ', 'MC_HZZ',     
           # 'MC'
]
years = ['2016']

# Utility functions
def skipsample(filename):
    if '201' in filename: return True
    if filename == 'data.root': return True
    if(filename.startswith(('ggTo4e', 'ggTo2e2mu', 'ggTo4mu', 'ggZZ', 'ZZZ', 'WZZ', 'WWZ', 'TTWW', 'TTZZ'))):
        return True
    if filename.split('.')[-1] != 'root':
        print(">>> strange sample:", filename)
        return True
    return False

def isVarSystematic(variable):
    return variable[:3] == 'SYS'

# Setup
if(os.path.isdir(outputDir)):
    call(['rm', '-r', outputDir])  # clean up the directory so that old files don't clutter it
   
# samples = set()  # deduce from files found
n_retrieved = 0
not_retrieved = []

# Start

# Output nominal
# schema: <year>/<region>.root -> <variable>/<sample>
# example: 2016/SR4P.root      -> mZZ/ZZTo4l

# Output systematics
# schema: <year>/<region>.root -> <variable>/<sample>_CMS_<syst>(Up|Down)
# example: 2016/SR4P.root      -> mZZ/ZZTo4l_CMS_QCDScale-muRUp

for year in years:
    path_out = os.path.join(outputDir, year)
    makedirs_ok(path_out)
    
    for region in regions:
        path_in  = os.path.join(inputDir , year, analyzer+'_'+region)
        if(not os.path.isdir(path_in)):
            print("WARN: Skipping non-existent dir:", path_in)
            continue

        samples_region = set([ d.rstrip('.root') for d in os.listdir(path_in) if not skipsample(d) ])
        print('INFO: region={} samples:'.format(region), samples_region)
        # samples.update(samples_region)
        files_in = {}  # mapping: sample <str> --> file <ROOT.TFile>

        # retrieve the full list of variables in this region. Open input files but don't close them
        variables_region = set()
        for sample in samples_region:
            files_in[sample] = ROOT.TFile(os.path.join(path_in, sample+'.root'))
            variables_region.update([ key.GetName() for key in files_in[sample].GetListOfKeys() if key.ReadObj().Class().InheritsFrom("TH1")])

        # Add data
        if(not blinded):
            try:
                data_obs = ROOT.TFile(os.path.join(path_in, 'data.root'))
            except OSError as e:
                print('WARN: While opening {}, caught'.format(data_obs.GetName()), e)
            else:
                samples_region.add('data_obs')
                files_in['data_obs'] = data_obs

        # Write to output
        with TFileContext(os.path.join(path_out, region+'.root'), "RECREATE") as fout:
            for variable in variables_region:
                if(not isVarSystematic(variable)):
                   continue 
                split = variable.split('_')
                var_name = split[1]
                syst = split[2]
                if(syst == 'central'):
                    skipIfData = False
                    out_name = '{sample}'.format(sample='%s')
                else:
                    skipIfData = True
                    direction = split[3]
                    if(direction == 'Dn'): direction='Down'
                    out_name = '{sample}_CMS_{syst}{direction}'.format(sample='%s', syst=syst, direction=direction)

                subdir = fout.Get(var_name)  # e.g. mZZ, mZZG
                if(not subdir):
                    subdir = fout.mkdir(var_name)
                subdir.cd()  # cd into this TDirectory
                
                for sample in samples_region:
                    if(sample == 'data_obs' and skipIfData):
                        continue
                    h = files_in[sample].Get(variable)
                    if(h):
                        h.SetName(out_name %(sample))
                        h.Write()
                        n_retrieved += 1
                    else:
                        not_retrieved.append([files_in[sample].GetName(), variable])

        del samples_region
        for _, handler in files_in.items():
            handler.Close()


for e in not_retrieved[:10]:
    print('ERROR: Could not get {:32s} from file {}'.format(e[1], e[0]))
if(len(not_retrieved) > 10):
    print('    + {:d} more'.format(len(not_retrieved) - 10))

print('INFO: Retrieved and wrote {:d} histograms'.format(n_retrieved))
