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
from argparse import ArgumentParser

# The configuration
regions = ['SR4P', 'CR3P1F' , 'CR2P2F' , 'SR4P_1L', 'SR4P_1P', 'CR4P_1F', 'CR4L',    
           'SR3P', 'CR110'  , 'CR101'  , 'CR011'  , 'CR100'  , 'CR001'  , 'CR010', 'CR000', 'SR3P_1L', 'SR3P_1P', 'CR3P_1F', 'CRLFR', 'CR3L',
           'SR2P', 'SR2P_1L', 'SR2P_1P', 'CR2P_1F'
           # 'SR_HZZ', 'CR2P2F_HZZ', 'CR3P1F_HZZ', 'CR_HZZ', 'MC_HZZ',     
           # 'MC'
]

parser = ArgumentParser()
parser.add_argument('-y', '--year'     , default='2016')
parser.add_argument(      '--blind'    , dest='unblind', action='store_false', default=True)
parser.add_argument('-i', '--inputdir' , default='results')
parser.add_argument('-o', '--outputdir', default='histogramsForCombine')
parser.add_argument('-A', '--analyzer' , default='VVGammaAnalyzer')
parser.add_argument('-r', '--regions'  , default=regions, nargs='+', choices=regions)
parser.add_argument('-v', '--verbose'  , dest='verbosity', action='count', default=1, help='Increase the verbosity level')
parser.add_argument('-q', '--quiet'    , dest='verbosity', action='store_const', const=0)
args = parser.parse_args()

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

# Output nominal
# schema: <year>/<region>.root -> <variable>/<sample>
# example: 2016/SR4P.root      -> mZZ/ZZTo4l

# Output systematics
# schema: <year>/<region>.root -> <variable>/<sample>_CMS_<syst>(Up|Down)
# example: 2016/SR4P.root      -> mZZ/ZZTo4l_CMS_QCDScale-muRUp

if __name__ == '__main__':
    # Setup
    if(os.path.isdir(args.outputdir)):
        call(['rm', '-r', args.outputdir])  # clean up the directory so that old files don't clutter it

    # samples = set()  # deduce from files found
    n_retrieved = 0
    not_retrieved = []

    # Start

    path_out = os.path.join(args.outputdir, args.year)
    makedirs_ok(path_out)
    
    for region in args.regions:
        path_in  = os.path.join(args.inputdir , args.year, args.analyzer+'_'+region)
        if(not os.path.isdir(path_in)):
            print("WARN: Skipping non-existent dir:", path_in)
            continue

        samples_region = set([ d.rstrip('.root') for d in os.listdir(path_in) if not skipsample(d) ])
        print('INFO: region={} samples:'.format(region), samples_region)

        files_in = {}  # mapping: sample <str> --> file <ROOT.TFile>

        # retrieve the full list of variables in this region. Open input files but don't close them
        variables_region = set()
        for sample in samples_region:
            files_in[sample] = ROOT.TFile(os.path.join(path_in, sample+'.root'))
            variables_region.update([ key.GetName() for key in files_in[sample].GetListOfKeys() if key.ReadObj().Class().InheritsFrom("TH1") and isVarSystematic(key.GetName())])

        # print('DEBUG: in', region, 'the variables are:', variables_region)
        # Add data
        if(args.unblind):
            try:
                data_obs = ROOT.TFile(os.path.join(path_in, 'data.root'))
            except OSError as e:
                print('WARN: While opening {}, caught'.format(data_obs.GetName()), e)
            else:
                files_in['data_obs'] = data_obs

        # Write fake_photons
        if(not 'fake_photons' in files_in):
            with TFileContext(os.path.join(path_in, 'fake_leptons.root'), 'RECREATE') as fFakePh:
                for variable in variables_region:
                    split = variable.split('_')
                    var_name = split[1]
                    syst = split[2]
                    if(var_name.endswith('failReweight')):
                        h = files_in['data_obs'].Get(variable)
                        new_name = var_name.replace('failReweight', 'loose')
                        h.SetName( h.GetName().replace(var_name, new_name) )
                        h.Write()
            try:
                fFakePh = ROOT.TFile(os.path.join(path_in, 'data.root'))
            except OSError as e:
                print('WARN: While opening {}, caught'.format(fFakePh.GetName()), e)
            else:
                files_in['fake_photons'] = fFakePh

        # Write to output
        with TFileContext(os.path.join(path_out, region+'.root'), "RECREATE") as fout:
            for variable in variables_region:
                split = variable.split('_')
                var_name = split[1]
                syst = split[2]
                if(var_name.endswith('failReweight')):
                    continue
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
                
                for sample, file_in in files_in.items():
                    if(sample == 'data_obs' and skipIfData):
                        continue
                    h = file_in.Get(variable)
                    if(h):
                        h.SetName(out_name %(sample))
                        h.Write()
                        n_retrieved += 1
                    else:
                        not_retrieved.append([files_in[sample].GetName(), variable])

        del samples_region
        for _, handler in files_in.items():
            handler.Close()

if(args.verbosity >= 1):
    files_prob = set([e[0] for e in not_retrieved])
    for file_prob in files_prob:
        problems = [e[1] for e in not_retrieved if e[0] == file_prob]
        print('WARN: From file {:60.60s} could not retrieve {:2d} plots'.format(file_prob, len(problems)))

if(args.verbosity >= 2):
    hists_prob         = { e[1] for e in not_retrieved }
    hists_prob_central = { e for e in hists_prob if e.endswith('central') }
    hists_prob_updn    = { e.rstrip('_Up').rstrip('_Down') for e in hists_prob if e.endswith(('Up', 'Down')) }
    for hist_prob in hists_prob_central:
        problems = { e[0] for e in not_retrieved if e[1] == hist_prob }
        print('WARN: Histogram {:40.40s} was missing from {:2d} files'.format(hist_prob, len(problems)), end='')
        if(len(problems) < 10):
            print(':', *[f.split('/')[-1] for f in problems])
        else:
            print()
    for hist_prob in hists_prob_updn:
        problems = { e[0] for e in not_retrieved if e[1].startswith(hist_prob) }
        print('WARN: Histogram {:40.40s} was missing from {:2d} files'.format(hist_prob+'(Up/Down)', len(problems)), end='')
        if(len(problems) < 10):
            print(':', *[f.split('/')[-1] for f in problems])
        else:
            print()

print('INFO: Total missing plots:', len(not_retrieved))
print('INFO: Retrieved and wrote {:d} histograms'.format(n_retrieved))
