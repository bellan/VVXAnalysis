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
import logging


# Utility functions
def skipsample(filename):
    if '201' in filename: return True
    if filename == 'data.root': return True
    if filename in ('ggTo4l', 'ttXY', 'triboson'): return True
    if filename.split('.')[-1] != 'root':
        logging.error("strange sample:", filename)
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

def main():
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
    parser.add_argument('-r', '--regions'  , default=['SR4P'], nargs='+', choices=regions)
    parser.add_argument(      '--remake-fake-photons', action='store_true', help='Force to recreate the fake_photons.root file from data.root')
    parser.add_argument('-v', '--verbose'  , dest='verbosity', action='count', default=1, help='Increase the verbosity level')
    parser.add_argument('-q', '--quiet'    , dest='verbosity', action='store_const', const=0)
    parser.add_argument('--log', dest='loglevel', type=str.upper, metavar='LEVEL', default='WARNING')
    args = parser.parse_args()

    if(not hasattr(logging, args.loglevel)):
        raise ValueError('Invalid log level: %s' % args.loglevel)
    logging.basicConfig(format='%(levelname)s:%(module)s:%(funcName)s: %(message)s', level=getattr(logging, args.loglevel))

    # Setup
    ok_retrieved  = []
    not_retrieved = []

    # Start
    path_out = os.path.join(args.outputdir, args.year)
    makedirs_ok(path_out)
    
    for region in args.regions:
        path_in  = os.path.join(args.inputdir , args.year, args.analyzer+'_'+region)
        if(not os.path.isdir(path_in)):
            logging.warning('Skipping non-existent dir: '+path_in)
            continue

        samples_region = set([ d.rstrip('.root') for d in os.listdir(path_in) if not skipsample(d) ])
        logging.info('INFO: region=%s samples: %s', region, samples_region)

        files_in = {}  # mapping: sample <str> --> file <ROOT.TFile>

        # retrieve the full list of variables in this region. Open input files but don't close them
        variables_set = set()
        for sample in samples_region:
            files_in[sample] = ROOT.TFile(os.path.join(path_in, sample+'.root'))
            variables_set.update([ key.GetName() for key in files_in[sample].GetListOfKeys() if key.ReadObj().Class().InheritsFrom("TH1") and isVarSystematic(key.GetName())])
        variables_region = sorted(variables_set)

        # logging.debug('in {} 'the variables are: {}'.format(region, variables_region))
        # Add data
        if(args.unblind):
            try:
                data_obs = ROOT.TFile(os.path.join(path_in, 'data.root'))
            except OSError as e:
                logging.warning('While opening %s, caught %s', data_obs.GetName(), e)
            else:
                files_in['data_obs'] = data_obs

        # Write fake_photons
        if((not 'fake_photons' in files_in) or args.remake_fake_photons):
            logging.info('recreating fake_photons file: %s', os.path.join(path_in, 'fake_photons.root'))
            with TFileContext(os.path.join(path_in, 'fake_photons.root'), 'RECREATE') as fFakePh:
                fFakePh.cd()
                for variable in variables_region:
                    split = variable.split('_')
                    var_name = split[1]

                    if('loose' in var_name):
                        reweight_name = var_name.replace('loose', 'failReweight')
                        full_name = variable.replace(var_name, reweight_name)
                        h = files_in['data_obs'].Get(full_name)
                        if(h):
                            h.SetName( h.GetName().replace(reweight_name, var_name) )
                            h.Write()
        try:
            fFakePh = ROOT.TFile(os.path.join(path_in, 'fake_photons.root'))
        except OSError as e:
            logging.warning('While opening %s, caught %s', fFakePh.GetName(), e)
        else:
            files_in['fake_photons'] = fFakePh

        # Write to output
        with TFileContext(os.path.join(path_out, region+'.root'), "RECREATE") as fout:
            for variable in variables_region:
                split = variable.split('_')
                var_split = split[1].split('-')
                var_name = var_split[0]
                prompt = '-'+var_split[1] if len(var_split) > 1 else ''
                syst = split[2]
                if(var_name.endswith('failReweight')):
                    continue
                if(syst == 'central'):
                    skipIfData = False if len(var_split) == 1 else True
                    out_name = '{sample}{prompt}'.format(sample='%s', prompt=prompt)
                else:
                    skipIfData = True
                    direction = split[3]
                    out_name = '{sample}{prompt}_CMS_{syst}{direction}'.format(sample='%s', prompt=prompt, syst=syst, direction=direction)

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
                        ok_retrieved.append( {'file':files_in[sample].GetName(), 'variable':variable})
                    else:
                        not_retrieved.append({'file':files_in[sample].GetName(), 'variable':variable})

        del samples_region
        for _, handler in files_in.items():
            handler.Close()

    if(args.verbosity >= 1):
        files_prob = { e['file'] for e in not_retrieved }  # set()
        max_len = max([len(f) for f in files_prob])
        format_str = 'From file {:%d.%ds} could not retrieve {:d}/{:d} plots' % (max_len, max_len)
        for file_prob in sorted(files_prob):
            problems = [e['variable'] for e in not_retrieved if e['file'] == file_prob]
            good     = [e['variable'] for e in  ok_retrieved if e['file'] == file_prob]
            logging.warning(format_str.format(file_prob, len(problems), len(problems)+len(good)))

    if(args.verbosity >= 2):
        hists_prob         = { e['variable'] for e in not_retrieved }
        hists_prob_central = { e for e in hists_prob if e.endswith('central') }
        hists_prob_updn    = { e.rstrip('_Up').rstrip('_Down') for e in hists_prob if e.endswith(('Up', 'Down')) }
        for hist_prob in sorted(hists_prob_central):
            problems = [ e['file'] for e in not_retrieved if e['variable'] == hist_prob ]
            good     = [ e['file'] for e in  ok_retrieved if e['variable'] == hist_prob ]
            msg = 'Histogram {:40.40s} was missing {:2d}/{:2d} times'.format(hist_prob            , len(problems), len(good)+len(problems))
            if(len(problems) < 10):
                msg += ': '+' '.join([f.split('/')[-1] for f in problems])
            logging.warning(msg)

        for hist_prob in sorted(hists_prob_updn):
            problems = [ e['file'] for e in not_retrieved if e['variable'].startswith(hist_prob) ]
            good     = [ e['file'] for e in  ok_retrieved if e['variable'].startswith(hist_prob) ]
            msg = 'Histogram {:40.40s} was missing {:2d}/{:2d} times'.format(hist_prob+'(Up/Down)', len(problems), len(good)+len(problems))
            if(len(problems) < 10):
                msg += ': '+' '.join([f.split('/')[-1] for f in problems])
            logging.warning(msg)

    logging.info('Retrieved and wrote {:d} histograms. {:d} were missing. Total: {:d}'.format(len(ok_retrieved), len(not_retrieved), len(ok_retrieved)+len(not_retrieved)))

if __name__ == '__main__':
    main()
