#!/usr/bin/env python

#############################################################################################
# Get histograms from an analyzer's results (file: sample, inside has histograms: variable) #
# and write it in a form usable by Combine (file: variable, inside histograms: sample       #
# Wuthor A. Mecca                                                                           #
#############################################################################################

from __future__ import print_function
import os, sys
import ROOT
from plotUtils import TFileContext, makedirs_ok
from plotUtils23 import retrieve_bin_edges, InputDir
from subprocess import call
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
import logging
from utils23 import lumi_dict
import re


# Utility functions
def skipsample(filename):
    if '201' in filename: return True
    if filename == 'data.root': return True
    if filename in ('ggTo4l', 'ttXY', 'triboson'): return True
    if filename.split('.')[-1] != 'root':
        logging.error("strange sample:", filename)
        return True
    if re.search(r'part\d', filename):
        logging.info('Skipping %s', filename)
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

def write_fake_photons(fFakePh, data_obs, variables): # <TFile>, <TFile>, <iterable> of <str>
    logging.info('recreating fake_photons file: %s', fFakePh.GetName())
    fFakePh.cd()
    for variable in variables:
        split = variable.split('_')
        var_name = split[1]

        if('loose' in var_name):
            reweight_var_name = var_name.replace('loose', 'failReweight')
            reweight_split = [split[0]] + [reweight_var_name] + split[2:]
            reweight = '_'.join(reweight_split)
            logging.debug('variable=%40s  reweight=%40s', variable, reweight)
            h = data_obs.Get(reweight)
            if(h):
                h.SetName(variable)
                h.Write()
            else:
                logging.warning('No failReweight histogram in data for variable %s --> could not retrieve %s', var_name, reweight)
    logging.info('Wrote fake photons file to %s', fFakePh.GetName())

def get_TH1keys_from_file(tfhandle):
    return [ key.GetName()
             for key in tfhandle.GetListOfKeys()
             if ROOT.TClass(key.GetClassName()).InheritsFrom("TH1") and isVarSystematic(key.GetName()) ]

def main():
    # The configuration
    regions = ['SR4P', 'CR3P1F' , 'CR2P2F' , 'SR4P_1L', 'SR4P_1P', 'CR4P_1F', 'CR4L',
               'SR3P', 'CR110'  , 'CR101'  , 'CR011'  , 'CR100'  , 'CR001'  , 'CR010', 'CR000', 'SR3P_1L', 'SR3P_1P', 'CR3P_1F', 'CRLFR', 'CR3L',
               'SR2P', 'SR2P_1L', 'SR2P_1P', 'CR2P_1F'
               # 'SR_HZZ', 'CR2P2F_HZZ', 'CR3P1F_HZZ', 'CR_HZZ', 'MC_HZZ',
               # 'MC'
    ]

    parser = ArgumentParser(formatter_class=ArgumentDefaultsHelpFormatter)
    parser.add_argument('-y', '--year'     , default='2018', choices=lumi_dict.keys())
    parser.add_argument(      '--blind'    , action='store_true', help='Do not write data_obs in output files')
    parser.add_argument('-i', '--inputdir' , default='results', help='Top level directory where the results of analyzers are stored')
    parser.add_argument('-o', '--outputdir', default='histogramsForCombine', help='Output location')
    parser.add_argument('-A', '--analyzer' , default='VVGammaAnalyzer', help='Name of the analyzer, used to compose the path of the input files')
    parser.add_argument('-r', '--regions'  , default=['SR4P'], nargs='+', choices=regions, metavar='REGION', help=' ')
    parser.add_argument(      '--remake-fake-photons', action='store_true', help='Force to recreate the fake_photons.root file from data.root')
    parser.add_argument('-v', '--verbose'  , dest='verbosity', action='count', default=1, help='Increase the verbosity level')
    parser.add_argument('-q', '--quiet'    , dest='verbosity', action='store_const', const=0)
    parser.add_argument('--log', dest='loglevel', metavar='LEVEL', default='WARNING', help='Level for the python logging module. Can be either a mnemonic string like DEBUG, INFO or WARNING or an integer (lower means more verbose).')
    args = parser.parse_args()
    args.unblind = not args.blind

    loglevel = args.loglevel.upper() if not args.loglevel.isdigit() else int(args.loglevel)
    logging.basicConfig(format='%(levelname)s:%(module)s:%(funcName)s: %(message)s', level=loglevel)

    # Setup
    ok_retrieved  = []
    not_retrieved = []

    # Start
    path_out = os.path.join(args.outputdir, args.year)
    makedirs_ok(path_out)
    
    for region in args.regions:
        path_in  = InputDir(basedir=args.inputdir, year=args.year, region=region, analyzer=args.analyzer).path()
        if(not os.path.isdir(path_in)):
            logging.warning('Skipping non-existent dir: '+path_in)
            continue

        samples_region = set([ d.rstrip('.root') for d in os.listdir(path_in) if not skipsample(d) ])
        logging.info('INFO: region=%s samples: %s', region, samples_region)

        files_in = {}  # mapping: sample <str> --> file <ROOT.TFile>

        # Add data
        if(args.unblind):
            try:
                data_obs = ROOT.TFile(os.path.join(path_in, 'data.root'))
            except OSError as e:
                logging.error('While opening %s, caught %s\nIt will be skipped and data_obs will NOT appear in the output', data_obs.GetName(), e)
            else:
                files_in['data_obs'] = data_obs

        # Try to open fake_photons
        fake_photons_fname = os.path.join(path_in, 'fake_photons.root')
        try:
            fFakePh = ROOT.TFile(fake_photons_fname)
        except OSError as e:
            logging.info('While opening %s, caught %s', fake_photons_fname, e)
        else:
            files_in['fake_photons'] = fFakePh

        # Write fake_photons
        if((not 'fake_photons' in files_in) or args.remake_fake_photons):
            logging.info('Recreating fake_photons: %s', fake_photons_fname)
            with TFileContext(fake_photons_fname, 'RECREATE') as fFakePh:
                variables_data = get_TH1keys_from_file(files_in['data_obs'])
                write_fake_photons(fFakePh, data_obs=files_in['data_obs'], variables=variables_data)
            files_in['fake_photons'] = ROOT.TFile(fake_photons_fname)

        # If the program was run just to remake fake_photons, exit now
        if(args.remake_fake_photons):
            logging.info('Remade fake_photons, now exiting')
            for _, handle in files_in.items(): handle.Close()
            return 0

        # Open all the MC file handles and retrieve the full list of variables in this region
        variables_set = set()
        for sample in samples_region:
            files_in[sample] = ROOT.TFile(os.path.join(path_in, sample+'.root'))
            variables_set.update(get_TH1keys_from_file(files_in[sample]))
        variables_region = sorted(variables_set)

        logging.info('in %s there are %d variables', region, len(variables_region))
        logging.debug('in {} the variables are: {}'.format(region, variables_region))

        # Sometimes the yield in data.root may be 0. In this case we must insert an empty histogram with the appropriate xaxis
        xbins_dict = dict()

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

                ordered_files_in = [[k, v] for k,v in files_in.items() if k != 'data_obs']
                ordered_files_in.append(['data_obs', files_in['data_obs']])
                for sample, file_in in ordered_files_in:
                    if(sample == 'data_obs' and skipIfData):
                        continue
                    h = file_in.Get(variable)
                    if(h):
                        h.SetName(out_name %(sample))
                        h.Write()
                        ok_retrieved.append( {'file':file_in.GetName(), 'variable':variable})
                        if(syst == 'central'):  # Save bin edges in case we need it
                            xbins = retrieve_bin_edges(h.GetXaxis())
                            xbins_dict.setdefault(variable, xbins)
                        h.SetDirectory(0)
                        del h
                    else:
                        not_retrieved.append({'file':file_in.GetName(), 'variable':variable})
                        if(sample == 'data_obs'):  # data_obs must not be missing
                            logging.error('data_obs (%s) is missing "%s" - replacing with empty histogram', file_in.GetName(), variable)
                            xbins = xbins_dict[variable]
                            h = ROOT.TH1F(out_name %(sample), '', len(xbins) - 1, xbins)
                            h.Write()
                            h.SetDirectory(0)
                            del h

        logging.debug('Closing files')
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
            if(len(problems)/2 < 10):
                msg += ': '+' '.join({f.split('/')[-1] for f in problems})
            logging.warning(msg)

    logging.info('Retrieved and wrote {:d} histograms. {:d} were missing. Total: {:d}'.format(len(ok_retrieved), len(not_retrieved), len(ok_retrieved)+len(not_retrieved)))
    return 0

if __name__ == '__main__':
    exit(main())
