#!/usr/bin/env python2

######################################################################################################################################################
# Produce tables with info about the number of events in histograms                                                                                  #
# For simplicity, one must explicitly edit the __main__ function to select which tables to print                                                     #
#                                                                                                                                                    #
# Author: A. Mecca (alberto.mecca@cern.ch)                                                                                                           #
######################################################################################################################################################

from __future__ import print_function
from os import path
from copy import deepcopy
from math import log10, floor
import pandas as pd
import ROOT
from plotUtils import TFileContext, getSamplesByRegion


def _getPlot(filename, var):
    try:
        with TFileContext(filename) as tf:
            retrieved = deepcopy(tf.Get(var))
    except (IOError, OSError) as e:
        # ROOT prints an error message on its own
        retrieved = None
    else:
        if(retrieved is None):
            print('WARN: plot "{}" not present in file "{}"'.format(var, filename))
    return retrieved


def getPlot(var, sample, region, inputdir='results', year='2016', analyzer='VVGammaAnalyzer'):
    filename = path.join(inputdir, year, '{:s}_{:s}'.format(analyzer, region), '{:s}.root'.format(sample))  # '{inputdir}/{year}/{analyzer}_{region}/{sample}.root'
    return _getPlot(filename, var)


def getPlots_added(var, samples, region, inputdir='results', year='2016', analyzer='VVGammaAnalyzer'):
    # Do getPlot() on multiple files and Add() the results
    result = None
    for sample in samples:
        h = getPlot(var, sample, region, inputdir, year, analyzer)
        if(h is not None):
            if(result is None):
                result = h
            else:
                result.Add(h)
    return result


def _my_formatter(df, column):
    the_max = df[column].max()
    decimals = max(0, 3 - floor( log10(the_max) ))
    string = '{:.%df}'%(decimals)
    # print('>>>', the_max, '-', decimals, '-', string)
    return string.format


def table1Plot(var, region, efficiencyType='cutflow'):
    samplesFullMC = getSamplesByRegion(region, 'pow', 'fullMC')
    samplesFromCR = getSamplesByRegion(region, 'pow', 'fromCR')
    backgrsFullMC = [ b for b in samplesFullMC if set(b['files']).isdisjoint(set(['WZGTo3LNuG', 'ZZGTo4LG', 'WWW'])) ]
    backgrsFromCR = [ b for b in samplesFromCR if set(b['files']).isdisjoint(set(['WZGTo3LNuG', 'ZZGTo4LG', 'WWW'])) ]
    
    hData   = getPlot(var, 'data', region)
    if region in ['SR4P', 'CR3P1F' , 'CR2P2F' , 'SR4P_1L', 'SR4P_1P', 'CR4P_1F', 'CR4L']:
        signal_sample = 'ZZGTo4LG'
    elif region in ['SR3P', 'CR110'  , 'CR101'  , 'CR011'  , 'CR100'  , 'CR001'  , 'CR010', 'CR000', 'SR3P_1L', 'SR3P_1P', 'CR3P_1F', 'CRLFR', 'CR3L']:
        signal_sample = 'WZGTo3LNuG'
    elif region in ['SR2P', 'SR2P_1L', 'SR2P_1P', 'CR2P_1F']:
        print('ERROR: sample for 2L not available yet')
        return
    else:
        print('ERROR: don\'t know what is the signal for region', region)
        return
    
    hSignal = getPlot(var, signal_sample, region)
    if(hSignal is None):
        print('ERROR: no signal plot: var={}, sample={}, region={}'.format(var, signal_sample, region))
        return
    
    hBackMC = None
    for bkg in backgrsFullMC:
        h = getPlots_added(var, bkg['files'], region)
        if(h):
            if(hBackMC is None):
                hBackMC = h
            else:
                hBackMC.Add(h)
            # print('INFO: Added MC:', bkg['files'])

    hBackCR = None
    for bkg in backgrsFromCR:
        h = getPlots_added(var, bkg['files'], region)
        if(h):
            if(hBackCR is None):
                hBackCR = h
            else:
                # print('hBackCR =', hBackCR, '  h =', h)
                hBackCR.Add(h)
            # print('INFO: Added from CR:', bkg['files'])
    
    infoTable = {}
    infoTable.setdefault(var, [])
    for b in range(1, hSignal.GetNbinsX()+1):
        labelSignal = hSignal.GetXaxis().GetBinLabel(b)
        labelBackMC = hBackMC.GetXaxis().GetBinLabel(b)
        labelBackCR = hBackCR.GetXaxis().GetBinLabel(b)
        if(labelBackMC != labelSignal):
            print('WARN: label fo bin {:d} differs from signal and fullMC background'.format(b))
            continue
        if(labelBackCR != labelSignal):
            print('WARN: label fo bin {:d} differs from signal and MC+CR background'.format(b))
            continue
    
        if efficiencyType == 'cutflow':
            if(len(infoTable[var]) == 0):
                eSignal = 1.
            else:
                eSignal = nSignal/infoTable[var][0]['signal']
            nSignal = hSignal.GetBinContent(b)
            nBackMC = hBackMC.GetBinContent(b)
            nBackCR = hBackCR.GetBinContent(b)
            nData   = hData  .GetBinContent(b)
            s_fullMC = nSignal/nBackMC
            s_backCR = nSignal/nBackCR
            label = labelSignal
        elif efficiencyType == 'cut_low':  # cut events lower than this bin
            nSignal = hSignal.Integral(b, -1)
            nBackMC = hBackMC.Integral(b, -1)
            nBackCR = hBackCR.Integral(b, -1)
            nData   = hData  .Integral(b, -1)
            eSignal = nSignal/hSignal.Integral(0, -1)
            s_fullMC = nSignal/nBackMC
            s_backCR = nSignal/nBackCR
            label = '{:g}-{:g}'.format(hSignal.GetBinLowEdge(b), hSignal.GetBinLowEdge(b+1))
            
        infoTable[var].append({ 'label':label, 'signal':nSignal, 'back fullMC':nBackMC, 'back CR+MC':nBackCR,
                                'S/fullMC': s_fullMC, 'S/(MC+CR)':s_backCR, 'eSignal': eSignal, 'data':nData})
    
    
    for var, table in infoTable.items():
        title = '{} - {}'.format(var, region)
        nPadding = (84 - len(title) - 2)/2
        print('-'*nPadding, title, '-'*nPadding)
        dictionary = { r['label']: {k:v for k,v in r.items() if k != 'label'} for r in table }
        df = pd.DataFrame.from_dict(dictionary, orient='index')
        df = df[['signal', 'eSignal', 'data', 'back fullMC', 'back CR+MC', 'S/fullMC', 'S/(MC+CR)']]
        df = df.reindex([r['label'] for r in table])
        
        formatters={c:_my_formatter(df,c) for c in df.columns}
        formatters.update({'eSignal': lambda x : '{:.1f}%'.format(100*x)})
        print(df.to_string(formatters=formatters))
        
        # print('         LABEL          |   Sig    | effS [%] | Bkg(MC)  |  MC+CR   |   S/MC   | S/MC+CR  |')
        # for row in table:
        #     numbers = [ '{:3.3g}'.format(f) for f in
        #                 [row['signal'], 100*row['eSignal'], row['backMC'], row['backCR'], row['S/MC'], row['S/MC+CR']]
        #               ]
        #     print( '{:23.23s}\t| {:8.8s} | {:>8.8s} | {:>8.8s} | {:>8.8s} | {:>8.8s} | {:>8.8s} |'.format(*[row['label']]+numbers) )
        
        print('-'*84)


def tableRegions(var, integration_range=[1,-1]):
    infoTable = {}
    for region in ['SR3P', 'CR110', 'CR101', 'CR011', 'CR100', 'CR010', 'CR001', 'CR000']:
        samplesFullMC = getSamplesByRegion(region, 'pow', 'fullMC')
        samplesFromCR = getSamplesByRegion(region, 'pow', 'fromCR')
        backgrsFullMC = [ b for b in samplesFullMC if set(b['files']).isdisjoint(set(['WZGTo3LNuG', 'ZZGTo4LG', 'WWW'])) ]
        backgrsFromCR = [ b for b in samplesFromCR if set(b['files']).isdisjoint(set(['WZGTo3LNuG', 'ZZGTo4LG', 'WWW'])) ]
        
        nSignal = getPlot(var, 'WZGTo3LNuG', region).Integral(*integration_range)
        nData   = getPlot(var, 'data'      , region).Integral(*integration_range)
        histsMC = [getPlots_added(var, b['files'] , region) for b in backgrsFullMC]
        histsCR = [getPlots_added(var, b['files'] , region) for b in backgrsFromCR]
        nBackMC = sum(h.Integral(*integration_range) for h in histsMC if h is not None)
        nBackCR = sum(h.Integral(*integration_range) for h in histsCR if h is not None)
        
        infoTable[region] = {'signal': nSignal, 'data':nData, 'back fullMC': nBackMC, 'back CR+MC': nBackCR}

    title = '{} (range: {})'.format(var, integration_range.__repr__())
    nPadding = (45 - len(title) - 2)/2
    print('-'*nPadding, title, '-'*nPadding)
    df = pd.DataFrame.from_dict(infoTable)
    with pd.option_context('display.float_format', '{:+3.1e}'.format):
        print(df.transpose())
    # print(' | '.join(['{:9.9s}']*4).format('REGION', ' signal', ' bkg MC', '  MC+CR'))
    # for region, row in infoTable.items():
    #     print(' | '.join( ['{:9.9s}'] + ['{: >+9.3g}']*3 ).format(region, row['signal'], row['backMC'], row['backCR']))
    print('-'*45)

if __name__ == '__main__':
    tableRegions('AAA_cuts', integration_range=[1,1])
    table1Plot('AAA_cuts'  , 'SR3P')
    table1Plot('lll_mass'  , 'SR3P', efficiencyType='cut_low')
    # table1Plot('AAA_cuts'  , 'SR3P')
    # table1Plot('AAA_cuts'  , 'CR110')
    # table1Plot('WZ_cutflow', 'CR110')
