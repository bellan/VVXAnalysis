#!/usr/bin/env python2

from __future__ import print_function
from os import path
from copy import deepcopy
# from json import dump
import ROOT
from plotUtils import TFileContext, getSamplesByRegion, getPlot


print('REGION  |  ZW   |  l1   |  l2   |  l3   | contrib')
for region in ['CR000', 'CR001', 'CR010', 'CR100', 'CR011', 'CR101', 'CR110']:
    h_ZW = getPlot('debug3L_ZW_FRSF', 'data', region)
    h_l1 = getPlot('debug3L_l1_FRSF', 'data', region)
    h_l2 = getPlot('debug3L_l2_FRSF', 'data', region)
    h_l3 = getPlot('debug3L_l3_FRSF', 'data', region)

    f_ZW = h_ZW.GetMean()
    f_l1 = h_l1.GetMean()
    f_l2 = h_l2.GetMean()
    f_l3 = h_l3.GetMean()
    n_ZW = f_ZW * h_ZW.GetEntries()

    # print('{:7s} |{:+.3f}\t|{:+.3f}\t|{:+.3f}\t|{:+.3f}\t|{:+6.0f}'.format(region, f_ZW, f_l1, f_l2, f_l3, n_ZW))
    def getSign(n):
        if  (n>0): return '+'
        elif(n<0): return '-'
        else     : return '0'
    print('{:7s} |   {}   |   {}   |   {}   |   {}   |{:+6.0f}'.format(region, getSign(f_ZW), getSign(f_l1), getSign(f_l2), getSign(f_l3), n_ZW))


# h_001    = getPlot('debug3L_ZW_mass'   , 'data', 'CR001')
# h_010    = getPlot('debug3L_ZW_mass'   , 'data', 'CR010')
# h_100    = getPlot('debug3L_ZW_mass'   , 'data', 'CR100')

# h_000_w1 = getPlot('debug3L_ZW_mass_w1', 'data', 'CR000')
# h_000_w2 = getPlot('debug3L_ZW_mass_w2', 'data', 'CR000')
# h_000_w3 = getPlot('debug3L_ZW_mass_w3', 'data', 'CR000')

# h_001.Add(h_000_w1, -1)
# h_010.Add(h_000_w2, -1)
# h_100.Add(h_000_w3, -1)


# h_000 = getPlot('debug3L_ZW_mass_w1', 'data', 'CR000')
# h_001 = getPlot('debug3L_ZW_mass_w1', 'data', 'CR001')
# h_010 = getPlot('debug3L_ZW_mass_w1', 'data', 'CR010')
# h_100 = getPlot('debug3L_ZW_mass_w1', 'data', 'CR100')
# h_011 = getPlot('debug3L_ZW_mass_w1', 'data', 'CR011')
# h_101 = getPlot('debug3L_ZW_mass_w1', 'data', 'CR101')
# h_110 = getPlot('debug3L_ZW_mass_w1', 'data', 'CR110')
