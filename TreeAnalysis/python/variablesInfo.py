#!/usr/bin/env python2

#####################################################################################################
# Single entry point to get a list of plots produced by the analysis, which can be used for data/MC #
# The idea is that each analysis mantains its own file, and this file remains relatively stable     #
#                                                                                                   #
# Author: A. Mecca (alberto.mecca@cern.ch)                                                          #
#####################################################################################################

from variablesInfo_VVGamma import getVarInfo_VVGamma
from variablesInfo_VVX import getVarInfo_VVX
from variablesInfo_ZZ import VarInfo_zz
from variablesInfo_VBS import VarInfo_vbs


def getVariablesInfo(analyzer, region):
    if   analyzer.startswith("ZZ"     ):
        VarInfo = VarInfo_zz
    elif analyzer.startswith("VVX"    ):
        VarInfo = getVarInfo_VVX(region)
    elif analyzer.startswith("VVGamma"):
        VarInfo = getVarInfo_VVGamma(region)
    else                               :
        VarInfo = getVarInfo_vvx(region)
    return VarInfo


if __name__ == '__main__':
    # test: print the varInfo
    from sys import argv
    from json import dumps
    
    if(len(argv) > 2):
        analyzer = argv[1]
        region   = argv[2]
    else:
        analyzer = 'VVGamma'
        region   = 'SR4P'
    
    print("TEST: analyzer =", analyzer, ", region =", region)
    VarInfo = getVariablesInfo(analyzer, region)
    print(dumps(VarInfo, indent=2))
