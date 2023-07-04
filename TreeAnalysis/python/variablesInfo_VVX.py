#!/usr/bin/env python

#####################################################################################################
# The basic dictionary of plots which most analyzer produce.                                        #
# This is used directly in variablesInfo, but can be also imported in each variablesInfo_[ANALYZER] #
#                                                                                                   #
# Author: A. Mecca (alberto.mecca@cern.ch)                                                          #
#####################################################################################################

def getVarInfo_VVX(region):
    VarInfo_VVX = {
        "AAA_cuts"  : {'title':'Cuts', 'text':True, 'unblind':True, 'logy':True, 'ymin':1},
        'channel_lep':{'title':'lepton flavour', 'text':True, 'unblind':True}
    }
    return VarInfo_VVX
