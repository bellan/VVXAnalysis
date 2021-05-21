#!/bin/sh

##############################################################
#  Creates the necessary symlinks in TreeAnalysis/python     #
#  and in Producers/test/pyFragments if they're not present  #
#                                                            #
#  Author: A. Mecca (amecca@cern.ch)                         #
##############################################################

guide="https://github.com/bellan/VVXAnalysis#packages-for-a-multi-boson-final-state-analysis"

VVXpython=TreeAnalysis/python
ZZpython=ZZAnalysis/AnalysisStep/python
samInfo=readSampleInfo.py

VVXpyF=Producers/test/pyFragments
ZZpyF=ZZAnalysis/AnalysisStep/test/prod/pyFragments
pyF_files="LHEProbabilities_GG_BKG_MCFM.py LHEProbabilities_VBF_BSI_0PM_H125_PhantomMCFM.py RecoProbabilities.py"


# Test if VVXAnalysis was correctly downloaded
[ -d $VVXpyF ] || { printf "%s folder in VVXanalysis not found.\nTry git pull and/or follow instructions at %s\n" $$VVXpyF $guide >&2 ; exit 1 ; }
[ -d $VVXpython ] || { printf "%s folder in VVXanalysis not found.\nTry git pull and/or follow instructions at %s\n" $VVXpython $guide >&2 ; exit 1 ; }

# Test if ZZAnalysys is where we expect it to be
[ -d ../ZZAnalysis ] || { printf "ZZAnalysis is needed. Follow instructions at %s\n" $guide >&2 ; exit 2 ; }


# Linking pyFragments
fileNotFound=0
for file in $pyF_files ; do
    [ -f ../$ZZpyF/$file ] || { printf "%s not found in ZZAnalysis\n" $file  >&2 ; fileNotFound=3 ; continue ; }
    if [ -e $VVXpyF/$file ] ; then  # There's already something: it should be a symlink to a file in ZZAnalysis
        [ ! -L $VVXpyF/$file ] && { printf "$VVXpyF/$file already exists but is not a symlink" >&2 ; }
    else
        echo "Linking $VVXpyF/$file --> ../../../../$ZZpyF/$file"
	ln -s ../../../../$ZZpyF/$file $VVXpyF/
    fi
    [ -L $VVXpyF/$file ] && [ -e $VVXpyF/$file ] || printf "Error: did not create a working symlink for $file\n" >&2
done
[ $fileNotFound -eq 0 ] || exit $fileNotFound

# Linking ReadSampleInfo.py
[ -f ../$ZZpython/$samInfo ] || { printf "%s not found in ZZAnalysis\n" $samInfo >&2 ; exit 3 ; }
if [ -e $VVXpython/$samInfo ] ; then  # There's already something: it should be a symlink to a file in ZZAnalysis
    [ ! -L $VVXpython/$samInfo ] && { printf "$VVXpython/$samInfo already exists but is not a symlink" >&2 ; }
else
    echo "Linking $VVXpython/$samInfo --> ../../../../$ZZpython/$samInfo"
    ln -s ../../../../$ZZpython/$samInfo $VVXpython/
fi
[ -L $VVXpython/$samInfo ] && [ -e $VVXpython/$samInfo ] || printf "Error: did not create a working symlink for $samInfo\n" >&2
