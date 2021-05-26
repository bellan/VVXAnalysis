#!/bin/bash

################################################################
#  Loops through subdirectories to find jobs-containing ones   #
#  (they contain a condor.sub file), then checks if there are  #
#  still unfinished jobs. Otherwise, it hadds all the chunks   #
#  in AAAOK using haddChunks.py                                #
#  It also synchronizes the hadd-ed file to a folder in eos    #
#                                                              #
#  Author: A. Mecca  (alberto.mecca@cern.ch)                   #
################################################################


eosdir="/eos/user/a/amecca/samples/May2021"
dirs=$(find . -maxdepth 4 -name condor.sub | grep -oP ".+(?=/condor.sub)" | grep -oP "(?<=\./).+" )

dirsExist () {
    for d in $@ ; do [ -d $d ] || return 1 ; done
}

for dir in $dirs ; do
    (
	cd $dir
	printf "%s \t--> " $dir
	if ls | grep -q Chunk ; then
	    echo "Still not done"
	else
	    samples=$(ls AAAOK/ | grep -oP ".+(?=_Chunk)" | uniq)
	    echo $samples | sed "s/\n/ /g" | xargs printf "%s "
	    dirsExist $(echo "$samples" | sed "s:^:AAAOK/:g") && echo "(Already hadd-ed)" || haddChunks.py AAAOK
	    for sampl in $samples ; do
		[ -e AAAOK/$sampl/ZZ4lAnalysis.root ] || exit 1
		stripDir=$(echo $dir | sed -r "s:/[^/]+$::g")
		destFile="$eosdir/$stripDir/$sampl.root"
		[ -e $destFile ] || rsync -au --progress AAAOK/$sampl/ZZ4lAnalysis.root $destFile
	    done
	fi
    )
done
