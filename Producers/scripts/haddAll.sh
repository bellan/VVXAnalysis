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

campaign=Dec2022
initial=$(printf "%.1s" $USER)
eosdir="/eos/user/$initial/$USER/samples/$campaign"

source _findJobDirs.sh ; jobdirs=$(findJobDirs $@)  # All subfolders of the arguments containing a file named condor.sub

dirsExist () {
    for d in $@ ; do [ -d $d ] || return 1 ; done
}

for dir in $jobdirs ; do
    (
	cd $dir
	printf "%s \t--> " $dir
	if ls | grep -q Chunk ; then
	    j_done=$( [ -d AAAOK ] && ls AAAOK/ | grep Chunk | wc -l || echo "0" )
	    j_todo=$(ls        | grep Chunk | wc -l)
	    printf "Still not done (%d/%d)\n" $j_done $(($j_done + $j_todo))
	else
	    samples=$(ls AAAOK/ | grep -oP ".+(?=_Chunk)" | uniq)
	    echo $samples | sed "s/\n/ /g" | xargs printf "%s "
	    dirsExist $(echo "$samples" | sed "s:^:AAAOK/:g") && echo "(Hadd-ed)" || haddChunks.py AAAOK
	    for sampl in $samples ; do
		[ -e AAAOK/$sampl/ZZ4lAnalysis.root ] || { echo "Tried to hadd but no rootfile has been created" 1>&2 ; exit 1 ; }
		stripDir=$(echo $dir | sed -r "s:/[^/]+$::g")
		mkdir -p "$eosdir/$stripDir"
		destFile="$eosdir/$stripDir/$sampl.root"
		[ -e $destFile ] || { rsync -au --progress AAAOK/$sampl/ZZ4lAnalysis.root $destFile && echo "Copied to $destFile" ; }
	    done
	fi
    )
done
