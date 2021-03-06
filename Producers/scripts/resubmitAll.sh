#!/bin/sh

###############################################################
#  Loops through subdirectories to find jobs-containing ones  #
#  (they contain a condor.sub file) and resubmits everything  #
#  using resubmit_Condor.csh                                  #
#                                                             #
#  Author: A. Mecca  (alberto.mecca@cern.ch)                  #
###############################################################

dirs=$(find . -maxdepth 4 -name condor.sub | grep -oP ".+(?=/condor.sub)" | grep -oP "(?<=\./).+")

for d in $dirs ; do
    (
	printf "%s \t--> " $d 
	cd $d && ls | grep -q Chunk &&  # Resubmit only if there are still chunks to do
	echo "DO" && resubmit_Condor.csh || echo "Done"
    )
done
