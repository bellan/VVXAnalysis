#!/bin/sh

###############################################################
#  Loops through subdirectories to find jobs-containing ones  #
#  (they contain a condor.sub file) and resubmits everything  #
#  using resubmit_Condor.csh                                  #
#                                                             #
#  Author: A. Mecca  (alberto.mecca@cern.ch)                  #
###############################################################

source _findJobDirs.sh ; jobdirs=$(findJobDirs $@)  # All subfolders of the arguments containing a file named condor.sub

for d in $dirs ; do
    (
	printf "%s \t--> " $d 
	cd $d && ls | grep -q Chunk &&  # Resubmit only if there are still chunks to do
	echo "DO" && cleanup.csh && { resubmit_Condor.csh || echo "Done" ; }
    )
done
