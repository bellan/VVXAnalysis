#!/bin/sh

###################################################################
#  Creates a readable summary of the output of checkProd.csh for  #
#  each of the subfolders in the folders passed as arguments      #
#                                                                 #
#  Author: A. Mecca (amecca@cern.ch)                              #
###################################################################

basedir=$(pwd -P)

if [ $# -lt 1 ] ; then
    jobdirs=$(find . -maxdepth 4 -name "*Chunk*" -prune -name condor.sub | grep -oP ".+(?=/condor.sub$)" | sed "s|^\./||g" )
else
    jobdirs=$(find $@ -maxdepth 4 -name "*Chunk*" -prune -name condor.sub | grep -oP ".+(?=/condor.sub$)" | sed "s|^\./||g" )
fi

#[ $# -lt 1 ] && echo "Usage: ${0##*/} FOLDER[S]" && exit 1
tempfile=$(mktemp)
tempfail=$(mktemp)

get_cause() {
    case $1 in
	-1) echo "Error return without specification" ;;
	1) echo "Hangup (POSIX)" ;;
	2) echo "Interrupt (ANSI)" ;;
	4) echo "Illegal instruction (ANSI)" ;;
	5) echo "Trace trap (POSIX)" ;;
	6) echo "Abort (ANSI) or IOT trap (4.2BSD)" ;;
	7) echo "BUS error (4.2BSD)" ;;
	8) echo "Floating point exception (ANSI)" ;;
	9) echo "killed, unblockable (POSIX) kill -9" ;;
	11) echo "segmentation violation (ANSI) (most likely user application crashed)" ;;
	15) echo "Termination (ANSI)" ;;
	24) echo "Soft CPU limit exceeded (4.2 BSD)" ;;
	25) echo "File size limit exceeded (4.2 BSD)" ;;
	30) echo "Power failure restart (System V.)" ;;
	50) echo "Required application version not found at the site." ;;
	64) echo "I/O error: cannot open data file (SEAL)" ;;
	65) echo "End of job from user application (CMSSW)" ;;
	66) echo "Application exception" ;;
	73) echo "Failed writing to read-only file system" ;;
	81) echo "Job did not find functioning CMSSW on worker node." ;;
	84) echo "Some required file not found; check logs for name of missing file." ;;
	85) echo "Job failed to open local and fallback files." ;;
	90) echo "Application exception" ;;
	92) echo "Job failed to open local and fallback files." ;;
	126) echo "Permission problem or command is not an executable" ;;
	127) echo "Command not found" ;;
	129) echo "Hangup (POSIX)" ;;
	132) echo "Illegal instruction (ANSI)" ;;
	133) echo "Trace trap (POSIX)" ;;
	134) echo "Abort (ANSI) or IOT trap (4.2 BSD) (most likely user application crashed)" ;;
	135) echo "Bus error (4.2 BSD)" ;;
	136) echo "Floating point exception (e.g. divide by zero)" ;;
	137) echo "SIGKILL; likely an unrelated batch system kill." ;;
	139) echo "Segmentation violation (ANSI)" ;;
	143) echo "Termination (ANSI)(or incorporate in the msg text)" ;;
	147) echo "Error during attempted file stageout." ;;
	151) echo "Error during attempted file stageout." ;;
	152) echo "CPU limit exceeded (4.2 BSD)" ;;
	153) echo "File size limit exceeded (4.2 BSD)" ;;
	155) echo "Profiling alarm clock (4.2 BSD)" ;;
	*)   echo "unknown" ;;
	    esac
}

for jobdir in $jobdirs ; do
    echo "----- $jobdir -----"
    prodSummary.sh "$jobdir" | tee $tempfile | sed "s/^/\t/g"
    #grep -q TODO $tempfile && totToDo=$(( $totToDo + $(grep -oP "(?<=TODO: )\d+" $tempfile) ))
    lines=$(grep failed $tempfile)
    [ -n "$lines" ] && echo "$lines" | sed "s/--> \tfailed, exit status = //g" >> $tempfail
    #printf "%d\t%d\n" $(echo "$lines" | grep -oP "\d+$") $(echo "$lines" | grep -oP "\d+(?= -->)") >> $tempfail
done

echo
echo "Total jobs to do:" $(cut -d " " -f 1 $tempfail | paste -s -d+ | bc)

keys=$(cut -d " " -f 2 $tempfail | sort | uniq)
for key in $keys ; do
    jobs_k=$(grep -oP "\d+(?= $key)" $tempfail | paste -s -d+ | bc)
    #printf "%d --> %d\n" $key "$jobs_k"
    printf "\t%d:\t%d \t%s\n" $key $jobs_k "$(get_cause $key)"
done

rm -f $tempfile $tempfail
