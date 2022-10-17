#!/bin/sh

###################################################################
#  Creates a readable summary of the output of checkProd.csh for  #
#  each of the subfolders in the folders passed as arguments      #
#                                                                 #
#  Author: A. Mecca (amecca@cern.ch)                              #
###################################################################

source _findJobDirs.sh ; jobdirs=$(findJobDirs $@)  # All subfolders of the arguments containing a file named condor.sub

#[ $# -lt 1 ] && echo "Usage: ${0##*/} FOLDER[S]" && exit 1
tempfile=$(mktemp)
tempfail=$(mktemp)

get_cause() {
    case $1 in
	-1)  echo "Error return without specification" ;;
	1)   echo "Hangup (POSIX)" ;;
	2)   echo "Interrupt (ANSI)" ;;
	4)   echo "Illegal instruction (ANSI)" ;;
	5)   echo "Trace trap (POSIX)" ;;
	6)   echo "Abort (ANSI) or IOT trap (4.2BSD)" ;;
	7)   echo "BUS error (4.2BSD)" ;;
	8)   echo "Floating point exception (ANSI)" ;;
	9)   echo "killed, unblockable (POSIX) kill -9" ;;
	11)  echo "segmentation violation (ANSI) (most likely user application crashed)" ;;
	15)  echo "Termination (ANSI)" ;;
	24)  echo "Soft CPU limit exceeded (4.2 BSD)" ;;
	25)  echo "File size limit exceeded (4.2 BSD)" ;;
	30)  echo "Power failure restart (System V.)" ;;
	50)  echo "Required application version not found at the site." ;;
	64)  echo "I/O error: cannot open data file (SEAL)" ;;
	65)  echo "End of job from user application (CMSSW)" ;;
	66)  echo "Application exception" ;;
	73)  echo "Failed writing to read-only file system" ;;
	81)  echo "Job did not find functioning CMSSW on worker node." ;;
	84)  echo "Some required file not found; check logs for name of missing file." ;;
	85)  echo "Error reading file" ;;
	86)  echo "Fatal ROOT error (probably while reading file)" ;;
	90)  echo "Application exception" ;;
	92)  echo "Failed to open file" ;;
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
	152|-152) echo "CPU limit exceeded (4.2 BSD)" ;;
	153|-153) echo "File size limit exceeded (4.2 BSD)" ;;
	154|-154) echo "You cancelled the job" ;;
	155) echo "Profiling alarm clock (4.2 BSD)" ;;
	0)   echo "Still running (or unknown failure)" ;;
	*)   echo "unknown" ;;
    esac
}

totToDo=0
for jobdir in $jobdirs ; do
    echo "----- $jobdir -----"
    prodSummary.sh "$jobdir" | tee $tempfile | sed "s/^/\t/g"
    grep -q TODO $tempfile && totToDo=$(( $totToDo + $(grep -oP "(?<=TODO: )\d+" $tempfile) ))
    lines=$(grep failed $tempfile)
    [ -n "$lines" ] && echo "$lines" | sed -e "s/--> failed, exit status = //g" -e "s/^\t//g" -e "s/^ \+//g" >> $tempfail
    running=$(grep "still running" $tempfile | grep -oP "\d+(?= -->)")
    [ -n "$running" ] && echo "$running 0" >> $tempfail
done

printf "\nTotal jobs to do: %d\n" $totToDo  # $(cut -d " " -f 1 $tempfail | paste -s -d+ | bc)

todoKnown=0
keys=$(cut -d " " -f 2 $tempfail | sort | uniq)

for key in $keys ; do
    jobs_k=$(grep -oP "\d+(?= $key)" $tempfail | paste -s -d+ | bc)
    todoKnown=$(( $todoKnown + $jobs_k ))
    printf "\t%d:\t%d \t%s\n" $key $jobs_k "$(get_cause $key)"
done
todoUnknown=$(( $totToDo - $todoKnown )) 
[ $todoUnknown -gt 0 ] && printf "\t%s:\t%d \t%s\n" "-" $todoUnknown "No info"

rm -f $tempfile $tempfail
