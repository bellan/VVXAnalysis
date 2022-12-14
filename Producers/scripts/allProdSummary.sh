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

trap "rm -- $tempfile $tempfail" EXIT

exitcode_descr() {
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
    grep ":" $tempfile | tr -d "\t" | sed "s/: /:/g" >> $tempfail
done

echo

done_n=$(grep -oP "\d+(?=:DONE)" $tempfail | paste -s -d+ | bc)
todo_n=$(grep -oP "\d+(?=:TODO)" $tempfail | paste -s -d+ | bc)
printf "Total jobs DONE: %d\n" $done_n
printf "Total jobs TODO: %d\n" $todo_n

todoKnown=0
# The "key" is the full description given by checkProd.csh, e.g.:
# "failed, exit status = 92 (failed to open file)"
# "still running (or unknown failure)"
while read -r key ; do
    lines_in_tempfail="$(grep "$key" $tempfail)"
    jobs_k=$(echo "$lines_in_tempfail" | grep -oP "\d+(?=:)" | paste -s -d+ | bc)
    todoKnown=$(( $todoKnown + $jobs_k ))
    # If checkProd.csh tried to guess the meaning of the exit status, correct it with the ones we know from https://twiki.cern.ch/twiki/bin/view/CMSPublic/StandardExitCodes
    status_n=$(echo "$lines_in_tempfail" | head -n 1 | grep -oP "(?<=failed, exit status = )\d+")
    if [ -n "$status_n" ] ; then
	description=$(exitcode_descr $status_n)
    else
	description="$key"
	status_n=0
    fi
    printf "\t%d:\t%d \t%s\n" $jobs_k $status_n "$description"
done <<EOF
$(cut -d ":" -f 2- $tempfail | sort | uniq | grep -v "TODO\|DONE")
EOF

todoUnknown=$(( $totToDo - $todoKnown )) 
[ $todoUnknown -gt 0 ] && printf "\t%s:\t%d \t%s\n" $todoUnknown "-" "No info"
