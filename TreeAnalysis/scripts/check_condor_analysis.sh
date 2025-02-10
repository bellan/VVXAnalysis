#!/usr/bin/sh

set -u
set -o pipefail

_DONE=0
_SUBM=1
_WARN=2
_FAIL=3
_UNKN=4
_ERROR=5

show_help(){ cat <<EOF 
Usage: ${0##*/} DIR"
    Check the production status of HTCondor jobs for EventAnalyzer in DIR
      0: DONE
      1: SUBMITTED
      2: WARNING
      3: FAILED
      4: UNKNOWN
EOF
}

# analyzer=VVGammaAnalyzer
OPTIND=1
while getopts "hA:" opt; do
    case $opt in
	h)
	    show_help
	    exit 0
	    ;;
	# A)
	#     analyzer=$OPTARG
	#     ;;
	*)
	    show_help >&2
	    exit 1
	    ;;
    esac
done
shift "$((OPTIND-1))"

[ $# -ge 1 ] && proddir="$1" || proddir=.  #{ show_help >&2 ; exit 1 ; }
[ -d $(realpath "$proddir") ] || { echo "$proddir does not exist or is not a directory" >&2 ; exit 2 ; }

cd "$proddir"

test_log(){
    local chunk="$1"
    local logfile="$(find $chunk/log -name '*.log')"
    [ -e "$logfile" ] || return $_UNKN
    [ $(wc -l $logfile | cut -d " " -f 1) -lt 3 ] && { echo "INFO: submitted" ; return $_SUBM ; }
    grep -q SYSTEM_PERIODIC_REMOVE "$logfile" && { echo "ERROR: $chunk SYSTEM_PERIODIC_REMOVE" && return $_FAIL ; }
    return $_UNKN
}

test_chunk(){
    # Return 0 if chunk is ok
    local chunk="$1"
    local logstatus=$_UNKN
    [ -e "$chunk" ] || { echo "ERROR: $chunk does not exist" ; return $_ERROR ; }
    if ! [ -e "$chunk"/exitStatus.txt ] ; then
	test_log "$chunk" ; logstatus=$?
	[ $logstatus -eq $_UNKN ] && \
	    echo "ERROR: $chunk missing exitStatus.txt"
	return $logstatus
    fi
    local status=$(cat "$chunk"/exitStatus.txt)
    [ $status -eq 0 ] || { echo "ERROR: $chunk exit status = $status" ; return $_FAIL ; }
    [ -e "$chunk"/results ] || { echo "ERROR: $chunk missing results" ; return $_FAIL ; }

    test_results "$chunk"/results || return $?
}

test_results(){
    local resultsdir="$1"
    local size=0
    for regiondir in "$resultsdir"/*/* ; do
	rootfile=$(find "$regiondir" -type f -name "*.root")
	[ -n "$rootfile" ] || { echo "WARNING: $regiondir no rootfile (no events of this sample in this region?)" ; return $_WARN ; }
	# rootfile="$(echo $rootfile | head -n1)"  # in the remote hypothesis that there is more than one
	size=$(stat -c "%s" "$rootfile")
	[ $size -gt 10000 ] || { echo "ERROR: $rootfile is too small ($size)" ; return $_FAIL ; }
    done
    return $_DONE
}

outDB=status.csv
# Backup old DB
[ -e $outDB ] && mv $outDB $outDB.bak
# Open DB on fd 3
exec 3>$outDB

for sampledir in $(find . -maxdepth 2 -mindepth 2 -type d) ; do
    sample="$(basename $sampledir)"
    # printf "%s (%s) " "$sampledir" "$sample"
    singlechunk=$sampledir/${sample}
    if [ -e $singlechunk ] ; then
	test_chunk $singlechunk
	printf "%s,%d\n" "$singlechunk" $? >&3
    else
	for chunk in $sampledir/${sample}_Chunk* ; do
	    test_chunk $chunk
	    printf "%s,%d\n" "$chunk" $? >&3
	done
    fi
done
