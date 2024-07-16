#!/usr/bin/sh

set -o pipefail

show_help(){ cat <<EOF 
Usage: ${0##*/} DIR"
    Check the production status of HTCondor jobs for EventAnalyzer in DIR
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
    # Return 0 if the cause of failure was identified
    local chunk="$1"
    local logfile="$(find $chunk/log -name '*.log')"
    [ -e "$logfile" ] || return 1
    grep -q SYSTEM_PERIODIC_REMOVE "$logfile" && { echo "ERROR: $chunk SYSTEM_PERIODIC_REMOVE" && return 0 ; }
    return 1
}

test_chunk(){
    # Return 0 if chunk is ok
    local chunk="$1"
    [ -e "$chunk" ] || { echo "ERROR: $chunk does not exist" ; return 1 ; }
    if ! [ -e "$chunk"/exitStatus.txt ] ; then
	test_log "$chunk" && return 1
	echo "ERROR: $chunk missing exitStatus.txt"
	return 1
    fi
    local status=$(cat "$chunk"/exitStatus.txt)
    [ $status -eq 0 ] || { echo "ERROR: $chunk exit status = $status" ; return 1 ; }
    [ -e "$chunk"/results ] || { echo "ERROR: $chunk missing results" ; return 1 ; }

    test_results "$chunk"/results || return
}

test_results(){
    local resultsdir="$1"
    local size=0
    for regiondir in "$resultsdir"/*/* ; do
	rootfile=$(find "$regiondir" -type f -name "*.root")
	[ -n "$rootfile" ] || { echo "WARNING: $regiondir no rootfile (no events of this sample in this region?)" ; return 0 ; }
	# rootfile="$(echo $rootfile | head -n1)"  # in the remote hypothesis that there is more than one
	size=$(stat -c "%s" "$rootfile")
	[ $size -gt 10000 ] || { echo "ERROR: $rootfile is too small ($size)" ; return 1 ; }
    done
    return 0
}

for sampledir in $(find . -maxdepth 2 -mindepth 2 -type d) ; do
    sample="$(basename $sampledir)"
    # printf "%s (%s) " "$sampledir" "$sample"
    singlechunk=$sampledir/${sample}
    if [ -e $singlechunk ] ; then
	test_chunk $singlechunk
    else
	for chunk in $sampledir/${sample}_Chunk* ; do
	    test_chunk $chunk
	done
    fi
done
