#!/bin/sh

set -e
set -u
set -o pipefail

show_help(){ cat <<EOF
Usage: ${0##*/} DIR"
    hadd parts of result files of an EventAnalyzer in DIR
EOF
}

haddOpt="-k -f"
dryrun=false
OPTIND=1
while getopts "hd" opt; do
    case $opt in
	h)
	    show_help
	    exit 0
	    ;;
	d)
	    dryrun=true
	    ;;
	\?)
	    show_help 1>&2
	    exit 1
	    ;;
    esac
done
shift "$((OPTIND-1))"

[ $# -eq 1 ] || { show_help 1>&2 ; exit 1 ; }
$dryrun && HADD="echo hadd" || HADD=hadd

# Target files WITHOUT the ".root" suffix
targets="$(find $1 -mindepth 3 -type f -name '*_part*.root' | sed -r 's:_part.+$::' | sort | uniq)"

for target in $targets ; do
    $HADD $haddOpt $target.root ${target}_part*.root \
	|| continue
	# && rm ${target}_part*.root \
done
