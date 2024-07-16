#!/bin/sh

set -e
set -u
set -o pipefail

show_help(){ cat <<EOF
Usage: ${0##*/} SOURCE [DEST]
    Gather the results and logs of Condor analysis jobs
    and copies them to DEST (default: .)
EOF
}

OPTIND=1
while getopts "h" opt; do
    case $opt in
	h)
	    show_help
	    exit 0
	    ;;
	\?)
	    show_help 1>&2
	    exit 1
	    ;;
    esac
done
shift "$((OPTIND-1))"

[ $# -ge 1 ] || { show_help 1>&2 ; exit 1 ; }
source="$1"
dest="${2:-.}"

echo "INFO: copying result rootfiles..."
find "$source" -mindepth 4 -maxdepth 5 -type d -name results | sed s:$:/: | xargs -I {} rsync -a --info=progress2 {} "$dest"

echo "INFO: copying logs..."
destlogdir="$dest"/logdir
mkdir -p "$destlogdir"
logs="$(find $source -mindepth 4 -maxdepth 5 -type f -name run.log)"

for log in $logs ; do
    localname=$(echo $log | sed -r "s:^$source/?::")
    year="$(echo $localname | cut -d / -f 1)"
    sample="$(echo $localname | cut -d / -f 3)"
    destpath="$destlogdir"/"${sample}"_"${year}".log
    cp --no-clobber "$log" "$destpath"
done

echo "INFO: done"

# Compress log dir
tar -c -z --remove-files -f "$destlogdir".tar.gz "$destlogdir"
