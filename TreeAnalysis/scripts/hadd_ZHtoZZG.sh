#!/bin/sh

set -e
set -u

show_help(){
cat <<EOF
Usage: ${0##*/} [-d] DEST_DIR SOURCE_DIR
    Create symlinks in DEST_DIR to the results in SRC_DIR

    -d      Dry-run: print commands that would be executed
EOF
}

dryrun=false
OPTIND=1
while getopts "dh" opt; do
    case $opt in
	h)  show_help
	    exit 0 ;;
	d)  dryrun=true ;;
	*)  show_help >&2
	    exit 1
	    ;;
    esac
done
shift "$((OPTIND-1))"

[ $# -eq 1 ] || { show_help >&2 ; exit 1 ; }
top="$1"
$dryrun && EXEC="echo" || EXEC=""

for d in $(find "$top" -mindepth 2 -maxdepth 2 -type d -name "VVGammaAnalyzer_*") ; do
    [ -e $d/signal.root ] && continue

    $EXEC hadd $d/signal.root $d/ZZGTo4LG.root $d/ZHtoZZG.root || break
done
