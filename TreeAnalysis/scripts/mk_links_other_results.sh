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

[ $# -eq 2 ] || { show_help >&2 ; exit 1 ; }
dsttop="$1"
srctop="$2"
$dryrun && EXEC="echo" || EXEC=""

ndone=0
nskip=0
for f in $(find "$srctop" -mindepth 3 -maxdepth 3 -type f -name "*.root") ; do
    relname=${f#*/}
    reldir=$(dirname $relname)
    dstdir=$dsttop/$reldir
    dstpath=$dsttop/$relname

    [ -e $dstpath ] && {
	nskip=$(($nskip+1))
	continue
    }

    $EXEC mkdir -p -v $dstdir || break
    $EXEC ln -s ../../../$f $dstdir || break
    ndone=$(($ndone+1))
done
