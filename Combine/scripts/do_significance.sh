#!/bin/sh

show_help(){ cat <<EOF
Usage: ${0##*/} [-u] CARD"
Convert CARD to workspace and assess the impact of the nuisance
parameters with fits

    -u      unblind
    -r VAL  set the expected limit to VAL; ignored if -u is set
EOF
}

unblind=0
mu=1
OPTIND=1
while getopts "hur:" opt; do
    case $opt in
	h)
	    show_help
	    exit 0
	    ;;
	u)
	    unblind=1
	    ;;
	r)
	    [ "$OPTARG" -eq "$OPTARG" ] || exit 2
	    mu=$OPTARG
	    ;;
	*)
	    show_help >&2
	    exit 1
	    ;;
    esac
done
shift "$((OPTIND-1))"

[ $# -eq 1 ] || { show_help >&2 ; exit 1 ; }
[ -f $1 ] || { echo "$1 does not exist or is not a file" >&2 ; exit 2 ; }

fit_options=
if [ $unblind -eq 0 ] ; then
    fit_options="$fit_options --expectSignal=$mu -t -1"
fi

combine -M Significance $fit_options $1
