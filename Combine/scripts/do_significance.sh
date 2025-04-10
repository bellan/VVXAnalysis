#!/bin/sh

show_help(){ cat <<EOF
Usage: ${0##*/} [-h] [-q] [-u] CARD"
Convert CARD to workspace and assess the impact of the nuisance
parameters with fits

    -h      show help and exit
    -q      filter the output of Combine and print only the result
    -u      unblind
    -r VAL  set the expected limit to VAL; ignored if -u is set
EOF
}

quiet=false
unblind=0
mu=1
OPTIND=1
while getopts "hqur:" opt; do
    case $opt in
	h)
	    show_help
	    exit 0
	    ;;
	q)
	    quiet=true
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

card="$1"
fit_options=
if [ $unblind -eq 0 ] ; then
    fit_options="$fit_options --expectSignal=$mu -t -1"
fi

# Use a global variable to hold the status of Combine, to circumvent
# the absence of $PIPESTATUS in POSIX sh
combinestatus=1 # default to fail; set it to 0 only if it completes succesfully
run(){
    combine -M Significance $fit_options "$card"
    combinestatus=$?
    return $combinestatus
}

$quiet && { run | grep "^Significance" ; exit $combinestatus ; } || run
