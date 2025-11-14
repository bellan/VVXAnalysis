#!/bin/sh

show_help(){ cat <<EOF 
Usage: ${0##*/} [-u] CARD"
Run FitDiagnostics on CARD and plot correlation matrix
EOF
}

print_error(){ printf "%s failed (%d).\n" "$1" $? ; exit 2 ; }

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

card="$(realpath $1)"
cardname="${1##*/}"
cardname="${cardname%.txt}"
nameNoAutoMC="${cardname}_noAutoMC"

fit_options="--robustFit=1 --saveNormalizations --saveShapes --saveWithUncertainties --plots"
if [ $unblind -eq 0 ] ; then
    fit_options="$fit_options --expectSignal=$mu -t -1"
fi

mkdir -p $cardname && cd $cardname || exit 1

# remove autoMC stats so they do not make the plot unreadable
sed -r '/^[^ ]+ autoMCStats( [0-9]+)+$/d' $card > "${nameNoAutoMC}.txt"

combine -M FitDiagnostics ${fit_options} "${nameNoAutoMC}.txt" || print_error "FitDiagnostics"

# Make a nice pdf with the correlation matrix
plotCorrMatrix.py fitDiagnosticsTest.root          && { mv covariance_fit_s.png covariance_fit_s_$cardname.png; mv covariance_fit_s.pdf covariance_fit_s_$cardname.pdf; }
plotCorrMatrix.py --b-only fitDiagnosticsTest.root && { mv covariance_fit_b.png covariance_fit_b_$cardname.png; mv covariance_fit_b.pdf covariance_fit_b_$cardname.pdf; }

mv -v fitDiagnosticsTest.root fitDiagnostics_"$nameNoAutoMC".root
