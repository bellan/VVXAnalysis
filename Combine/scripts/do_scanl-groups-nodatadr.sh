#!/bin/sh

show_help(){ cat <<EOF 
Usage: ${0##*/} [-u] CARD"
Convert CARD to workspace and fit the signal strength while 
sequentially freezing the parameter groups.

    -u      unblind
    -r VAL  set the expected limit to VAL; ignored if -u is set
    -n VAL  set the number of points in the scan (default: 50)
EOF
}

print_error(){ printf "%s failed (%d).\n" "$1" $? ; exit 2 ; }

unblind=0
mu=1
minr=0
maxr=2
npoints=50
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
	n)
	    [ "$OPTARG" -eq "$OPTARG" ] || exit 2  # Ensure that OPTARG is a number
	    npoints=$OPTARG
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

mkdir -p $cardname && cd $cardname || exit 1

fit_options="--saveWorkspace --algo grid --setParameterRanges r=$minr,$maxr --robustFit=1 --points=$npoints"
if [ $unblind -eq 0 ] ; then
    fit_options="$fit_options --expectSignal=$mu -t -1"
    outname="scan_expected_$cardname"
else
    outname="scan_observed_$cardname"
fi

echo "### text2workspace ###"
text2workspace.py "$card" -o workspace.root || print_error "text2workspace"

echo "### Performing initial fit ###"
combine workspace.root -M MultiDimFit $fit_options -n .postfit || print_error "Initial fit"
rootpattern=higgsCombine.%s.MultiDimFit.mH120.root
rootpostfit=$(printf $rootpattern postfit)

echo "### Fitting total ###"
combine $rootpostfit -M MultiDimFit --snapshotName MultiDimFit $fit_options -n .total || print_error "Fit total"

echo "### Fitting freeze_lumi ###"
combine $rootpostfit -M MultiDimFit --snapshotName MultiDimFit $fit_options -n .freeze_lumi --freezeNuisanceGroups lumi || print_error "Fit freeze_lumi"

echo "### Fitting freeze_lumi_theory ###"
combine $rootpostfit -M MultiDimFit --snapshotName MultiDimFit $fit_options -n .freeze_lumi_theory --freezeNuisanceGroups lumi,theory || print_error "Fit freeze_lumi_theory"

echo "### Fitting freeze_all ###"
combine $rootpostfit -M MultiDimFit --snapshotName MultiDimFit $fit_options -n .freeze_all --freezeParameters allConstrainedNuisances || print_error "Fit freeze_all"

echo "### Plotting ###"
plot1DScan.py $(printf $rootpattern total) --main-label "Total Uncert." --others \
    $(printf $rootpattern freeze_lumi):"freeze lumi":3 \
    $(printf $rootpattern freeze_lumi_theory):"freeze lumi+theory":7 \
    $(printf $rootpattern freeze_all):"stat only":2 \
    --output $outname --y-max 10 --y-cut 40 --breakdown "lumi,theory,rest,stat" \
    || print_error "Plotting"

# mv breakdown.png "$outname.png"
# mv breakdown.pdf "$outname.pdf"
