#!/bin/sh

show_help(){ cat <<EOF
Usage: ${0##*/} [-u] CARD"
Convert CARD to workspace and assess the Goodness of Fit

    -r VAL  set the expected limit to VAL; ignored if -u is set
    -n VAL  set the number of points in the scan (default: 50)
EOF
}

print_error(){ printf "%s failed (%d).\n" "$1" $? ; exit 2 ; }

mu=1
ntoys=1000
OPTIND=1
while getopts "hr:n:" opt; do
    case $opt in
	h)
	    show_help
	    exit 0
	    ;;
	r)
	    mu=$OPTARG
	    ;;
	n)
	    [ "$OPTARG" -eq "$OPTARG" ] || exit 2  # Ensure that OPTARG is a number
	    ntoys=$OPTARG
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

card=$(realpath $1)
echo $card
cardname=${1##*/}
cardname=${cardname%.txt}
original_dir=$(pwd -P)
fit_options="-m 125 --robustFit 1 --cminApproxPreFitTolerance 0.01"
outname="gof_unblind_$cardname"

mkdir -p $cardname && cd $cardname || { echo "Unable to make dir $cardname" 1>&2 ; exit 2 ; }

echo "### text2workspace ###"
text2workspace.py -v 0 -o workspace.root ${card} || print_error "text2workspace"

echo "### Fit the data ###"
combine -M GoodnessOfFit -d workspace.root --algo saturated -n .goodnessOfFit_data || print_error "Fit data"

echo "### Fit the toys ###"
combine -M GoodnessOfFit -d workspace.root --algo saturated -t $ntoys -s -1 --toysFrequentist -n .goodnessOfFit_toys || print_error "Fit toys"

echo "### Collect the toys ###"
combineTool.py -M CollectGoodnessOfFit --input higgsCombine.goodnessOfFit_data*.root higgsCombine.goodnessOfFit_toys*.root -o gof.json || "Collect toys"

fitted_mass=$(python - <<EOF
import json
with open('gof.json') as f: d = json.load(f)
print(d.keys()[0])
EOF
)
echo "fitted_mass:" $fitted_mass

echo "### Plotting GoF ###"
plotGof.py gof.json --statistic saturated --mass $fitted_mass -o gof_$outname
