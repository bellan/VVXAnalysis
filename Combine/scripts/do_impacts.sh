#!/bin/sh

set -e
set -u

show_help(){ cat <<EOF
Usage: ${0##*/} [-u] CARD"
Convert CARD to workspace and assess the impact of the nuisance
parameters with fits

    -u      unblind (once for pulls, twice for pulls and signal strength)
    -r VAL  set the expected limit to VAL; ignored if -u is set
    -n VAL  set the number of points in the scan (default: 50)
EOF
}

print_error(){ printf "%s failed (%d).\n" "$1" $? ; exit 2 ; }

unblind=0
mu=1
# npoints=50
OPTIND=1
while getopts "hur:" opt; do
    case $opt in
	h)
	    show_help
	    exit 0
	    ;;
	u)
	    unblind=$(($unblind+1))
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

card=$(realpath $1)
echo $card
cardname=${1##*/}
cardname=${cardname%.txt}
original_dir=$(pwd -P)

T2Woptions="-m 125 --for-fits --no-wrappers --optimize-simpdf-constraints=cms --X-pack-asympows --use-histsum -v 0 -o workspace.root"
fit_options="-M Impacts -d workspace.root -m 125 -v 0 --rMin -1 --rMax 3 --robustFit 1"
plot_options=""
exclude="--exclude rgx{^prop_bin}"

if [ $unblind -eq 0 ] ; then
    fit_options="$fit_options --expectSignal=$mu -t -1"
    plot_options="$plot_options --blind"
    outname="impacts_expected_$cardname"
elif [ $unblind -eq 1 ] ; then
    plot_options="$plot_options --blind"
    outname="impacts_unblind_$cardname"
else
    outname="impacts_observed_$cardname"
fi

initial_options="$fit_options"
robust_options="$fit_options $exclude"
extract_options="$fit_options $exclude"


mkdir -p $cardname && cd $cardname || { echo "Unable to make dir $cardname" 1>&2 ; exit 2 ; }

echo "### text2workspace ###"
text2workspace.py $T2Woptions ${card} || print_error "text2workspace"

echo "### Performing initial fit ###"
combineTool.py $initial_options --doInitialFit || print_error "Initial fit"

echo "### Performing robust fit ###"
combineTool.py $robust_options --doFits || print_error "Robust fit"

echo "### Extracting impacts ###"
combineTool.py $extract_options -o $outname.json || print_error "Impacts extraction"

echo "### Plotting impacts ###"
fix_postfit_pull.py --log info $outname.json
plotImpacts.py -i $outname.json -o ${outname} $plot_options || print_error "Plotting impacts"

echo "### Convert to png ###"
convert -density 300 ${outname}.pdf -trim ${outname}.png || print_error "Conversion to PNG"
