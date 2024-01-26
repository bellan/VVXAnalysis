#!/bin/sh

show_help(){ cat <<EOF
Usage: ${0##*/} [-u] CARD"
Convert CARD to workspace and assess the impact of the nuisance
parameters with fits

    -u      unblind
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
	    unblind=1
	    ;;
	r)
	    [ "$OPTARG" -eq "$OPTARG" ] || exit 2
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
fit_options="-m 125 --robustFit 1"
if [ $unblind -eq 0 ] ; then
    fit_options="$fit_options --expectSignal=$mu -t -1"
    outname="impacts_expected_$cardname"
else
    outname="impacts_observed_$cardname"
fi

mkdir -p $cardname && cd $cardname || { echo "Unable to make dir $cardname" 1>&2 ; exit 2 ; }

echo "### text2workspace ###"
text2workspace.py ${card} -o workspace.root || print_error "text2workspace"

echo "### Performing initial fit ###"
combineTool.py -M Impacts -d workspace.root $fit_options --doInitialFit || print_error "Initial fit"

echo "### Performing robust fit ###"
combineTool.py -M Impacts -d workspace.root $fit_options --doFits || print_error "Robust fit"

echo "### Extracting impacts ###"
combineTool.py -M Impacts -d workspace.root -m 125 -o impacts.json || print_error "Impacts extraction"

echo "### Plotting impacts ###"
plotImpacts.py -i impacts.json -o ${outname} || print_error "Plotting impacts"

echo "### Convert to png ###"
convert -density 300 ${outname}.pdf -trim ${outname}.png || print_error "Conversion to PNG"
