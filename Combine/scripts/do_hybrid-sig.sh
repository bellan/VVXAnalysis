#!/bin/sh

set -e
set -u

show_help(){ cat <<EOF
Usage: ${0##*/} [-u] CARD"
Convert CARD to workspace compute the significance as recommended in
https://cms-analysis.github.io/HiggsAnalysis-CombinedLimit/latest/part3/commonstatsmethods/#computing-significances-with-toys

    -u      unblind (once for pulls, twice for pulls and signal strength)
    -r VAL  set the expected signal strength to VAL; ignored if -u is set
EOF
}

print_error(){ printf "%s failed (%d).\n" "$1" $? ; exit 2 ; }

dryrun=false
unblind=false
mu=1
outdir=.
OPTIND=1
while getopts "hdur:o:" opt; do
    case $opt in
	h)
	    show_help
	    exit 0
	    ;;
	d)
	    dryrun=true
	    echo "WARN: some commands redirect their output to files"
	    ;;
	o)
	    outdir=$OPTARG
	    ;;
	u)
	    unblind=true
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
$dryrun && EXEC=echo || EXEC=""

card=$(realpath $1)
echo $card
cardname=${1##*/}
cardname=${cardname%.txt}
original_dir=$(pwd -P)

common_options="-M HybridNew --LHCmode LHC-significance"
gentoys_options="$common_options --saveToys --fullBToys --saveHybridResult" # -T toys -i iterations
signif_options="$common_options --readHybridResult"
if $unblind ; then
    prefix=observed
else
    prefix=expected
    gentoys_options="$gentoys_options --expectSignal $mu"
fi

workdir=${outdir}/${prefix}/${cardname}
echo "INFO: working in $workdir"
$EXEC mkdir -p $workdir && $EXEC cd $workdir || { echo "Unable to make dir $workdir" 1>&2 ; exit 2 ; }

echo "### Generating toys ###"
echo "DEBUG: flags = $gentoys_options -s \$seed"
for seed in $(seq 1 16); do
    [ -e higgsCombineTest.HybridNew.mH120.$seed.root ] && continue
    echo "# seed $seed #"
    $EXEC combine $gentoys_options -s $seed $card &> toys$seed.log &
done
wait || print_error "Generating toys"

echo "### Hadd-ing toys ###"
[ -f merged_toys.root ] || \
    $EXEC hadd -f merged_toys.root higgsCombineTest.HybridNew.mH120.*.root || print_error "Hadd-ing toys" 

echo "### Significance (unblind: $unblind) ###"
echo "DEBUG: flags = $signif_options --toysFile=merged_toys.root"
if $unblind ; then
    $EXEC combine $signif_options --toysFile=merged_toys.root $card > signif.log # || print_error "Observed signif"
    $EXEC sed -n "/Significance/,//p" signif.log > significance.txt
else
    for expect in 0.16 0.5 0.84 ; do
	{
	    echo "DEBUG:         --expectedFromGrid=$expect"
	    $EXEC combine $signif_options --expectedFromGrid=$expect --toysFile=merged_toys.root $card > signif_$expect.log
	    $EXEC sed -n "/Significance/,//p" signif_$expect.log > significance_$expect.txt
	} &
    done
fi
wait
