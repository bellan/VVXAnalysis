#!/bin/sh

set -e
set -u

show_help(){ cat <<EOF
Usage: ${0##*/} [-h] [-d] RESULTS_DIR"
    Produce datacards and histograms for Combine from the results of an EventAnalyzer
    This script must be called from VVXAnalysis/TreeAnalysis/

    -h   Show help and exit
    -d   Dry-run: print the commands that would be executed
EOF
}

[ $(pwd -P | rev | cut -d / -f 1-2 | rev) = VVXAnalysis/TreeAnalysis ] || { show_help >&2 ; exit 1 ; }

dryrun=false
OPTIND=1
while getopts "hd" opt; do
    case $opt in
        h)
            show_help
            exit 0
            ;;
        d)
            dryrun=true
            ;;
        *)
            echo "Unknown option \"$opt\"" >&2
            show_help >&2
            exit 1
            ;;
    esac
done
shift "$((OPTIND-1))"

[ $# -eq 1 ] || { show_help >&2 ; exit 1 ; }
res_dir_full="${1%/}"
$dryrun && EXEC=echo || EXEC=""

loglevel=error
Run2years="2016preVFP 2016postVFP 2017 2018"

printf "### Hadd result data ###\n"
$EXEC scripts/hadd_result_data.sh "$res_dir_full"

printf "### Hadd result parts ###\n"
$EXEC scripts/hadd_result_parts.sh "$res_dir_full"

printf "### Split triboson ###\n"
$EXEC python/split_triboson.py --log info -i "$res_dir_full"

for res_dir in "$res_dir_full" "${res_dir_full}"_triboson ; do
    res_name="${res_dir#results_}"
    for year in $Run2years ; do
	printf "\n### fake_photons.root %s year=%s ###\n" "$res_name" $year
	$EXEC python/prepareHistoCombine.py --remake-fake-photons --log $loglevel -r SR4P -y $year -i "${res_dir}"

	printf "\n### Systematics JSON %s year=%s ###\n" "$res_name" $year
        $EXEC python/systematics.py -A VVGammaAnalyzer --log $loglevel -r SR4P -y $year -i "${res_dir}"

        printf "\n### Prepare Histo Combine %s year=%s ###\n" "$res_name" $year
	$EXEC python/prepareHistoCombine.py --log $loglevel -r SR4P -y $year -i "${res_dir}" -o histogramsForCombine_"${res_name}"

	printf "\n### Make datacards %s year=%s ###\n" "$res_name" $year
	# for strategy in $(cat combine/strategiesSR4P) ; do
	for strategy in combine/*.json ; do
	    $EXEC python/produceDataCard_VVGamma.py -r SR4P --log $loglevel -y $year -i histogramsForCombine_"${res_name}" $strategy
	done
    done
    carddir=combine/cards_"${res_name}"
    [ -e "$carddir" ] && {
	echo "INFO: moving existing $carddir"
	$EXEC mv -v --no-clobber "$carddir" "${carddir}-old"
    }
    $EXEC mv --no-clobber combine/cards "$carddir"

    # Keep copies of systematics JSONs
    json_dir="data/systematics/${res_name}"
    $EXEC mkdir -p "$json_dir"
    $EXEC cp --no-clobber data/systematics_*.json "$json_dir/"
    printf "\n### Done %s ###\n" "$res_name"
done
