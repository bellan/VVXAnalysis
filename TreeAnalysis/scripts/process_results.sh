#!/bin/sh

set -e
set -u

show_help(){ cat <<EOF 
Usage: ${0##*/} RESULTS_DIR"
    Produce datacards and histograms for Combine from the results of an EventAnalyzer
    This script must be called from VVXAnalysis/TreeAnalysis/
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
if $dryrun ; then
    split_triboson="echo python/split_triboson.py"
    prepareHistoCombine="echo python/prepareHistoCombine.py"
    systematics="echo python/systematics.py"
    produceDataCard_VVGamma="echo python/produceDataCard_VVGamma.py"
    MV="echo mv"
    CP="echo cp"
    MKDIR="echo mkdir"
else
    split_triboson=python/split_triboson.py
    prepareHistoCombine=python/prepareHistoCombine.py
    systematics=python/systematics.py
    produceDataCard_VVGamma=python/produceDataCard_VVGamma.py
    MV=mv
    CP=cp
    MKDIR=mkdir
fi
loglevel=error

Run2years="2016preVFP 2016postVFP 2017 2018"

printf "### Split triboson ###\n"
python/split_triboson.py --log info -i res_dir_aeta2p4-leptonplots

for res_dir in "$res_dir_full" ; do #"${res_dir_full}"_triboson ; do
    res_name="${res_dir#results_}"
    for year in $Run2years ; do
	printf "\n### fake_photons.root %s ###\n" "$res_name"
	$prepareHistoCombine --remake-fake-photons --log $loglevel -r SR4P -y $year -i "${res_dir}"

	printf "\n### Systematics JSON %s ###\n" "$res_name"
        $systematics -A VVGammaAnalyzer --log $loglevel -r SR4P -y $year -i "${res_dir}"

        printf "\n### Prepare Histo Combine %s ###\n" "$res_name"
	$prepareHistoCombine --log $loglevel -r SR4P -y $year -i "${res_dir}" -o histogramsForCombine_"${res_name}"

	printf "\n### Make datacards %s ###\n" "$res_name"
	for strategy in $(cat combine/strategiesSR4P) ; do
	    $produceDataCard_VVGamma -r SR4P --log $loglevel -y $year -i histogramsForCombine_"${res_name}" combine/SR4P_$strategy.json
	done
    done
    $MV --no-clobber combine/cards combine/cards_"${res_name}"

    # Keep copies of systematics JSONs
    json_dir="data/systematics/${res_name}"
    $MKDIR -p "$json_dir"
    $CP --no-clobber data/systematics_*.json "$json_dir/"
    printf "\n### Done %s ###\n" "$res_name"
done
