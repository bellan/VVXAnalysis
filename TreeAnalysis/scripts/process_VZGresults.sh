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
if $dryrun ; then
    prepareHistoCombine="echo python/prepareHistoCombine.py"
    VZGsystematics="echo python/VZGsystematics.py"
    produceDataCard_VZG="echo python/produceDataCard_VZG.py"
    MV="echo mv"
    CP="echo cp"
    MKDIR="echo mkdir"
else
    prepareHistoCombine=python/prepareVZGHistoCombine.py
    VZGsystematics=python/VZGsystematics.py
    produceDataCard_VZG=python/produceDataCard_VZG.py
    MV=mv
    CP=cp
    MKDIR=mkdir
fi
loglevel=error

Run2years="2016preVFP","2016postVFP","2017","2018"

for res_dir in "$res_dir_full" ; do
    res_name="${res_dir#results_}"
    for year in $Run2years ; do
	printf "\n### Systematics JSON %s year=%s ###\n" "$res_name" $year
        $VZGsystematics -A VZGAnalyzer --log $loglevel -r SR2P -y $year -i "${res_dir}"

        printf "\n### Prepare Histo Combine %s year=%s ###\n" "$res_name" $year
	$prepareHistoCombine --log $loglevel -r SR2P -y $year -i "${res_dir}" -o histogramsForCombine_"${res_name}"

	printf "\n### Make datacards %s year=%s ###\n" "$res_name" $year

	#for strategy in $(cat combine/strategiesSR2P) ; do
	for strategy in combine/*.json ; do
	    $produceDataCard_VZG -r SR2P --log $loglevel -y $year -i histogramsForCombine_"${res_name}" $strategy # combine/SR2P_$strategy.json
	done
    done
    carddir=combine/cards_"${res_name}"
    [ -e "$carddir" ] && { echo "INFO: moving existing $carddir" ; $MV -v --no-clobber "$carddir" "${carddir}-old" ; }
    $MV --no-clobber combine/cards "$carddir"

    # Keep copies of systematics JSONs
    json_dir="data/VZGsystematics/${res_name}"
    $MKDIR -p "$json_dir"
    $CP --no-clobber data/VZGsystematics_*.json "$json_dir/"
    printf "\n### Done %s ###\n" "$res_name"
done
