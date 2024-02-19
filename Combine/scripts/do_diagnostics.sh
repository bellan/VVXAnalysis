#!/bin/sh

show_help(){ cat <<EOF 
Usage: ${0##*/} [-u] CARD"
Run FitDiagnostics on CARD and extract some results. Currently diffNuisances.py
and mlfitNormsToText.py
EOF
}

print_error(){ printf "%s failed (%d).\n" "$1" $? ; exit 2 ; }

OPTIND=1
while getopts "hn:" opt; do
    case $opt in
	h)
	    show_help
	    exit 0
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
combine_testdir=${CMSSW_BASE?}/src/HiggsAnalysis/CombinedLimit/test

mkdir -p $cardname && cd $cardname || exit 1

combine -M FitDiagnostics --robustFit=1 --saveNormalizations --saveWithUncertainties "$card" || print_error "FitDiagnostics"

python ${combine_testdir}/diffNuisances.py --all fitDiagnosticsTest.root -f latex > diffNuisances_$cardname.tex || print_error "diffNuisances"
python ${combine_testdir}/mlfitNormsToText.py -u fitDiagnosticsTest.root > fitNorms_$cardname.txt || print_error "mlfitNormsToText"
cut -c 42- fitNorms_$cardname.txt \
    | sed -r -e 's/-{3,}/\\midrule/' \
    -e "s/(   {2,})/\1\& /g" \
    -e 's/ +& *$/ \\\\/g' \
    -e 's:\+/-:\\pm:g' \
> fitNorms_$cardname.tex || print_error "fitNorms post-processing"
