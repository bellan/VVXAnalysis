#!/bin/bash

set -e
set -u
set -o pipefail

show_help(){ cat <<EOF
Usage: ${0##*/} DIR"
    hadd data result files of an EventAnalyzer in DIR to produce a data.root
    for each region (in each year), and fake_leptons.root in SR4P and SR3P.
EOF
}

haddOpt="-k -f"
OPTIND=1
while getopts "h" opt; do
    case $opt in
	h)
	    show_help
	    exit 0
	    ;;
	\?)
	    show_help 1>&2
	    exit 1
	    ;;
    esac
done
shift "$((OPTIND-1))"

[ $# -eq 1 ] || { show_help 1>&2 ; exit 1 ; }

basenames(){
    # read stdin if present, arguments ($@) otherwise
    [ -p /dev/stdin ] && input="$(cat -)" || input=$@
    for arg in $input; do
	basename "$arg"
    done
}

# Data
for yeardir in $(find "$1" -maxdepth 1 -mindepth 1 -type d -name "20*" | sed s:^./::g); do
    echo "INFO: year = " $(basename "$yeardir")
    dirs=$(find $yeardir -maxdepth 1 -mindepth 1 -type d)
    analyzer=$(basenames $dirs | cut -d _ -f 1 | sort | head -n 1) # TODO loop on multiple analyzers
    regions=$(basenames $dirs | cut -d _ -f 2)

    # Hadd data streams to create data.root in each region of each year
    for dir in $dirs; do
    	data_samples=$(ls $dir/*.root | grep -P '((Single|Double)(Ele|Mu|EG)|EGamma|Mu(on)?EG)' | sort)
        echo hadd $haddOpt $dir/data.root $(printf "$dir/%s " ${data_samples})  # printf reuses the format string to consume all of its arguments
    done
    echo

    # fake_leptons in SR4P
    if echo "$regions" | grep -q -e CR3P1F -e CR2P2F ; then
        echo hadd $haddOpt $yeardir/${analyzer}_SR4P/fake_leptons.root \
    	    $yeardir/${analyzer}_CR3P1F/data.root \
    	    $yeardir/${analyzer}_CR2P2F/data.root
    fi
    echo

    # fake_leptons in SR3P
    if echo "$regions" | grep -q -P 'CR[10][10][10]' ; then
        echo hadd $haddOpt $yeardir/${analyzer}_SR3P/fake_leptons.root \
    	    $yeardir/${analyzer}_CR110/data.root \
    	    $yeardir/${analyzer}_CR101/data.root \
    	    $yeardir/${analyzer}_CR011/data.root \
    	    $yeardir/${analyzer}_CR100/data.root \
    	    $yeardir/${analyzer}_CR010/data.root \
    	    $yeardir/${analyzer}_CR001/data.root \
    	    $yeardir/${analyzer}_CR000/data.root
    fi
    echo
done
