#!/bin/sh

[ $# -eq 1 ] || { echo "Usage: $0 DATACARD" >&2 ; exit 2 ; }
[ -f $1 ] || { echo "$1 does not exist or is not a file" >&2 ; exit 2 ; }
card=$(realpath $1)
echo $card
cardname=${1##*/}
cardname=${cardname%.txt}
original_dir=$(pwd -P)

mkdir -p $cardname && cd $cardname || exit 1

echo "### text2workspace ###"
text2workspace.py ${card} -o workspace.root || exit 1

echo "### Performing initial fit ###"
combineTool.py -M Impacts -d workspace.root -m 125 --doInitialFit --robustFit 1 -t -1 --expectSignal=1 || exit 1

echo "### Performing robust fit ###"
combineTool.py -M Impacts -d workspace.root -m 125 --robustFit 1 --doFits -t -1 --expectSignal=1 || exit 1

echo "### Extracting impacts ###"
combineTool.py -M Impacts -d workspace.root -m 125 -o impacts.json || exit 1

echo "### Plotting impacts ###"
cd ${original_dir}
plotImpacts.py -i ${cardname}/impacts.json -o ${cardname}_impacts || exit 1

echo "### Convert to png ###"
convert -density 300 ${cardname}_impacts.pdf -trim ${cardname}_impacts.png || exit 1
