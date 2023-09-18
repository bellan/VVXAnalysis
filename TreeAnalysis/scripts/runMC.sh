#!/bin/sh -xu

haddOpt="-k -f"  # "-j $(grep processor /proc/cpuinfo | wc -l)"
rmOpt="-r -f"
analyzer=VVGammaAnalyzer # VVXAnalyzer
year=2018
nevents=-1

regions="SR4P;CR3P1F;CR2P2F;SR3P;CR110;CR101;CR011;CR100;CR010;CR001;CR000;SR2P;CRLFR"
options="--external-cross-section"  # --internal-cross-section --nofr


if [ ! -d logdir ] ; then
    # echo "logdir created"
    mkdir logdir
fi

while getopts "y:r:" option ; do
    case "$option" in
	y) year=$OPTARG ;;
	r) regions=$OPTARG ;;
    esac
done

echo -e "--- Analyses on MC ---\n"
git log -1 --oneline HEAD

make || exit

mc_samples=$(ls samples/MC/$year | grep -v $year | grep -oP ".+(?=\.root)" | grep -v "WLLGTo2L2j_5f_LO\|ZZTo4l_M1ToInf" | sort)

for sample in $mc_samples ; do
    ./python/run.py $analyzer $sample -r $regions -y $year -n $nevents -d samples/MC $options >logdir/${sample}_${year}.log 2>&1 &
done

wait

for region in $(echo $regions | tr ';' ' ') ; do
    (
	cd results/$year/${analyzer}_${region}
	hadd $haddOpt ggTo4l.root ggTo4e_Contin_MCFM701.root ggTo2e2mu_Contin_MCFM701.root ggTo4mu_Contin_MCFM701.root
	hadd $haddOpt triboson.root ZZZ.root WZZ.root WWZ.root
	hadd $haddOpt ttXY.root TTWW.root TTZZ.root
    )
done
