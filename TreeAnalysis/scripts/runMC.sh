#!/bin/sh -x

prefixcsvfile="../Producers/python/samples_"
suffixcsvfile="_MC.csv"
haddOpt="-k -f"  # "-j $(grep processor /proc/cpuinfo | wc -l)"
rmOpt="-r -f"
analyzer=VVGammaAnalyzer # VVXAnalyzer
years="2016"
nevents=-1

regions="SR4P;CR2P2F;CR3P1F;SR3P;CR000;CR001;CR010;CR011;CR100;CR101;CR110;SR2P"

options="--external-cross-section" #--internal-cross-section --fpw  --nofr


if [ ! -d logdir ] ; then
    # echo "logdir created"
    mkdir logdir
fi

echo -e "--- Analyses on MC ---\n"

make || exit

for year in $years ; do
    csvfile=$prefixcsvfile${year}UL$suffixcsvfile
    mc_samples=$(ls samples/MC/$year | grep -v $year | grep -oP ".+(?=\.root)" | grep -v "WLLGTo2L2j_5f_LO\|ZZTo4l_M1ToInf" | sort)

    for sample in $mc_samples ; do
 	# echo "Running $analyzer on $sample in region $regions with $csvfile $nevents $options"
	./python/run.py $analyzer $sample -r $regions -c $csvfile -n $nevents -d samples/MC $options &> logdir/${sample}_${year}.log &
    done
    
    wait
    
    for region in $(echo $regions | tr ';' ' ') ; do
	(
	    cd results/$year/${analyzer}_${region}
    	    hadd $haddOpt ggTo4l.root ggTo4e_Contin_MCFM701.root ggTo2e2mu_Contin_MCFM701.root ggTo4mu_Contin_MCFM701.root
    	    [ -L ggZZ.root ] || ln -s ggTo4l.root ggZZ.root
	    hadd $haddOpt triboson.root ZZZ.root WZZ.root WWZ.root
	    hadd $haddOpt ttXY.root TTWW.root TTZZ.root
	)
    done
    
done

wait