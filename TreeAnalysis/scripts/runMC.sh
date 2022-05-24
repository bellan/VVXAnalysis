#!/bin/sh

prefixcsvfile="../Producers/python/samples_"
suffixcsvfile="_MC.csv"
haddOpt="-k -f"  # "-j $(grep processor /proc/cpuinfo | wc -l)"
rmOpt="-r -f"
analyzer=VVXAnalyzer  # VVGammaAnalyzer
years="2016"
nevents=-1

#regions="SR4P CR3P1F CR2P2F SR4P_1L SR4P_1P CR4P_1F SR3P CR110 CR101 CR011 CR100 CR001 CR010 CR000 SR3P_1L SR3P_1P CR3P_1F CRLFR SR2P SR2P_1L SR2P_1P CR2P_1F"
regions="SR4P;CR2P2F;CR3P1F;SR3P;CR000;CR001;CR010;CR011;CR100;CR101;CR110;SR2P"

options="--internal-cross-section --nofr --fpw"


if [ ! -d logdir ] ; then
    echo "logdir created"
    mkdir logdir
fi

# echo "--- Analyses on MC ---"
# echo ""

for year in $years ; do
    csvfile=$prefixcsvfile${year}UL$suffixcsvfile
    mc_samples=$(ls samples/MC/$year | grep -v $year | grep -oP ".+(?=\.root)" | grep -v "WLLGTo2L2j_5f_LO\|ZZTo4l_M1ToInf" | sort)

    for sample in $mc_samples ; do
 	echo "Running $analyzer on $sample in region $regions with $csvfile $nevents $options"
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
