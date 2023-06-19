#!/bin/bash -xu

haddOpt="-k -f"
rmOpt="-r -f"
analyzer=VVGammaAnalyzer  # VVXAnalyzer
years="2018"

regions="SR4P;CR3P1F;CR2P2F;SR3P;CR110;CR101;CR011;CR100;CR010;CR001;CR000;SR2P;CRLFR"
options=""  # "--nofr"  # 

make || exit

# ~~~~~ Analysis on data
echo -e "--- Analyses on Data ---\n"

for year in $years ; do
    data_samples=$(ls samples/Data/$year | grep $year | grep -oP ".+(?=\.root)" | sort)
    for sample in ${data_samples} ; do
        ./python/run.py $analyzer $sample -y $year -r $regions $options -d samples/Data &> logdir/${sample}.log &
    done

    wait

    for region in $(echo $regions | tr ';' ' ') ; do
        folder="results/$year/${analyzer}_${region}"
        (
            cd $folder
            hadd $haddOpt data.root $(printf "%s.root " ${data_samples})  # printf reuses the format string to consume all of its arguments
        )
    done
    
    # Hadd lepton CRs only if they've been regenerated
    if echo $regions | grep -q CR3P1F || echo $regions | grep -q CR2P2F ; then
        hadd -k -f results/$year/${analyzer}_SR4P/fake_leptons.root \
            results/$year/${analyzer}_CR3P1F/data.root \
            results/$year/${analyzer}_CR2P2F/data.root
    fi
    
    if echo $regions | grep -q CR110 || echo $regions | grep -q CR101 || echo $regions | grep -q CR011 || \
       echo $regions | grep -q CR100 || echo $regions | grep -q CR010 || echo $regions | grep -q CR001 || echo $regions | grep -q CR000 ; then
        hadd -k -f results/$year/${analyzer}_SR3P/fake_leptons.root \
            results/$year/${analyzer}_CR110/data.root \
            results/$year/${analyzer}_CR101/data.root \
            results/$year/${analyzer}_CR011/data.root \
            results/$year/${analyzer}_CR100/data.root \
            results/$year/${analyzer}_CR010/data.root \
            results/$year/${analyzer}_CR001/data.root \
            results/$year/${analyzer}_CR000/data.root
    fi
done
