#!/bin/bash -x


prefixcsvfile="../Producers/python/samples_"
suffixcsvfile="_Data.csv"
haddOpt="-k -f -j 10"  # "-j $(grep processor /proc/cpuinfo | wc -l)"
rmOpt="-r -f"
analyzer=VVGammaAnalyzer  # VVXAnalyzer
years="2016"

regions="SR4P;CR3P1F;CR2P2F;SR3P;CR110;CR101;CR011;CR100;CR010;CR001;CR000;SR2P"  # "CRLFR"  # 
options=""  # "--nofr"  # 

eras="Bver2 Chipm Dhipm Ehipm Fhipm F G H"  # Bver1

make || exit

# ~~~~~ Analysis on data
echo -e "--- Analyses on Data ---\n"

for year in $years ; do
    csvfile=$prefixcsvfile${year}UL$suffixcsvfile

    for era in $eras ; do			
    	./python/run.py $analyzer 2016$era -y $year -r $regions -c $csvfile $options -d samples/Data &> logdir/data_${year}${era}.log &
    done

    wait

    for region in $(echo $regions | tr ';' ' ') ; do
    	folder="results/$year/${analyzer}_${region}"
    	(
    	    cd $folder
    	    for era in $eras ; do
    		ls | grep -q ".\+${year}${era}.root" && hadd $haddOpt ${year}${era}.root $(ls | grep ".\+${year}${era}.root") || echo "Warning: no files for ${year}${era}.root"
    	    done
    	    hadd $haddOpt data.root $(printf "${year}%s.root " $eras)  # printf reuses the format string to consume all of its arguments
    	)
    done
    
    hadd -k -f results/$year/${analyzer}_SR4P/fake_leptons.root \
	results/$year/${analyzer}_CR3P1F/data.root \
	results/$year/${analyzer}_CR2P2F/data.root
    
    hadd -k -f results/$year/${analyzer}_SR3P/fake_leptons.root \
	results/$year/${analyzer}_CR110/data.root \
	results/$year/${analyzer}_CR101/data.root \
	results/$year/${analyzer}_CR011/data.root \
	results/$year/${analyzer}_CR100/data.root \
	results/$year/${analyzer}_CR010/data.root \
	results/$year/${analyzer}_CR001/data.root \
	results/$year/${analyzer}_CR000/data.root
done
