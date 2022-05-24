#!/bin/bash


prefixcsvfile="../Producers/python/samples_"
suffixcsvfile="_Data.csv"
haddOpt="-k -f"  # "-j $(grep processor /proc/cpuinfo | wc -l)"
rmOpt="-r -f"
analyzer=VVXAnalyzer  # VVGammaAnalyzer
years="2016"

regions="SR4P;CR2P2F;CR3P1F;SR3P;CR000;CR001;CR010;CR011;CR100;CR101;CR110;SR2P"
options="--nofr --fpw"

eras="Bver1 Bver2 Chipm Dhipm Ehipm Fhipm F G H"

make || exit

# ~~~~~ Analysis on data
echo "--- Analyses on Data ---"
echo "  "

for year in $years ; do
    csvfile=$prefixcsvfile${year}UL$suffixcsvfile

    for era in $eras ; do			
	echo $analyzer -y $year era $era -r $regions -c $csvfile $options -d samples/Data
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
done


: '
for year in 2016 2017 2018 ; do
		# samples and csv file
		samples=$(find ./samples/$year/ -type f \( -name "Double*.root" -or -name "Single*.root" -or -name "Muon*.root" -or -name "EGamma*.root" \) | grep -oP "[^/]+\.root" | grep -oP "[^\.]+(?=\.root)" | sort | uniq)
		csvfile=$prefixcsvfile$year$suffixcsvfile
				
		for sam in $samples ; do
    		for reg in SR CR2P2F CR3P1F ; do    				
						echo "python/run.py ZZjjAnalyzer $sam -y $year -r $reg $csvfile"
    		done
    
    		echo "python/run.py WZAnalyzer $sam -y $year -r SR3L $csvfile"
		done
done


# ~~~~~ Hadd-ing results
echo "  " 
echo "--- Hadd-ing data results ---"
echo "  " 

cd results

for year in 2016 2017 2018 ; do
		cd $y
		
		# array of samples and analyses_regions
		samples=$(find ./samples/$year/ -type f \( -name "Double*.root" -or -name "Single*.root" -or -name "Muon*.root" -or -name "EGamma*.root" \) | grep -oP "[^/]+\.root" | grep -oP "[^\.]+(?=\.root)" | sort | uniq)
		an_reg=$(echo $(ls -I "*_CR") | tr ' ' '\n' | sort -u | tr '\n' ' ')

		for an in $an_reg ; do
    		# remove prefious file
    		rm $rmOpt $an/data.root
    
    		# create path to hadd the updated results
    		for sam in $samples ; do
						path+=$(echo $an/$sam.root )
						path+=$(echo " " )
    		done
    		#echo $path
    
    		# new hadd-ed file
    		hadd $haddOpt $an/data.root $path
    
    		# restore the path variable
    		path=		
    
    		# remove hadd-ed dataset
    		for sam in $samples ; do
						rm $rmOpt $an/$sam.root
    		done
    
    		echo "  "
    		echo "--- + ---"
    		echo "  "
		done
done

cd ../..
'
