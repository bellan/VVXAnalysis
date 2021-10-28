#!/bin/bash

##################################################
#  Executes the analysis on data samples         #
#                                                #
#  Author: E. Racca (eleonora.racca@cern.ch)     #
##################################################


# ~~~~~ Setting options
prefixcsvfile="-c ../Producers/python/samples_"
suffixcsvfile="_Data.csv"
haddOpt="-k -f -j $(grep processor /proc/cpuinfo | wc -l)"
rmOpt="-r -f"


# ~~~~~ Analysis on data
echo "--- Analyses on Data ---"
echo "  " 

for year in 2016 2017 2018 ; do
		# csv file
		csvfile=$prefixcsvfile$year$suffixcsvfile
		
 		for reg in SR CR2P2F CR3P1F ; do    				
				python/run.py ZZjjAnalyzer data -y $year -r $reg $csvfile
 		done
    
 		python/run.py WZAnalyzer data -y $year -r SR3L $csvfile
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
