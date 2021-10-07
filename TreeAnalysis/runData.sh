#!/bin/bash

##################################################
#  Executes the analysis on data samples         #
#                                                #
#  Author: E. Racca (eleonora.racca@cern.ch)     #
##################################################


# ~~~~~ Setting options
csvfile="-c ../Producers/python/samples_2016_Data.csv"
haddOpt="-k -f -j $(grep processor /proc/cpuinfo | wc -l)"
rmOpt="-r -f"


# ~~~~~ Existing samples
samples=$(find ./samples/2016/ -type f \( -name "Double*.root" -or -name "Single*.root" -or -name "Muon*.root" \) | grep -oP "[^/]+\.root" | grep -oP "[^\.]+(?=\.root)" | sort | uniq)


# ~~~~~ Analysis on data
echo "--- Analysis on Data ---"
echo "  " 

for sam in $samples ; do
    for reg in SR CR2P2F CR3P1F ; do
	python/run.py ZZjjAnalyzer $sam -r $reg $csvfile
    done
    
    python/run.py WZAnalyzer $sam -r SR3L $csvfile
done


# ~~~~~ Hadd-ing results
echo "  " 
echo "--- Hadd-ing data results ---"
echo "  " 

cd results/2016

# array of analyses and regions names
an_reg=$(echo $(ls) | tr ' ' '\n' | sort -u | tr '\n' ' ')

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

cd ../..
