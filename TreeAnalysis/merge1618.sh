#!/bin/bash

##################################################
#  Merges results from all the years             #
#                                                #
#  Author: A. Mecca (alberto.mecca@cern.ch)      #
#  Co-author: E. Racca (eleonora.racca@cern.ch)  #
##################################################

# setting options
haddOpt="-k -f -j $(grep processor /proc/cpuinfo | wc -l)"
rmOpt="-r -f"


cd results


# array of analyses and regions names
exist_an_reg=$(echo $(ls 2016) $(ls 2017) $(ls 2018) | tr ' ' '\n' | sort -u | tr '\n' ' ')
#echo $exist_an_reg


# hadd-ing ggTo4l results and removing previous files
for y in 2016 2017 2018 ; do
	for an_reg in $exist_an_reg ; do
		#echo "+++" $y "+++" $an_reg
		[ -d $y/$an_reg ] || continue
		rm $rmOpt $y/$an_reg/gg.root
		rm $rmOpt $y/$an_reg/WZ.root
		
		hadd $haddOpt $y/$an_reg/gg.root $(find $y/$an_reg -type f -name "ggTo*.root")
	done
done


# array of results file names
exist_sam=$(find ./201*/ -type f -name "*.root" | grep -oP "[^/]+\.root" | grep -oP "[^\.]+(?=\.root)" | grep -v "ggTo" | sort | uniq)
#echo $exist_sam


# removing previous files
rm $rmOpt 1618


# hadd-ing results by year and by analysis
for an_reg in $exist_an_reg ; do
	mkdir -p 1618/$an_reg
 
	for sam in $exist_sam ; do
		hadd $haddOpt 1618/$an_reg/$sam.root 2016/$an_reg/$sam.root 2017/$an_reg/$sam.root 2018/$an_reg/$sam.root
	done
done


cd ..
