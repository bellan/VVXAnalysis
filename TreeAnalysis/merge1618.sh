#!/bin/bash

##################################################
#  Merges results from all the years             #
#                                                #
#  Author: A. Mecca (alberto.mecca@cern.ch)      #
#  Co-author: E. Racca (eleonora.racca@cern.ch)  #
##################################################

haddOpt="-k -j $(grep processor /proc/cpuinfo | wc -l)"

cd results

exist_an_reg=$(echo $(ls 2016) $(ls 2017) $(ls 2018) | uniq)
#echo $exist_an_reg

for y in 2016 2017 2018 ; do
	for an_reg in $exist_an_reg ; do
		#echo "+++" $y "+++" $an_reg
		[ -d $y/$an_reg ] || continue
		hadd $haddOpt $y/$an_reg/gg.root $(find $y/$an_reg -type f -name "ggTo*.root")
	done
done

exist_sam=$(find . -type f -name "*.root" | grep -oP "[^/]+\.root" | grep -oP "[^\.]+(?=\.root)" | grep -v "ggTo" | sort | uniq)
#echo $exist_sam

for an_reg in $exist_an_reg ; do
	mkdir -p 1618/${an_reg}
	for sam in $exist_sam ; do
		hadd $haddOpt 1618/${an_reg}/${sam}.root 2016/${an_reg}/${sam}.root 2017/${an_reg}/${sam}.root 2018/${an_reg}/${sam}.root
	done
done

cd ..
