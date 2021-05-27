#!/bin/bash

##################################################
#  Merges results from all the years             #
#                                                #
#  Author: A. Mecca (alberto.mecca@cern.ch)      #
#  Co-author: E. Racca (eleonora.racca@cern.ch)  #
##################################################


exist_an_reg=$(echo $(ls results/2016) $(ls results/2017) $(ls results/2018) | uniq)
#echo $exist_an_reg

exist_sam=$(find results -type f -name "*.root" | grep -oP "[^/]+\.root" | grep -oP "[^\.]+(?=\.root)" | sort | uniq)
#echo $exist_sam

for an_reg in $exist_an_reg ; do
	mkdir -p results/1618/${an_reg}
	for sam in $exist_sam ; do
		hadd -k -f results/1618/${an_reg}/${sam}.root results/2016/${an_reg}/${sam}.root results/2017/${an_reg}/${sam}.root results/2018/${an_reg}/${sam}.root
	done
done
