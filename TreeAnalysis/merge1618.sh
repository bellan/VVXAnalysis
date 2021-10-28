#!/bin/bash

##################################################
#  Merges results from all the years             #
#                                                #
#  Author: A. Mecca (alberto.mecca@cern.ch)      #
#  Co-author: E. Racca (eleonora.racca@cern.ch)  #
##################################################

# ~~~~~ Setting options
haddOpt="-k -f -j $(grep processor /proc/cpuinfo | wc -l)"
rmOpt="-r -f"


cd results


# ~~~~~ Removing previous files
echo "  " 
echo "--- Removing elder files ---"
echo "  " 

rm $rmOpt 1618

for y in 2016 2017 2018 ; do
	rm $rmOpt $y/ZZjjAnalyzer_CR
done


# array of analyses and regions names
exist_an_reg=$(echo $(ls 2016) $(ls 2017) $(ls 2018) | tr ' ' '\n' | sort -u | tr '\n' ' ')


# changing names of some of the results as the original samples have different names across years despide being the same process
for an_reg in $exist_an_reg ; do
	if test -f 2016/$an_reg/ZZ4lJJ.root ; then
		rm $rmOpt 2016/$an_reg/ZZJJTo4L.root
		mv -f 2016/$an_reg/ZZ4lJJ.root 2016/$an_reg/ZZJJTo4L.root
	fi
	
	if test -f 2018/$an_reg/ZZTo4lext1.root ; then
		rm $rmOpt 2018/$an_reg/ZZTo4l.root
		mv -f 2018/$an_reg/ZZTo4lext1.root 2018/$an_reg/ZZTo4l.root	
	fi
done


# ~~~~~ Hadd-ing results
# hadd-ing ggTo4l results and removing hadd-ed files
echo "  " 
echo "--- Hadd-ing gg -> ZZ results ---"
echo "  " 

for y in 2016 2017 2018 ; do	
	for an_reg in $exist_an_reg ; do
		
		#echo "--- + ---" $y " -> " $an_reg
		[ -d $y/$an_reg ] || continue
		rm $rmOpt $y/$an_reg/gg.root
		rm $rmOpt $y/$an_reg/WZ.root
		
		hadd $haddOpt $y/$an_reg/gg.root $(find $y/$an_reg -type f -name "ggTo*.root")
	done
done

# array of results file names
exist_sam=$(find ./201*/ -type f -name "*.root" | grep -oP "[^/]+\.root" | grep -oP "[^\.]+(?=\.root)" | grep -v "ggTo" | grep -v "Double" | grep -v "Single" | grep -v "Muon" | grep -v "EGamma" | grep -v "20" | sort | uniq)
#echo $exist_sam


# hadd-ing results by year and by analysis
echo "  " 
echo "--- Hadd-ing results into 1618 ---"
echo "  " 

for an_reg in $exist_an_reg ; do
	mkdir -p 1618/$an_reg
 
	for sam in $exist_sam ; do
		hadd $haddOpt 1618/$an_reg/$sam.root 2016/$an_reg/$sam.root 2017/$an_reg/$sam.root 2018/$an_reg/$sam.root
	done
done


# hadd-ing control regions
echo "  " 
echo "--- Creating CR for all years ---"
echo "  " 

for y in 2016 2017 2018 1618 ; do
	mkdir -p $y/ZZjjAnalyzer_CR

	for sam in $exist_sam ; do
		hadd $haddOpt $y/ZZjjAnalyzer_CR/$sam.root $y/ZZjjAnalyzer_CR3P1F/$sam.root $y/ZZjjAnalyzer_CR2P2F/$sam.root
	done
done


cd ..
