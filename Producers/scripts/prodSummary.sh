#!/bin/sh

################################################################
#  Creates a readable summary of the output of checkProd.csh.  #
#  It also counts the number of jobs already done              #
#  First argument is the job folder, which defaults to "."     #
#                                                              #
#  Author: A. Mecca (amecca@cern.ch)                           #
################################################################

#printf "Usage: ${0##*/} FOLDER\nFOLDER is the name of a subfolder where the job chunks are\n" && exit 1

print_done () { 
    [ -d ${1}/AAAOK ] && find ${1}/AAAOK -type d -name "*_Chunk*" | wc -l | xargs printf "Done: %d jobs\n"
}

[ $# -ne 1 ] && folder="." || folder=$1


output=$(cd ${folder} && checkProd.csh 2>/dev/null) #~/work/CMSSW_10_2_18/bin/slc7_amd64_gcc700/checkProd.csh
[ $? -ne 0 ] && { print_done $folder ;  exit ; }

[ -f $folder/condor.sub ] && grep "JobFlavour" $folder/condor.sub | xargs echo
keys=$(echo "$output" | cut -d : -f 2 | sort | uniq | sed "s/^ //g")

todo=$(ls $folder | grep Chunk | wc -l)
[ $todo -eq 0 ] && { print_done $folder && exit 0 ; }

printf "TODO: %d jobs\n" $(echo "$output" | wc -l)
echo "$keys" | while read -r key; do
        printf "\t%d --> \t%s\n" $(grep "$key" <<<"$output" | wc -l) "$(sed 's/ \+/ /g' <<<$key)"
done #< <(echo "$keys")

print_done $folder
