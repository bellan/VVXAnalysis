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
    [ -d ${1}/AAAOK ] && find -L ${1}/AAAOK -maxdepth 1 -type d -name "*_Chunk*" | wc -l | xargs printf "%d: DONE\n"
}

[ $# -ne 1 ] && folder="." || folder=$1


output=$(cd ${folder} && checkProd.csh 2>/dev/null) #~/work/CMSSW_10_2_18/bin/slc7_amd64_gcc700/checkProd.csh
[ $? -ne 0 ] && { print_done $folder ;  exit ; }

[ -f $folder/condor.sub ] && grep --color=never "JobFlavour" $folder/condor.sub
keys=$(echo "$output" | cut -d : -f 2 | sed "s/^ //g" | sort | uniq)

todo=$(ls $folder | grep Chunk | wc -l)
[ $todo -eq 0 ] && { print_done $folder && exit 0 ; }

printf "%d: TODO\n" $todo
[ -z "$output" ] && printf "\t%4d: %s\n" $todo "No info on jobs" && exit 0

echo "$keys" | while read -r key; do
        printf "\t%d:\t%s\n" $(grep "$key" <<<"$output" | wc -l) "$(sed 's/ \+/ /g' <<<$key)"
done #< <(echo "$keys")

print_done $folder
