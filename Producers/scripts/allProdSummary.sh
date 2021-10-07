#!/bin/sh

###################################################################
#  Creates a readable summary of the output of checkProd.csh for  #
#  each of the subfolders in the folders passed as arguments      #
#                                                                 #
#  Author: A. Mecca (amecca@cern.ch)                              #
###################################################################

[ $# -lt 1 ] && echo "Usage: ${0##*/} FOLDER[S]" && exit 1

for topdir in $@ ; do
    #echo "----- $topdir -----"
    dirs=$(find $topdir -maxdepth 1 -mindepth 1 -type d)
    for d in $dirs ; do
        echo "----- $d -----"
        prodSummary.sh "$d" | sed "s/^/\t/g"
    done
    echo ""
done
