#!/bin/sh

#######################################################
#  Modifies JobFlavour in the condor.sub file in the  #
#  folder passed as first argument                    #
#                                                     #
#  Author: A. Mecca (amecca@cern.ch)                  #
#######################################################

[ $# -ne 1 ] && echo "Usage: ${0##*/} NEW_FLAVOUR" && exit 1

flavours="espresso microcentury longlunch workday tomorrow testmatch nextweek"
new_fl=$1
echo $flavours | grep -q $new_fl || { echo "NEW_FLAVOUR must be one of:" $(echo $flavours | sed "s/ /, /g") 1>&2 ; exit 1 ; }

sed -i "s/\(JobFlavour\s\+\= \"\)[^\"]\+/\1$new_fl/" condor.sub
grep "JobFlavour" condor.sub | grep -P --color=always "(?<=\")[^\"]+(?=\")"
