#!/bin/bash

#########################################################################
# Hack to have each sample in its own folder when using batch_Condor.py #
#                                                                       #
# Author A. Mecca (alberto.mecca@cern.ch)                               #
#########################################################################

all=$@
outputdir="."
options=""

positional=""
while [ $# -gt 0 ] ; do
  case $1 in
    -o|--output-dir)
      if [ -n "$2" ] && [ ${2:0:1} != "-" ]; then
        outputdir=$2
        shift 2
      else
        echo "Warning: argument needed for $1. Using default: $outputdir"
        shift 1
      fi ;;
    -*|--*) # Stuff that batch_Condor.py will handle
      if [ -n "$2" ] && [ ${2:0:1} != "-" ]; then
        options="$options $1 $2"
        shift 2
      else
        options="$options $1"
        shift 1
      fi ;;
    *) # preserve positional arguments
      positional="$positional $1"
      shift ;;
  esac
done
# set positional arguments in their proper place
#eval set -- "$positional"

echo "positional = $positional"
echo "options = $options"
echo "outputdir = $outputdir"

inCSV=$(echo $positional | cut -d " " -f 1)
[ -z $inCSV ] && echo "Error: CSV file needed (positional arg)" >&2 && exit
positional=$(echo "$positional" | sed "s:\s*$inCSV\s*::")

year=$(echo "$inCSV" | grep -oP "\d{4}" )
CSVcontent=$(grep -P '^\s*[^#$]' $inCSV)
header=$(echo "$CSVcontent" | head -n 1)
lines=$(echo "$CSVcontent" | tail -n +2 | head -n 5)

tempfile=$(mktemp)

echo "$lines" | while read -r line
do
  sample=$(echo "$line" | cut -d "," -f 1)
  printf '%s\n%s\n' "$header" "$line" > $tempfile
  echo "batch_Condor.py $options -o $outputdir/$year/$sample $tempfile $positional"
  #batch_Condor.py $options -o $outputdir/$year/$sample $tempfile $positional
done

rm $tempfile

