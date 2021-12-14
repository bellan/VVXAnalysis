#!/bin/bash

all=$@
outputdir="."
extra=""

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
        extra="$extra $1 $2"
        shift 2
      else
        extra="$extra $1"
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
echo "extra = $extra"
echo "outputdir = $outputdir"

inCSV=$(echo $positional | cut -d " " -f 1)
[ -z $inCSV ] && echo "Error: CSV file needed (positional arg)" >&2 && exit

CSVcontent=$(grep -P '^\s*[^#$]' $inCSV)
header=$(echo "$CSVcontent" | head -n 1)
lines=$(echo "$CSVcontent" | tail -n +2 | head -n 5)
year=$(echo "$inCSV" | grep -oP "\d{4}" )

tempfile=$(mktemp)

echo "$lines" | while read -r line
do
  sample=$(echo "$line" | cut -d "," -f 1)
  printf '%s\n%s\n' "$header" "$line" > $tempfile
  #echo "batch_Condor.py $extra -o $outputdir/$year/$sample $positional"
  batch_Condor.py $extra -o $outputdir/$year/$sample $positional
done

rm $tempfile

