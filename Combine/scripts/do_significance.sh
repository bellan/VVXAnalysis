#!/bin/sh

[ $# -eq 1 ] || { echo "Usage: $0 DATACARD" >&2 ; exit 2 ; }
[ -f $1 ] || { echo "$1 does not exist or is not a file" >&2 ; exit 2 ; }

combine -M Significance -t -1 --expectSignal=1 $1
