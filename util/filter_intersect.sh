#!/bin/bash

if [ $# -lt 2 ]; then
    echo "USAGE: $0 hamr filter_bed" >&2
    exit 1
fi

# remove HAMR sites intersecting with sites in a BED file
./hamr2bed.awk $1 | \
  intersectBed -wb -a $2 -b - | \
  awk 'NR == FNR {filt[$7";"$8] = 1; next; }
       FNR==1 || !filt[$1";"$2] { print }' \
     - $1


