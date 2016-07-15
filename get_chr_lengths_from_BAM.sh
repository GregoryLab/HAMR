#!/bin/bash

#accept either argument as input or STDIN as input
[ $# -ge 1 -a -f "$1" ] && input="$1" || input="-"

#grab lengths file
samtools view -H $input | grep "@SQ" | sed -e 's/@SQ\t//' | sed -e 's/SN://' | sed -e 's/LN://' | sed -e 's/\s+/\t/' 

