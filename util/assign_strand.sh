#!/bin/bash

# takes HAMR output and a gene annotation BED (preferably exonic)
# and assigns strand based on transcripts at that site
# NOTE: discards site at sites of bidirectional transcription
#       and sites that don't fall within the gene BED intervals

if [ $# -lt 2 ]; then
    echo "USAGE: $0 hamr_file transcripts.bed" >&2
    exit 1
fi

hamr_file=$1
tx_bed=$2

cat $hamr_file | ./hamr2bed.awk | \
  intersectBed -wao -a - -b $tx_bed  | \
  awk '{FS="\t"; OFS="\t"}
       NR==FNR {
         if ($12 == "+") {x[$4] = or(x[$4], 1)}
         if ($12 == "-") {x[$4] = or(x[$4], 2)}
         if ($12 == ".") {x[$4] = 0 }
         next
       }
       FNR == 1 { print; next }
       {
         k = $1";"$2
         if (x[k] == 1) { $3 = "+" }
         else if (x[k] == 2) { $3 = "-" }
         else { next }
         print
       }
       ' - $hamr_file | \
awk 'BEGIN {
       FS="\t"; OFS="\t"
       comp["A"]="T"; comp["C"]="G"; comp["G"]="C"; comp["T"]="A";
     }
     FNR==1 {print; next}
     $3 == "+" {print; next}
     $3 == "-" {
       $3 = "-"
       $4 = comp[$4]
       tmp = $5; $5 = $8; $8 = tmp
       tmp = $6; $6 = $7; $7 = tmp
       print 
     }'

