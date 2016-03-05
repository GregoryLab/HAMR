#!/bin/bash

# call the nucleotide pattern for H0_1-type (SNP-like) sites

if [ $# -lt 1 ]; then
  echo "USAGE: $0 hamr_file"
  exit 1
fi

hamr_file=$1

awk 'BEGIN{ FS="\t"; OFS="\t";
       nucs["A"] = 1; nucs["C"] = 2;
       nucs["G"] = 3; nucs["T"] = 4;
       for(k in nucs) { nucs[nucs[k]] = k } }
     FNR == 1 { print $0"\ttype"; next }
      {
       ref = $4; max_cnt = -1; max_i=-1;
       for(i=1; i<=4; ++i) {
         count = $(i+4)
         if (i != nucs[ref] && count > max_cnt) {
           max_cnt = count;
           max_i = i
         }
       }
       $(1+NF) = (ref)">"(nucs[max_i])
       print
     }' $hamr_file

