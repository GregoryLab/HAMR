#!/bin/bash

# takes mismatch bed file and counts nucleotide
#   frequencies at each site

INBED=$1

awk '
function init() {
  counts["A"] = 0;
  counts["C"] = 0;
  counts["G"] = 0;
  counts["T"] = 0;
  counts["NR"] = 0;
}
function output() {
  if (counts["NR"]>0)
  print prev_chr, prev_bp, prev_str, prev_from_nuc,
     counts["A"], counts["C"], counts["G"],counts["T"],counts["NR"]
}
BEGIN {
  FS="\t"
  OFS="\t"
  init()
}
{
  loc = $1";"$2";"$6

  if (NR > 1 && (loc != prev_loc)) {
    output()
    init()
  }

  split($4,a,">")
  from_nuc = a[1]
  to_nuc = a[2]

  split($5,b,";")
  count = b[1]

  if (to_nuc == ".")
    to_nuc = from_nuc
  else
    counts["NR"] += count

  counts[to_nuc] += count

  prev_from_nuc = from_nuc
  prev_chr = $1
  prev_bp = $2
  prev_str = $6
  prev_loc = $1";"$2";"$6
}
END { output() }' ${INBED} | sort -k1,1 -k2,2n -k3,3

