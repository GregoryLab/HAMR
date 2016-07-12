#!/bin/awk -f

FNR > 1 {print $1"\t"$2"\t"(1+$2)"\t"$1";"$2"\t0\t"$3}
