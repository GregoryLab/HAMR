#!/usr/local/bin/python
# Converts first n and last m nucleotides Q scores to 0 (ASCII+33), such that HAMR ignores them
# version 1.0

import sys
import pysam
import argparse

# Command line arguments
parser = argparse.ArgumentParser(description='Overlaps bed feature file with gff3 annotations, and outputs summary')
parser.add_argument('input', help='input BAM file')
parser.add_argument('output', help='output BAM file')
parser.add_argument('--five_prime', '-5p', action='store', type=int, help='number of bases to ignore at 5` end of read')
parser.add_argument('--three_prime', '-3p', action='store', type=int, help='number of bases to '
                                                                           'ignore at 3` end of read')
args = parser.parse_args()

########

infile = pysam.AlignmentFile(args.input, "rb")
outfile = pysam.AlignmentFile(args.output, "wb", template=infile)
for read in infile.fetch():
    read_length = len(read.query_qualities)
    new_query_qualities = read.query_qualities
    for i in range(0, args.five_prime):
        new_query_qualities[i] = 0
    for j in range(read_length - args.three_prime, read_length):
        new_query_qualities[j] = 0
    read.query_qualities = new_query_qualities
    outfile.write(read)

infile.close()
outfile.close()
