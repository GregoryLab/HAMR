#!/bin/python2.7
#
#  Copyright (c) 2016 University of Pennsylvania
#
#  Permission is hereby granted, free of charge, to any person obtaining a
#  copy of this software and associated documentation files (the "Software"),
#  to deal in the Software without restriction, including without limitation
#  the rights to use, copy, modify, merge, publish, distribute, sublicense,
#  and/or sell copies of the Software, and to permit persons to whom the
#  Software is furnished to do so, subject to the following conditions:
#
#  The above copyright notice and this permission notice shall be included in
#  all copies or substantial portions of the Software.
#
#  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
#  OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
#  FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
#  AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
#  LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
#  FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
#  DEALINGS IN THE SOFTWARE.

import argparse
import subprocess
import os
import sys
import datetime
import time
import re
import os.path
import tempfile
import collections
import shutil

###############################################################################################################

parser = argparse.ArgumentParser(
    description="Calls differential modified positions across two HAMR analyses. Assumes stranded data")
parser.add_argument('a_tag', help='Short tag to identify analysis A in output files')
parser.add_argument('b_tag', help='Short tag to identify analysis B in output files')
parser.add_argument('out_file', help='Output file for comparisons')
parser.add_argument('threshold', help='HAMR-accessibility threshold')
parser.add_argument('--a_mods_datasets', '-a', action='store', dest='a_mods_filenames', nargs='+',
                    help='raw.txt files from analysis A')
parser.add_argument('--b_mods_datasets', '-b', action='store', dest='b_mods_filenames', nargs='+',
                    help='raw.txt files from analysis B')
args = parser.parse_args()

# parser.add_argument('--hamr_acc_threshold','-t',action='store', dest='hamr_acc_threshold', nargs='?',
# default='unspecified', help='(Integer) Manually specify the threshold for HAMR accessibility')


print("###Calling differential modified positions...")

# get the date and time
now = datetime.datetime.now()
datelist = [str(now.year), str(now.month), str(now.day), str(now.hour), str(now.minute), str(now.second),
            str(now.microsecond)]
rightnow = "_".join(datelist)

# Make a temporary directory for intermediate files
tmp_folder_name = './tmp_' + rightnow
subprocess.check_call(['mkdir', tmp_folder_name])

# define delimiter for hash keys
delimiter = ";"

##################################################################################################
# Define positions modified in all or any A dataset. Also compile modifications to dictionary of lists.

positionsAllA = collections.defaultdict(int)
positionsAnyA = collections.defaultdict(int)
modsA = collections.defaultdict(list)

for infn in args.a_mods_filenames:
    with open(infn) as f:
        next(f)  # skip header
        for line in f:
            line = line.rstrip().split("\t")
            if line[14] == "TRUE":
                position = line[0] + delimiter + line[1] + delimiter + line[2]
                positionsAllA[position] += 1
                positionsAnyA[position] += 1
    mods_fn = infn.replace('raw', 'mods')
    with open(mods_fn) as f:
        next(f)  # skip header
        for line in f:
            line = line.rstrip().split("\t")
            position = line[0] + delimiter + line[1] + delimiter + line[2]
            mod = line[15]
            modsA[position].append(mod)

# Delete any positions that aren't modified in all datasets
for position in positionsAllA.keys():
    if positionsAllA[position] < len(args.a_mods_filenames):
        del positionsAllA[position]

# Define positions modified in all or any B dataset. Also compile modifications to dictionary of lists.

positionsAllB = collections.defaultdict(int)
positionsAnyB = collections.defaultdict(int)
modsB = collections.defaultdict(list)

for infn in args.b_mods_filenames:
    with open(infn) as f:
        next(f)  # skip header
        for line in f:
            line = line.rstrip().split("\t")
            if line[14] == "TRUE":
                position = line[0] + delimiter + line[1] + delimiter + line[2]
                positionsAllB[position] += 1
                positionsAnyB[position] += 1
    mods_fn = infn.replace('raw', 'mods')
    with open(mods_fn) as f:
        next(f)  # skip header
        for line in f:
            line = line.rstrip().split("\t")
            position = line[0] + delimiter + line[1] + delimiter + line[2]
            mod = line[15]
            modsB[position].append(mod)

# delete any positions that aren't modified in all datasets
for position in positionsAllB.keys():
    if positionsAllB[position] < len(args.a_mods_filenames):
        del positionsAllB[position]

#############################################################################################################

# Define maximum coverage of positions in any A dataset

maxCovA = collections.defaultdict(int)

for infn in args.a_mods_filenames:
    with open(infn) as f:
        next(f)  # skip header
        for line in f:
            line = line.rstrip().split("\t")
            position = line[0] + delimiter + line[1] + delimiter + line[2]
            coverage = int(line[8]) + int(line[9])
            if coverage > maxCovA[position]:
                maxCovA[position] = coverage

# Define maximum coverage of positions in any B dataset

maxCovB = collections.defaultdict(int)

for infn in args.b_mods_filenames:
    with open(infn) as f:
        next(f)  # skip header
        for line in f:
            line = line.rstrip().split("\t")
            position = line[0] + delimiter + line[1] + delimiter + line[2]
            coverage = int(line[8]) + int(line[9])
            if coverage > maxCovB[position]:
                maxCovB[position] = coverage

#########################################################################################################

# Define positions modified in all A but none of B

positionsOnlyA = positionsAllA.copy()

for position in positionsAllA.keys():
    if position in positionsAnyB:
        del positionsAllA[position]

# Define positions modified in all B but none of A

positionsOnlyB = positionsAllB.copy()

for position in positionsAllB.keys():
    if position in positionsAnyA:
        del positionsAllB[position]

# Define positions modified in all A and all B

positionsAllAAllB = positionsAllA.copy()

for position in positionsAllAAllB.keys():
    if position not in positionsAnyB:
        del positionsAllAAllB[position]

# Define positions modified in either all A or all B (for purposes of output file)
positionsToOutput = set(positionsAllA.keys() + positionsAllB.keys())

#############################################################################################

# Compile coverages to a single file
print("finding differential mods...")
outfn = open(args.out_file, 'w')
outfn.write(
    'chr\tstart\tstrand\t' + args.a_tag + '.mods\t' + args.b_tag + '.mods\t' + args.a_tag
    + '.maxCov\t' + args.b_tag + '.maxCov' + '\t' + 'unique_to' + '\t' + 'both_HAMR_acc' + '\n')
total_differential = 0
pass_cov_differential = 0
threshold = int(args.threshold)
for position in positionsToOutput:
    (chromosome, start, strand) = position.split(delimiter)
    coverage_a = (maxCovA[position])
    coverage_b = (maxCovB[position])
    modification_a = ",".join(modsA[position])
    modification_b = ",".join(modsB[position])
    unique = ""
    both_HAMR_acc = ""
    if position in positionsOnlyA:
        unique = args.a_tag
        total_differential += 1
        if maxCovB[position] < threshold:
            both_HAMR_acc = "no"
        else:
            both_HAMR_acc = "yes"
            pass_cov_differential += 1
    if position in positionsOnlyB:
        unique = args.b_tag
        total_differential += 1
        if maxCovA[position] < threshold:
            both_HAMR_acc = "no"
        else:
            both_HAMR_acc = "yes"
            pass_cov_differential += 1
    if position in positionsAllAAllB:
        unique = "both"
        both_HAMR_acc = "yes"
    outfn.write(chromosome + '\t' + start + '\t' + strand + '\t' + modification_a + '\t' + modification_b + '\t' + str(
        coverage_a) + '\t' + str(coverage_b) + '\t' + unique + '\t' + both_HAMR_acc + '\n')
outfn.close()

#####################################################################################################

# remove tmp directory
shutil.rmtree(tmp_folder_name)

# conclusion message
print(
    str(pass_cov_differential) + " differential modifications on bases accessible to both libraries,\n of a total of "
    + str(total_differential)
    + " called without considering HAMR accessibility###\n\n------------------------------------\n")
exit()
