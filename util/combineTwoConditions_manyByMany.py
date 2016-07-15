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

###version 1.0

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
from distutils import spawn

BEDTOOLS=spawn.find_executable("bedtools")
if BEDTOOLS is None:
   print  "***ERROR: bedtools is not found"
   sys.exit("Please install bedtools or make sure it is in the PATH")

RSCRIPT=spawn.find_executable("Rscript")
if RSCRIPT is None:
   print "***ERROR: Rscript is not found"
   sys.exit("Please instal R / Rscript or make sure it is in the PATH")

####################################################################################################################################################################

parser = argparse.ArgumentParser(description="Calls differential modified positions across two HAMR analyses. Assumes stranded data")
parser.add_argument('--a_mods_datasets','-a',action='store', dest='a_mods_filenames', nargs='+', help='bed files from analysis A')
parser.add_argument('--b_mods_datasets','-b',action='store', dest='b_mods_filenames', nargs='+', help='bed files from analysis B')
parser.add_argument('--output', '-o', action='store', dest='output', nargs=1, required=True, help='Output directory')
args=parser.parse_args()

# define delimiter for hash keys
delimiter = ";"

####################################################################################################################################################################
### Define positions modified in all or any A dataset. Also compile modifications to dictionary of lists.

positions = {}
total_files = len(args.a_mods_filenames+args.b_mods_filenames)
mods = {}
weights = {}


i = -1
for infn in args.a_mods_filenames:
    i += 1
    with open(infn) as f:
        next(f) #skip header
        for line in f:
            fields = line.rstrip().split("\t")
            position = fields[0]+delimiter+fields[1]+delimiter+fields[5] #chr;start;strand
            mod = fields[4]
            if not mods.has_key(position):
                mods[position] = ["NA"] * total_files
            mods[position][i] = mod

for infn in args.b_mods_filenames:
    i += 1
    with open(infn) as f:
        next(f) #skip header
        for line in f:
            fields = line.rstrip().split("\t")
            position = fields[0]+delimiter+fields[1]+delimiter+fields[5] #chr;start;strand
            mod = fields[4]
            if not mods.has_key(position):
                mods[position] = ["NA"] * total_files
            mods[position][i] = mod

#assign weights to each element of pie chart
weights = [str(0.5/len(args.a_mods_filenames))]*len(args.a_mods_filenames)+[str(0.5/len(args.b_mods_filenames))]*len(args.b_mods_filenames) #all A samples will occupy half of pie chart
weights_output =  ",".join(weights)

#print to temporary bed6 file
combinedModsTemp = str(args.output[0])
combinedModsTempOut = open(combinedModsTemp,'w')
for position in sorted(mods):
    chromosome, start, strand = position.split(";")
    mods_output = ",".join(mods[position])
    stop = str(int(start)+1)
    combinedModsTempOut.write(chromosome+"\t"+start+"\t"+stop+"\t"+weights_output+"\t"+mods_output+"\t"+strand+"\n")

