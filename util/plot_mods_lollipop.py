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

### version 1.1
#updates since version 1.0: added pie chart functionality. See -pie optional flag

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


# command line arguments
parser = argparse.ArgumentParser(description="Plots modifications along features, using lollipop-style plots")
parser.add_argument('--mods','-m',action='store', dest='mods', nargs=1, required=True, help='bed (6-column) of modifications')
parser.add_argument('--bed12','-b',action='store', dest='bed12', nargs='+', required=True, help='transcriptome bed12. NOTE: Assumes column 4 is transcript id, and column 5 gene id')
parser.add_argument('--output_dir', '-o', action='store', dest='output_dir', nargs=1, required=True, help='Output directory')
parser.add_argument('--output_prefix', '-p', action='store', dest='output_prefix', nargs=1, required=True, help='Prefix for output files')
parser.add_argument('--multi_sample_pie', '-pie', action='store_true', dest='pie', help='print pie charts over each mod depending on which samples mod is observed. NOTE: requires mods.bed output from combineTwoConditions_manyByMany.py. Mods list will be plotted counterclockwise starting from 12 o`clock. Thus, sample 2 replicates should usually be in reverse order')
parser.add_argument('--separate_isoforms', '-s', action='store_true', dest='boundaries', help='print each isoform to a different page')
args=parser.parse_args()
output_folder = str(args.output_dir[0])
output_prefix = str(args.output_prefix[0])
mods_file = str(args.mods[0])

#locations of called scripts
hamr_dir=os.path.dirname(os.path.realpath(sys.argv[0]))
plot=hamr_dir+"/"+"plot_mods_lollipop.R" #R-script
plotPie=hamr_dir+"/"+"plot_mods_lollipop.twoSamplePie.R" #R-script

#Check for output directory and make it if neccessary
if  os.path.isdir(output_folder): #if no out dir then make one
    print "Existing output folder " + output_folder + " detected, will overwrite all internal files"
subprocess.check_call(['mkdir', '-p', output_folder])

#complies regex
temp = re


####################################################################################################################################################################
### Run through input
print "###Parsing bed12..."

infn=args.bed12[0]
output = output_folder+"/names.tmp"
outputOut = open(output, 'w')
with open(infn) as f:
    #initiate using first line
    firstline = f.readline()
    fields = firstline.rstrip().split("\t")
    current_gene = fields[4]
    bed12temp = output_folder+'/'+current_gene+".bed12.tmp"
    bed12tempOut = open(bed12temp,'w')
    bed12tempOut.write(firstline)
    #run through remaining lines
    for line in f:
        fields = line.rstrip().split("\t")
        gene = fields[4]
        #print to same bed12 if gene is the same as previous line
        if gene == current_gene:
            bed12tempOut.write(line)
        #one new gene in reached, 1) intersectBed to grab mods, 2) plot if any mods found, and 3) intiate new temp files
        else:
            #grab mods for previous gene
            bed12tempOut.close()
            modsTemp = output_folder+'/'+current_gene+".mods.tmp"
            modsTempOut = open(modsTemp,'w')
            subprocess.check_call([BEDTOOLS, 'intersect','-s','-u','-a',mods_file,'-b',bed12temp], stdout=modsTempOut)
            modsTempOut.close()
            #plot, only if mods lie in transcript of interest
            if os.stat(modsTemp).st_size != 0:
                pdf = output_folder+'/'+output_prefix+"."+current_gene+".pdf"
                outputOut.write(modsTemp+"\t"+bed12temp+"\t"+pdf+"\n")
            else:
                os.remove(modsTemp)
                os.remove(bed12temp)
            #intialize new temp file
            current_gene = str(gene)
            bed12temp = output_folder+'/'+current_gene+".bed12.tmp"
            bed12tempOut = open(bed12temp,'w')
            bed12tempOut.write(line)

outputOut.close()

###Plot each file in R
print "Generating plots..."
if args.pie:
    subprocess.check_call([RSCRIPT, plotPie, output])
else:
    subprocess.check_call([RSCRIPT, plot, output])
sys.exit()

####################################################################################################################################################################




