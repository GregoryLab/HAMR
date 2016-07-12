#!/bin/python2.7
#
#  Copyright (c) 2013 University of Pennsylvania
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

####################################################################################################################################################################

parser = argparse.ArgumentParser(description="Calls differential modified positions across two HAMR analyses. Assumes stranded data")
parser.add_argument('a_mods',help='A raw.txt file from analysis A')
parser.add_argument('b_mods',help='A raw.txt file from analysis B')
#parser.add_argument('a_bam',help='Mapped reads (BAM) file from analysis A. Must be sorted by position')
#parser.add_argument('b_bam',help='Mapped reads (BAM) file from analysis B. Must be sorted by position')
parser.add_argument('a_tag',help='Short tag to identify analysis A in output files')
parser.add_argument('b_tag',help='Short tag to identify analysis B in output files')
parser.add_argument('out_file',help='Output file for comparisons')
parser.add_argument('threshold',help='HAMR-accessibility threshold')
args=parser.parse_args()

print "###Calling differential modified positions..."
 
#get the date and time
now = datetime.datetime.now()
datelist = [str(now.year),str(now.month),str(now.day),str(now.hour),str(now.minute),str(now.second),str(now.microsecond)]
rightnow= "_".join(datelist)

# Make a temporary directory for intermediate files
tmp_folder_name = './tmp_'+rightnow
subprocess.check_call(['mkdir', tmp_folder_name])

####################################################################################################################################################################

# Convert HAMR output (.txt) to .bed
print "converting input to bed format..."
mods_file_a=(args.a_mods).replace("raw", "mods")
mods_file_b=(args.b_mods).replace("raw", "mods")
bed_file_a=tmp_folder_name+'/'+os.path.basename(mods_file_a).replace("txt", "bed")
bed_file_b=tmp_folder_name+'/'+os.path.basename(mods_file_b).replace("txt", "bed")
outfn=open(bed_file_a,'w')
subprocess.check_call(['awk', 'FNR > 1 {print $1"\t"$2"\t"(1+$2)"\t"$6";"$2"\t"$16"\t"$3}', mods_file_a],stdout=outfn)
outfn.close()
outfn=open(bed_file_b,'w')
subprocess.check_call(['awk', 'FNR > 1 {print $1"\t"$2"\t"(1+$2)"\t"$1";"$2"\t"$16"\t"$3}', mods_file_b],stdout=outfn)
outfn.close()

# Create a master list of mods in .bed format
print 'Merging mods to master list...'
merged_bed_file = tmp_folder_name+'/'+'merged.bed'
outfn=open(merged_bed_file,'w')
ps1 = subprocess.Popen(['cat', bed_file_a, bed_file_b], stdout=subprocess.PIPE)
ps2 = subprocess.Popen(['awk', 'FNR > 1 {print $1"\t"$2"\t"$3"\t"null"\t"null"\t"$6}'], stdin=ps1.stdout, stdout=subprocess.PIPE)
ps1.stdout.close()
#ps3 = subprocess.Popen(['sort', '-k1,1', '-k2,2n'], stdin=ps2.stdout, stdout=subprocess.PIPE)
ps3 = subprocess.Popen(['sortBed'], stdin=ps2.stdout, stdout=subprocess.PIPE)
ps2.stdout.close()
ps4 = subprocess.check_call(['uniq'], stdin=ps3.stdout, stdout=outfn)
ps3.stdout.close()

# From merged bed, define an indicator hash of all modified positions
positions = collections.OrderedDict()
infn=open(merged_bed_file,'r')
for line in infn:
    line = line.rstrip().split("\t")
    position = line[0]+'_'+line[1]+'_'+line[5]
    positions[position] = ""
infn.close()

####################################################################################################################################################################

# Grab modification identity at each modified position. If not modified, call "NA"
print "tabulating modification identity..."
mod_a = dict((position,"NA") for position in positions)
mod_b = dict((position,"NA") for position in positions)
infn=open(bed_file_a,'r')
for line in infn:
    line = line.rstrip().split("\t")
    position = line[0]+'_'+line[1]+'_'+line[5]
    mod = line[4]
    mod_a[position] = mod
infn.close()
#
infn=open(bed_file_b,'r')
for line in infn:
    line = line.rstrip().split("\t")
    position = line[0]+'_'+line[1]+'_'+line[5]
    mod = line[4]
    mod_b[position] = mod
infn.close()

####################################################################################################################################################################

# Grab coverages at modified positions from both raw HAMR outputs
print "tabulating coverages..."
cov_a = dict((position,"NA") for position in positions)
cov_b = dict((position,"NA") for position in positions)

with open(args.a_mods) as infn:
    next(infn) #skip header
    for line in infn:
        line = line.rstrip().split("\t")
        position = line[0]+'_'+line[1]+'_'+line[2]
        if position in positions:
            truefalse = line[14]
            coverage = int(line[8])+int(line[9])
            cov_a[position] = coverage
            #in_a[position] = truefalse
infn.close()
#
with open(args.b_mods) as infn:
    next(infn) #skip header
    for line in infn:
        line = line.rstrip().split("\t")
        position = line[0]+'_'+line[1]+'_'+line[2]
        if position in positions:
            truefalse = line[14]
            coverage = int(line[8])+int(line[9])
            cov_b[position] = coverage
            #in_b[position] = truefalse
infn.close()

####################################################################################################################################################################

# Compile coverages to a single file
print "finding differential mods..."
outfn=open(args.out_file,'w')
outfn.write('chr\tstart\tstrand\t'+args.a_tag+'.mod\t'+args.b_tag+'.mod\t'+args.a_tag+'.cov\t'+args.b_tag+'.cov'+'\t'+'unique_to'+'\t'+'both_HAMR_acc'+'\n')
total_differential = 0
pass_cov_differential = 0
threshold = int(args.threshold)
for position in positions:
    (chromosome, start, strand) = position.split('_')
    coverage_a = (cov_a[position])
    coverage_b = (cov_b[position])
    modification_a = mod_a[position]
    modification_b = mod_b[position]
    unique = ""
    both_HAMR_acc = ""
    if (modification_a is not "NA" and modification_b is "NA"):
        unique = args.a_tag
        total_differential += 1
    elif (modification_a is "NA" and modification_b is not "NA"):
        unique = args.b_tag
        total_differential += 1
    else:
        unique = "both"
    if (coverage_a is "NA" or coverage_b is "NA"):
        both_HAMR_acc = "no"
    elif (coverage_a < threshold or coverage_b < threshold):
        both_HAMR_acc = "no"
    else:
        both_HAMR_acc = "yes"
        pass_cov_differential += 1
    outfn.write(chromosome+'\t'+start+'\t'+strand+'\t'+modification_a+'\t'+modification_b+'\t'+str(coverage_a)+'\t'+str(coverage_b)+'\t'+unique+'\t'+both_HAMR_acc+'\n')
outfn.close()

####################################################################################################################################################################

#remove tmp directory
shutil.rmtree(tmp_folder_name)

#conclusion message
print str(pass_cov_differential)+" differential modifications on bases accessible to both libraries,\n of a total of "+str(total_differential)+" called without considering HAMR accessibility###\n\n------------------------------------\n"
exit()