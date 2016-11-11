#!/usr/local/bin/python
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

#version 


import sys
import os
import argparse
import subprocess
import re
import random
import math
import multiprocessing
import datetime
import collections
from distutils import spawn


#check python version
if sys.hexversion < 0x020700F0:
    print "Detected Python " + sys.version
    sys.exit("***ERROR: Must be using Python 2.7.x (recommended) or above")


#global binaries not standard in unix environment
BEDTOOLS=spawn.find_executable("bedtools")
if BEDTOOLS is None:
   print  "***ERROR: bedtools is not found"
   sys.exit("Please install bedtools or make sure it is in the PATH")

RSCRIPT=spawn.find_executable("Rscript")
if RSCRIPT is None:
   print "***ERROR: Rscript is not found"
   sys.exit("Please instal R / Rscript or make sure it is in the PATH")

#local scripts not standard in unix environment
hamr_util_dir=os.path.dirname(os.path.realpath(sys.argv[0]))
gff3_overlap_plot=hamr_util_dir+"/"+"gff3_overlap_plot.R"


#command line arguments
parser = argparse.ArgumentParser(description='Overlaps bed feature file with gff3 annotations, and outputs summary')
parser.add_argument('bed', help='bed file of features (e.g. modifications)')
parser.add_argument('gff3', help='gff3 file of features of interest')
parser.add_argument('config',help='configuration file with comma-separated list of features to tabulate')
parser.add_argument('--tag','-t', action='store', help='prefix tag for all output files')
parser.add_argument('--out_dir','-o', action='store', help="specify output dir")
args = parser.parse_args()
output_folder = str(args.out_dir)
output_prefix = str(args.tag)

#check features bed format (i.e. bed6, bed9, or bed12)
bed_colnum = subprocess.check_output(
                      ['awk', '{print NF; exit}',args.bed] )
bed_colnum = int(bed_colnum)
print 'bed'+str(bed_colnum)+' target file detected'


#Check for output directory and make it if neccessary
output_folder = re.sub('\/$', '', output_folder)
if  os.path.isdir(output_folder):
    print "Existing output folder " + output_folder + " detected, will overwrite files if they are already present"
subprocess.check_call(['mkdir', '-p', output_folder])


# initialize dictionary containing lists of gff3 features assigned to each bed feature
features_dict = collections.defaultdict(set)

# overlap bed and gff3, and parse through output
overlap_ps = subprocess.Popen([BEDTOOLS, 'intersect', '-s', '-wo', '-a', args.bed, '-b', args.gff3], stdout=subprocess.PIPE)
while True:
  line = overlap_ps.stdout.readline()
  if line != '':
    line_elements = line.rstrip().split("\t")
    position = line_elements[0]+'_'+line_elements[1]+'_'+line_elements[5]
    feature = line_elements[bed_colnum+2]
    features_dict[position].add(feature)
  else:
    break

# initialize dictionary of features of interest, define by config file
config_FH=open(args.config,'r')
featuresOfInterest=config_FH.readline().rstrip().split(",")
featuresOfInterestPlusAltIntrons = list(featuresOfInterest)
featuresOfInterestPlusAltIntrons.append("alt_intron")
featuresOfInterestPlusAltIntrons.append("multiple")
print featuresOfInterestPlusAltIntrons
featureCounts = collections.OrderedDict((key,int(0)) for key in featuresOfInterestPlusAltIntrons)



# loop through features dict and tabulate overlap. Output to summary file
summary_file=output_folder+'/'+output_prefix+".summary.txt"
summary_FH=open(summary_file,'w')
for position in features_dict.keys():
	currentFeatures = features_dict[position]
	currentFeaturesOfInterest = currentFeatures.intersection(set(featuresOfInterest))
	currentFeaturesOfInterestNoIntrons = currentFeaturesOfInterest - set(["intron"])
	if (len(currentFeaturesOfInterest)==1):
		featureCounts[list(currentFeaturesOfInterest)[0]] += 1
	elif (len(currentFeaturesOfInterestNoIntrons)>1):
		featureCounts["multiple"] += 1
	else:
		featureCounts["alt_intron"] += 1
	summary_FH.write(position+"\t"+",".join(currentFeaturesOfInterest)+"\n")
summary_FH.close()


# print histogram
histogram_file=output_folder+'/'+output_prefix+".histogram.txt"
histogram_FH=open(histogram_file,'w')

for feature,count in featureCounts.items():
	histogram_FH.write(feature+"\t"+str(count)+"\n")

histogram_FH.close()

# plot results in R
plot_file=output_folder+'/'+output_prefix+".plot.pdf"
subprocess.Popen([RSCRIPT, gff3_overlap_plot, histogram_file, plot_file])





	