#!/usr/local/bin/python
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

import sys

if sys.hexversion < 0x020700F0:
    print "Detected Python " + sys.version
    sys.exit("***ERROR: Must be using Python 2.7.x (recommended)")

import argparse
import subprocess
import os
import datetime
import time
import re
import os.path
from distutils import spawn

RSCRIPT=spawn.find_executable("Rscript")
if RSCRIPT is None:
   print "***ERROR: Rscript is not found"
   sys.exit("Please instal R / Rscript or make sure it is in the PATH")

SAMTOOLS=spawn.find_executable("samtools")
if SAMTOOLS is None:
   print  "***ERROR: samtools is not found"
   sys.exit("Please install samtools or make sure it is in the PATH")

# command line arguments
parser = argparse.ArgumentParser(description="Takes a bam file that has been sorted with redundant reads removed and generates a HAMR predicted_mods.txt output")
parser.add_argument('bam',help='A sorted bam file consisting of nonredundant reads')
parser.add_argument('genome_fas',help='Genome fasta file; WARNING: remember to index the reference using samtools faifx')
parser.add_argument('prediction_training_set',help='modification identity training set model file; .RData format')
parser.add_argument('output_folder',help='name of folder to put HAMR output')
parser.add_argument('out_prefix',help='Prefix for HAMR output')
parser.add_argument('min_qual',help='The minimum quality score of a read to be analyzed')
parser.add_argument('min_cov',help='The minimum coverage of a nucleotide to be analyzed')
parser.add_argument('seq_err',help='The percentage of mismatches based solely on sequencing error')
parser.add_argument('hypothesis',help='The hypothesis to be tested, either "H1" or "H4"')
parser.add_argument('max_p',help='The maximum p-value cutoff')
parser.add_argument('max_fdr',help='The maximum FDR cutoff')
parser.add_argument('refpercent',help='The percentage of reads that must match the reference nucleotide')
parser.add_argument('--target_bed', '-n', action='store', dest='target_bed', nargs='?', default='unspecified', help='Specifies genomic intervals for analysis; e.g. all mRNAs. If unspecified, defaults to whole genome')
parser.add_argument('--paired_ends','-pe',action='store_true',help='Use this tag to indicate paired-end sequencing')
parser.add_argument('--filter_ends','-fe',action='store_true',help='Exclude the first and last nucleotides of a read from the analysis')

args=parser.parse_args()


#Raise error if hypothesis has invalid value
if args.hypothesis != 'H1' and args.hypothesis != 'H4':
    raise ValueError('Hypothesis must be H1 or H4.')
 
#locations of C, Bash, and R scripts
hamr_dir=os.path.dirname(os.path.realpath(sys.argv[0]))
rnapileup=hamr_dir+"/"+"rnapileup" #C-script
filter_pileup=hamr_dir+"/"+"filter_pileup" #C-script
rnapileup2mismatchbed=hamr_dir+"/"+"rnapileup2mismatchbed" #C-script
mismatchbed2table=hamr_dir+"/"+"mismatchbed2table.sh" #Shell script
detect_mods_definite=hamr_dir+"/"+"detect_mods_definite.R" #R script
classify_mods=hamr_dir+"/"+"classify_mods.R" #Rscript

#get flags
pairedends=""
if (args.paired_ends):
    pairedends="--paired"

#Check for output directory and make it if neccessary
output_folder = re.sub('\/$', '', args.output_folder)
subprocess.check_call(['mkdir', '-p', output_folder])
if  os.path.isdir(args.output_folder): #if no out dir then make one
    print "existing output folder detected, will overwrite internal files unless program is broken..."

# make tmp directory if necessary
tmpDIR=output_folder + '/tmp'
subprocess.check_call(['mkdir', '-p', tmpDIR])



#get the date and time
now = datetime.datetime.now()
datelist = [str(now.year),str(now.month),str(now.day),str(now.hour),str(now.minute),str(now.second),str(now.microsecond)]
rightnow= "_".join(datelist)
rTag=tmpDIR + '/' + rightnow + '.mismatches' #date included in file

run_mode = "genome-wide"
if (args.target_bed != 'unspecified'):
   run_mode = 'targeted'
inputBAM=args.bam
print 'Analyzing %s (%s)' %(inputBAM, run_mode)
bamForAnalysis = inputBAM
if (args.target_bed != 'unspecified'):
	# extract alignments for the region(s) of interest
        target_bed = args.target_bed
        print 'Target BED is specified: ' + target_bed
        print 'Restricting BAM to regions in ' + target_bed
        inputBAMbasename=os.path.basename(inputBAM)
        bam_constrained = output_folder + '/' + re.sub('\.[^.]+$','.constrained.bam',inputBAMbasename)    
        fout=open(bam_constrained,'wb')
        subprocess.check_call([SAMTOOLS,'view','-b',inputBAM,'-L',target_bed],stdout=fout)
        fout.close()
        subprocess.check_call([SAMTOOLS,'index',bam_constrained])
        bamForAnalysis=bam_constrained

print "BAM for HAMR analysis: " + bamForAnalysis



print 'Running RNApileup ' + rnapileup
mmbed4=open(rTag+'.bed4','w')
subprocess.check_call([rnapileup,bamForAnalysis,args.genome_fas,pairedends],stdout=mmbed4)
mmbed4.close()

print 'Running filter_pileup...'
mmbed3=open(rTag+'.bed3','w')
input_mmbed3=rTag+'.bed4'
subprocess.check_call([filter_pileup,input_mmbed3,str(args.min_qual),str(int(args.filter_ends))],stdout=mmbed3)
mmbed3.close()

print ("Filter coverage...")
## this will output ALL sites with read depth >= min_cov!!
## this will be the total # of sites for HAMR analysis
mmbed2=open(rTag+'.bed2','w')
input_mmbed2=rTag+'.bed3'
subprocess.check_call(['awk','$4>=' + str(args.min_cov),input_mmbed2],stdout=mmbed2) 
mmbed2.close()



print 'Running rnapileup2mismatchbed...'
# convert pileups into BED file with entry corresponding to the observed (ref nuc) --> (read nucleotide) transitions
mmbed=open(rTag+'.bed','w')
input_mmbed=rTag+'.bed2'
subprocess.check_call([rnapileup2mismatchbed,input_mmbed],stdout=mmbed)
mmbed.close()

print "converting mismatch BED to nucleotide frequency table"
# mismatchbed2table outputs all sites with at least 1 non-ref nuc
final_bed_file=rTag+'.bed'
freq_table=rTag+'.txt'
txt_output=open(freq_table,'w')
subprocess.check_call([mismatchbed2table, final_bed_file],stdout=txt_output)
txt_output.close()

#print "filtering out sites based on non-ref/ref proportions"
# filter by:
#  min ref nuc pct
#  non-ref/ref > 1%
final_freq_table=rTag+'.final.txt'
min_ref_pct=args.refpercent
outf=open(final_freq_table,'w')
#subprocess.check_call(['awk','{cov=$5+$6+$7+$8;nonref=$9; ref=cov-nonref; if (ref/cov>=0.05 && nonref/ref>=0.01) print;}', freq_table],stdout=outf)
subprocess.check_call(['awk','{cov=$5+$6+$7+$8;nonref=$9; ref=cov-nonref; if (ref/cov>='+min_ref_pct+') print;}', freq_table],stdout=outf)
outf.close()

#OUTPUT steps

print "testing for statistical significance..."
last_tmp_file= final_freq_table #rTag+'.txt'
raw_file=output_folder+'/'+args.out_prefix+'.raw.txt'
outfn=open(raw_file,'w')
subprocess.check_call([RSCRIPT,detect_mods_definite,last_tmp_file,args.seq_err,args.hypothesis,args.max_p,args.max_fdr,args.refpercent],stdout=outfn)
outfn.close()

print "predicting modification identity..."
retOut=subprocess.check_output(['grep', '-c','TRUE',raw_file])
true_mods = int(retOut)
prediction_file=output_folder+'/'+args.out_prefix+'.mods.txt'
if (true_mods > 0):
    outfn=open(prediction_file,'w')
    subprocess.check_call([RSCRIPT,classify_mods,raw_file,args.prediction_training_set],stdout=outfn)
    outfn.close()
else:
	sys.exit("No HAMR modifications predicted, output will contain raw table only\nHAMR analysis complete\n\n------------------------------\n")

print "converting output to bed format..."
bed_file=output_folder+'/'+args.out_prefix+".mods.bed"
outfn=open(bed_file,'w')
subprocess.check_call(['awk', 'FNR > 1 {print $1"\t"$2"\t"(1+$2)"\t"$1";"$2"\t"$16"\t"$3}', prediction_file],stdout=outfn)
outfn.close()

threshold = int(args.min_cov)
print "calculating number of HAMR-accessible bases..."
# this is readily available from the filtered by min_cov pileup file
filt_pileup_file=rTag+'.bed2'
retOut=subprocess.check_output(['awk', 'END{print NR}',filt_pileup_file])
HAMR_accessible_bases = int(retOut) 


print "Sites analyzed (read depth>=%d): %d" % (threshold, HAMR_accessible_bases) 
print "Modification sites found: " + str(true_mods) 

mods_per_acc_bases_file=output_folder+'/'+args.out_prefix+".hamr_acc_bases.txt"
outfn=open(mods_per_acc_bases_file,'w')
mods_per_acc_bases = float(true_mods)/float(HAMR_accessible_bases)*1000000
outfn.write('sample\tmods\thamr_accessible_bases\tmods_per_million_accessible_bases\n')
outfn.write(args.out_prefix+'\t'+str(true_mods)+'\t'+str(HAMR_accessible_bases)+'\t'+str(mods_per_acc_bases)+'\n')
outfn.close()

#conclusion message
print "HAMR analysis complete\n\n------------------------------\n"
