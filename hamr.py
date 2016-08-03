#!/usr/local/bin/python
#
#  Copyright (c) 2013-2016 University of Pennsylvania
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

#version 3.2.1 - Bug fix to predicting modification identity... section

import sys

if sys.hexversion < 0x020700F0:
    print "Detected Python " + sys.version
    sys.exit("***ERROR: Must be using Python 2.7.x (recommended) or above")

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

BEDTOOLS=spawn.find_executable("bedtools")
if BEDTOOLS is None:
   print  "***ERROR: coverageBed is not found"
   sys.exit("Please install bedtools 2.2x or make sure it is in the PATH")

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
parser.add_argument('refproportion',help='The proportion of reads that must match the reference nucleotide (from 0 to 1)')
parser.add_argument('--target_bed', '-n', action='store', dest='target_bed', nargs='?', default='unspecified', help='Specifies genomic intervals for analysis; e.g. all mRNAs. If unspecified, defaults to whole genome')
parser.add_argument('--paired_ends','-pe',action='store_true',help='Use this tag to indicate paired-end sequencing')
parser.add_argument('--filter_ends','-fe',action='store_true',help='Exclude the first and last nucleotides of a read from the analysis')
parser.add_argument('--empirical_hamr_acc_threshold','-et',action='store_true',help='Calculate the treshold of HAMR accessibility empirically. Otherwise, assumes it is equal to min_cov (reasonable assumption at 10x or more')
parser.add_argument('--type_plot','-tp',action='store_true',help='Use this tag to include plots of predicted modification type')

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
detect_mods_definite=hamr_dir+"/"+"detect_mods.R" #R script
classify_mods=hamr_dir+"/"+"classify_mods.R" #Rscript
coverageBed_depth_histogram = hamr_dir+"/"+"coverageBed_depth_histogram.pl" #perl script
get_chr_lengths_from_bam = hamr_dir+"/"+"get_chr_lengths_from_BAM.sh" #Shell script
plot_mod_type = hamr_dir+"/"+"summ_mod_type_fromBed.R" #R script

#get flags
pairedends=""
if (args.paired_ends):
    pairedends="--paired"

#Check for output directory and make it if neccessary
output_folder = re.sub('\/$', '', args.output_folder)
if  os.path.isdir(args.output_folder): #if no out dir then make one
    print "Existing output folder " + output_folder + " detected, will overwrite all internal files"
subprocess.check_call(['mkdir', '-p', output_folder])

# make tmp directory if necessary
tmpDIR=output_folder + '/HAMR_temp'
subprocess.check_call(['mkdir', '-p', tmpDIR])

#get the date and time (unqiue string in case multiple HAMR runs output to same temp folder)
now = datetime.datetime.now()
datelist = [str(now.year),str(now.month),str(now.day),str(now.hour),str(now.minute),str(now.second),str(now.microsecond)]
rightnow= "_".join(datelist)
rTag=tmpDIR + '/' + rightnow + '.HAMR.' + args.out_prefix #date included in file
#rTag=tmpDIR + '/' + 'HAMR.' + args.out_prefix

###################################################################################################################################

##run HAMR
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

        #check target bed format (i.e. bed6, bed9, or bed12)
        target_bed_colnum = subprocess.check_output(
                              ['awk', '{print NF; exit}',target_bed] )
        target_bed_colnum = int(target_bed_colnum)
        print 'bed'+str(target_bed_colnum)+' target file detected'

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
rawpileup=rTag+'.pileup.raw'
frawpileup=open(rawpileup,'w')
subprocess.check_call([rnapileup,bamForAnalysis,args.genome_fas,pairedends],stdout=frawpileup)
frawpileup.close()

print 'Running filter_pileup...'
filteredpileup=rTag+'.pileup.filtered'
ffilteredpileup=open(filteredpileup,'w')
subprocess.check_call([filter_pileup,rawpileup,str(args.min_qual),str(int(args.filter_ends))],stdout=ffilteredpileup)
ffilteredpileup.close()

print ("Filter coverage...")
## this will output ALL sites with read depth >= min_cov
## this will be the total # of sites for HAMR analysis
filteredpileupcov=rTag+'.pileup.filtered.'+str(args.min_cov)
ffilteredpileupcov=open(filteredpileupcov,'w')
subprocess.check_call(['awk','$4>=' + str(args.min_cov),filteredpileup],stdout=ffilteredpileupcov) 
ffilteredpileupcov.close()

print 'Running rnapileup2mismatchbed...'
# convert pileups into BED file with entry corresponding to the observed (ref nuc) --> (read nucleotide) transitions
mismatchbed=rTag+'.mismatch.bed'
fmismatchbed=open(mismatchbed,'w')
subprocess.check_call([rnapileup2mismatchbed,filteredpileupcov],stdout=fmismatchbed)
fmismatchbed.close()

print "converting mismatch BED to nucleotide frequency table..."
# mismatchbed2table outputs all sites with at least 1 non-ref nuc
final_bed_file=mismatchbed
freq_table=rTag+'.freqtable.txt'
txt_output=open(freq_table,'w')
subprocess.check_call([mismatchbed2table, final_bed_file],stdout=txt_output)
txt_output.close()

print "filtering out sites based on non-ref/ref proportions..."
final_freq_table=rTag+'.freqtable.final.txt'
min_ref_pct=args.refproportion
outf=open(final_freq_table,'w')
subprocess.check_call(['awk','{cov=$5+$6+$7+$8;nonref=$9; ref=cov-nonref; if (ref/cov>='+min_ref_pct+') print;}', freq_table],stdout=outf)
outf.close()

#OUTPUT steps

print "testing for statistical significance..."
last_tmp_file= final_freq_table
raw_file=output_folder+'/'+args.out_prefix+'.raw.txt'
outfn=open(raw_file,'w')
subprocess.check_call([RSCRIPT,detect_mods_definite,last_tmp_file,args.seq_err,args.hypothesis,args.max_p,args.max_fdr,args.refproportion],stdout=outfn)
outfn.close()

print "predicting modification identity..."
#retOut=subprocess.check_output(['grep', '-c','TRUE',raw_file])
# if retOut is None:
#     true_mods = int(0)
# else:
#     true_mods = int(retOut)
ps1 = subprocess.Popen(('grep', 'TRUE', raw_file), stdout=subprocess.PIPE)
true_mods = subprocess.Popen(('wc', '-l'), stdin=ps1.stdout, stdout=subprocess.PIPE).communicate()[0].rstrip()
ps1.stdout.close()
true_mods=int(true_mods)
print true_mods
prediction_file=output_folder+'/'+args.out_prefix+'.mods.txt'
if (true_mods > 0):
    outfn=open(prediction_file,'w')
    subprocess.check_call([RSCRIPT,classify_mods,raw_file,args.prediction_training_set],stdout=outfn)
    outfn.close()
else:
    #tabulate number of HAMR accessible bases and then exit
    threshold = int(args.min_cov)
    filt_pileup_file=filteredpileupcov
    retOut=subprocess.check_output(['awk', 'END{print NR}',filt_pileup_file])
    HAMR_accessible_bases=int(retOut)
    #print summary file
    mods_per_acc_bases_file=output_folder+'/'+args.out_prefix+".hamr_acc_bases.txt"
    outfn=open(mods_per_acc_bases_file,'w')
    mods_per_acc_bases = 0
    outfn.write('sample\tmods\thamr_accessible_bases\tmods_per_million_accessible_bases\n')
    outfn.write(args.out_prefix+'\t'+str(true_mods)+'\t'+str(HAMR_accessible_bases)+'\t'+str(mods_per_acc_bases)+'\n')
    outfn.close()
    sys.exit("No HAMR modifications predicted, output will contain raw table and number of HAMR accessible bases only\nAssuming min_cov of %s as threshold for HAMR accessibility\nHAMR analysis complete\n\n------------------------------\n" % args.min_cov)

print "converting output to bed format..."
bed_file=output_folder+'/'+args.out_prefix+".mods.bed"
outfn=open(bed_file,'w')
subprocess.check_call(['awk', 'FNR > 1 {print $1"\t"$2"\t"(1+$2)"\t"$1";"$2"\t"$16"\t"$3}', prediction_file],stdout=outfn)
outfn.close()


print "calculating number of HAMR-accessible bases..."

if not args.empirical_hamr_acc_threshold:
    # tablulate number of bases at or above min_cov
    threshold = int(args.min_cov)
    filt_pileup_file=filteredpileupcov
    retOut=subprocess.check_output(['awk', 'END{print NR}',filt_pileup_file])
    HAMR_accessible_bases=int(retOut)
    #num_bases_with_min_cov=int(retOut) 
    #HAMR_accessible_bases=int(num_bases_with_min_cov)

else:
    # calculate threshold for HAMR-accessibility unless manually specified
    print "calculating threshold of HAMR-accessibility..."
    min_cov_file=output_folder+'/'+args.out_prefix+".min_cov.txt"
    outfn=open(min_cov_file,'w')    
    threshold=subprocess.check_output(
         ['awk','{if (NR==1) next; cov=$9+$10; if (NR==2) thres=cov; if (cov < thres) thres=cov; }END{ print thres+0;}', prediction_file])
    threshold=int(threshold)
    print "HAMR accessibility threshold="+str(threshold)
    outfn.write(args.out_prefix+'\t'+str(threshold)+'\n')
    outfn.close()

    # tablulate number of bases at or above threshold
    print "calculating number of HAMR-accessible bases...",
    hamr_acc_bases_file=output_folder+'/'+args.out_prefix+".hamr_acc_bases.txt"
    outfn=open(hamr_acc_bases_file,'w')
    HAMR_accessible_bases=subprocess.check_output(
            ['awk', '{if ($4>='+str(threshold)+') ++a}END{print a}', filteredpileup])
    HAMR_accessible_bases=int(HAMR_accessible_bases);
    outfn.write(args.out_prefix+'\t'+str(HAMR_accessible_bases)+'\n')
    outfn.close()
    print "HAMR accessible bases="+str(HAMR_accessible_bases)

#output mods within normalization universe (if -n specified). Else skip this step
if (args.target_bed != 'unspecified'): 
    print "Counting modifications in each feature of target bed file..."
    counts_file=output_folder+'/'+args.out_prefix+"."+"featureCounts.bedPlus1"
    positiveCounts_file=output_folder+'/'+args.out_prefix+"."+"positiveFeatureCounts.txt"
    infn=open(counts_file,'w')
    subprocess.check_call([BEDTOOLS,'intersect','-c','-s','-b', bed_file, '-a', args.target_bed],stdout=infn)
    infn.close()
    infn=open(positiveCounts_file,'w')
    awk_parameters = 'BEGIN {OFS="\t";FS="\t"} {if ($'+str(target_bed_colnum+1)+' > 0) print $0}'
    subprocess.check_call(['awk', awk_parameters, counts_file],stdout=infn)
    infn.close()    
else:
    pass

#plot modification identity (if -tp is specified)
if (args.type_plot):
    print "plotting summary of modification types..."
    plot_file = output_folder+'/'+args.out_prefix+".mod_type.pdf"
    subprocess.check_call(['Rscript', plot_mod_type, bed_file, plot_file])

#print summary file
mods_per_acc_bases_file=output_folder+'/'+args.out_prefix+".hamr_acc_bases.txt"
outfn=open(mods_per_acc_bases_file,'w')
mods_per_acc_bases = float(true_mods)/float(HAMR_accessible_bases)*1000000
outfn.write('sample\tmods\thamr_accessible_bases\tmods_per_million_accessible_bases\n')
outfn.write(args.out_prefix+'\t'+str(true_mods)+'\t'+str(HAMR_accessible_bases)+'\t'+str(mods_per_acc_bases)+'\n')
outfn.close()

#conclusion message
print "HAMR analysis complete\n\n------------------------------\n"
print "HAMR accessible sites analyzed (read depth>=%d): %d" % (threshold, HAMR_accessible_bases) 
print "Modification sites found: " + str(true_mods) 
# final_freq_table contains sites with ref / non-ref nucleotide mixtures
print "Sites used for analysis: %s" % final_freq_table
print "Statistical testing results: %s" % raw_file
print "Modification sites + predicted types saved to: %s" % prediction_file

