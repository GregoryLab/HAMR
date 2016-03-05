#!/bin/bash

### tRNA annotation pipeline v 4
## Requires:
# Infernal v1.0.2
# RFAM covariance model for tRNAs: tRNA.cm
# tRNAscan-SE 1.03

# needs:
# genome annot gff (genome.gff)
# chromInfo.txt

if [ $# -lt 2 ]; then
    echo "USAGE: $0 genome genomedir" > /dev/stderr
    exit 1
fi

genome=$1
genomedir=$2
flank=47     # how much surrounding sequence to look at

trna_code=`pwd`

rfam_dir=~/data/rfam/
trnascan_dir=~/src/tRNAscan-SE-1.3/

#genome=hsa19
#genomedir=~/data/genomes/hsa19

# $genomedir/trna_cmp4
# aligner=bowtie or bwa
# -> /indexes/$aligner

### Build tRNA annotation

mkdir -p $genomedir/trna_cmp4
origdir=`pwd`
cd $genomedir/trna_cmp4

# get tRNA loci from genome annotation gff
awk '$3 == "tRNA" || $3 == "mt-tRNA"' $genomedir/${genome}.gff > trna.gff
# use one-based coords in ID
awk '$3 == "tRNA" || $3 == "mt-tRNA" {
   id=$1"_"($4)"_"$7
   print $1"\t"($4-1)"\t"$5"\t"id"\t.\t"$7 
   }' $genomedir/${genome}.gff | \
     fastaFromBed -s -fi $genomedir/seq/${genome}.fa -bed - -name -fo trna.fa

# now do it separately for nuclear and mitochondrial tRNAs
awk '$3 == "tRNA"' $genomedir/${genome}.gff > trna_nuc.gff
# use one-based coords in ID
awk '$3 == "tRNA" {
   id=$1"_"($4)"_"$7
   print $1"\t"($4-1)"\t"$5"\t"id"\t.\t"$7 
   }' $genomedir/${genome}.gff | \
     fastaFromBed -s -fi $genomedir/seq/${genome}.fa -bed - -name -fo trna_nuc.fa

awk '$3 == "mt-tRNA"' $genomedir/${genome}.gff > trna_mito.gff
# use one-based coords in ID
awk '$3 == "mt-tRNA" {
   id=$1"_"($4)"_"$7
   print $1"\t"($4-1)"\t"$5"\t"id"\t.\t"$7 
   }' $genomedir/${genome}.gff | \
     fastaFromBed -s -fi $genomedir/seq/${genome}.fa -bed - -name -fo trna_mito.fa

# build fasta and gff for tRNA-CCAs
awk '/^>/ {print $0; next} {print $0"CCA"}' trna.fa > trna_cca.fa
awk 'BEGIN{OFS="\t"} $7=="+"{$5+=3;} $7=="-"{$4-=3} {print $0}' trna.gff > trna_cca.gff

# build indexes
echo "Building Bowtie index for tRNA-CCAs ..."
mkdir -p indexes/bowtie
cd indexes/bowtie
bowtie-build ../../trna_cca.fa trna_cca
cd ../..

mkdir -p indexes/bwa
cd indexes/bwa
bwa index -p trna_cca ../../trna_cca.fa
cd ../..

##
slopBed -i trna.gff -g $genomedir/chromInfo.txt -b 47 > trna_flank47.gff

## build more detailed annotation for each tRNA locus
# use InfeRNAl to align to consensus RFAM structure
cmsearch --toponly -g -o trna.infernal $rfam_dir/tRNA.cm trna_cca.fa
# extract structure alignments from Infernal output
ruby $trna_code/annot/infernal_to_table.rb trna.infernal > trna.straln

# make linear structures
echo "Generating linear structures..."
mkdir -p linstructs
# strip ">" from fasta hdr to get locus IDs
# then take highest scoring hit for each locus
sed -e 's/^>//' trna.straln | \
    sort -rgk7 | \
    awk '$1 != prev {print;} {prev=$1}' | \
    awk '{ print $9 > "linstructs/"$1".bracket"
           print $10 > "linstructs/"$1".seq" }'
for i in linstructs/*.bracket; do
    ruby $trna_code/annot/linear_struct.rb $i > linstructs/`basename $i .bracket`.linear
done

#export PERL5LIB=$PERL5LIB:$trnascan_dir
echo "Running tRNAscan on nuclear tRNAs ..."
tRNAscan-SE -q -Q -o trnascannuc_tab.txt -f trnascannuc_str.txt trna_nuc.fa
# do a further search for organellar tRNAs 
echo "Running tRNAscan on mt-tRNAs ..."
tRNAscan-SE -O -q -Q -o trnascanorg_tab.txt -f trnascanorg_str.txt trna_mito.fa

# take maximally scoring result for each locus
tail -n+4 trnascanorg_tab.txt | cat trnascannuc_tab.txt - > trnascan_tab.txt
cat trnascannuc_str.txt trnascanorg_str.txt > trnascan_str.txt

# find loci for which tRNAscan found nothing
cut -f1 trnascan_tab.txt | tail -n+4 | sed -e 's/^ *//; s/ *$//' | \
  sort > trnascan_found.txt
grep "^>" trna.fa | sed -e 's/^>//' | sort > all_trnas
comm -1 -3 trnascan_found.txt all_trnas > trnascan_missing.txt

# generate a table with information on each locus
ruby $trna_code/annot/annot_trna.rb trnascan_tab.txt > annot.txt

# call initiator tRNAs
ruby $trna_code/annot/annotate_initiator.rb .


###########

cd $origdir
