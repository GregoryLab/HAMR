# High Throughput Annotation of Modified Ribonucleotides

HAMR detects modified ribonucleotides based upon their ability to intefere with Watson-Crick base pairing. This leads to reverse transcriptase base misincorporations, which can be detected as mismatchea from the reference genome in cDNA-based RNA-seq libraries. HAMR tabulates these mismatches and tests for patterns of mismatches that cannot be explained by:
1. Sequencing errors
2. Single nucleotide polymorphisms (SNPs)
3. RNA editing 

HAMR's assumptions only hold true for a haploid or diploid genome. 

HAMR input consists of coordinate-sorted mapped reads in BAM format, and a reference genome in FASTA format. Reads should be mapped with at least 1 mismatch, and we suggest at most floor(0.06 * mean read length) mismatches. HAMR outputs predicted modification sites, and predicted modification identity.

Please cite the following papers when using HAMR:
1. [HAMR: high-throughput annotation of modified ribonucleotides](https://dx.doi.org/10.1261/rna.036806.112)
2. [Chemical Modifications Mark Alternatively Spliced and Uncapped Messenger RNAs in Arabidopsis](https://doi.org/10.1105/tpc.15.00591)

## Installing HAMR

HAMR includes pre-compiled binary programs.
If neccessary, to re-compile use

```
make clean
```

```
make
```

## Preparing Mapped Reads for HAMR
Spliced alignments (with N in CIGAR string) in BAM files will be ignored.
Alignments with insertion/deletions will be ignored. 
Only continuous, un-interrupted read alignments will be used.
Any spliced alignments should be resolved (e.g., by splitting) before
running HAMR, e.g., using GATK

```
java -jar GenomeAnalysisTK.jar -T SplitNCigarReads -R <genome.fa> -I <input.bam> -o <input.splitN.bam> -U ALLOW_N_CIGAR_READS
```

## Operating instructions

Python with version v2.7.x is preferred. Python 3+ is not yet supported. Rscript and samtools are required for running HAMR (make sure they are in searchable path).

### Usage

```
python hamr.py [align.bam] [genome.fa] [prediction_model_file] [output_dir] [output_prefix] [min_read_qual] [min_read_coverage] [seq_error_rate] [hypothesis] [max_p] [max_fdr] [min_ref_percent] [OPTIONS]
```

### Examples

#### download reference genome (e.g., hg19)

```
wget http://tesla.pcbi.upenn.edu/hamr/genomes/hg19_all_chr.fas
```

#### run HAMR in genome-wide mode

```
python hamr.py trial.human.bam  genomes/hg19_all_chr.fas models/euk_trna_mods.Rdata HAMRtest human 30 10 0.05 H4 0.01 0.05 0.05
```

#### run HAMR in target mode (BED-restricted)

```
python hamr.py trial.human.bam  genomes/hg19_all_chr.fas models/euk_trna_mods.Rdata HAMRtest human_region 30 10 0.05 H4 0.01 0.05 0.05 --target_bed region.human.bed
```

	
### REQUIRED arguments
	<align.bam>
		input BAM file
	<genome.fa>
		reference genome FASTA file (should be the same as the genome file used in mapping (input BAM). NOTE: remember to index the reference genome using samtools faidx before running HAMR.
	<min_read_qual>
		the minimum base calling quality score. All low-quality bases will removed from analysis
	<min_read_cov>
		the minimum read coverage of a genomic position for it to be analyzed
	<seq_error_rate>
		The expected percentage of mismatches based solely on sequencing error
	<hypothesis>
		The hypothesis to be tested, either "H1" or "H4" (see HAMR paper. H4 is recommended).
	<max_p>
		the maximum p-value cutoff. All sites with P-value><max_p> will be filtered out.
	<max_fdr>
		the maximum FDR cutoff
	<refpercent>
		the minimum percentage of reads that must match the reference nucleotide. All sites with reference read nucleotide proportion < <refpercent> will be filtered out
	<output_dir>
		directory for the HAMR output
        <output_prefix>
		prefix for HAMR output files
 

### OPTIONS
	--target_bed <targets.bed>
		BED file with regions of interest for HAMR analysis (this replaces the default genome-wide HAMR mode and restricts HAMR analysis to specific genomic regions listed in BED)
	--paired_ends, -pe
		indicates paired-end sequencing was used 
	--filter_ends,-fe
		excludes the first and last positions in the read from the analysis



## Copyright
	Copyright (c) 2013-2018 University of Pennsylvania

	Permission is hereby granted, free of charge, to any person obtaining a
	copy of this software and associated documentation files (the "Software"),
	to deal in the Software without restriction, including without limitation
	the rights to use, copy, modify, merge, publish, distribute, sublicense,
	and/or sell copies of the Software, and to permit persons to whom the
	Software is furnished to do so, subject to the following conditions:

	The above copyright notice and this permission notice shall be included in
	all copies or substantial portions of the Software.

	THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
	OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
	FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
	AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
	LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
	FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
	DEALINGS IN THE SOFTWARE.
