
HAMR takes in a coordinate-sorted BAM file with mapped reads. Reads should be mapped with m>0 mismatches: patterns of mismatches will be used to detect and predict identity of modifications. HAMR outputs predicted modification sites including identity of the modifications.

NOTE: Only continuous, un-interrupted read alignments will be used.
Any spliced alignments should be resolved (e.g., by splitting) before
running HAMR (spliced alignments in BAM (if any) will be ignored). Alignments with insertion/deletions will be ignored. 

Operating instructions:

To run HAMR pipeline:
python hamr.py [align.bam] [genome.fa] [prediction_model_file] [output_dir] [output_prefix] [min_read_qual] [min_read_coverage] [seq_error_rate] [hypothesis] [max_p] [max_fdr] [min_ref_percent] [OPTIONS]

Examples:
# genome-wide
python hamr.py trial.bam  genomes/Arabidopsis_thaliana.TAIR10.25.dna.genome.fa models/euk_trna_mods.Rdata HAMRtest arabidopsis  30 10 0.01 H4 0.01 0.05 0.05

# BED-restricted
python hamr.py trial.bam  genomes/Arabidopsis_thaliana.TAIR10.25.dna.genome.fa models/euk_trna_mods.Rdata HAMRtest arabidopsis_bed 30 10 0.01 H4 0.01 0.05 0.05 --target_bed region.bed 	
	
NOTE: Python with version v2.7.x is preferred. Rscript and samtools are required for running HAMR (make sure they are in searchable path).


REQUIRED arguments:
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
		The hypothesis to be tested, either "H1" or "H4" (see HAMR
paper. H4 is recommended).
	<max_p>
		the maximum p-value cutoff. All sites with P-value><max_p>
will be filtered out.
	<max_fdr>
		the maximum FDR cutoff
	<refpercent>
		the minimum percentage of reads that must match the reference nucleotide. All sites with reference read nucleotide proportion < <refpercent> will be filtered out
	<output_dir>
		directory for the HAMR output
         
        <output_prefix>
		prefix for HAMR output files
 

OPTIONS:
	--target_bed <file.bed>
		BED file with regions of interest for HAMR analysis (this replaces the default genome-wide HAMR mode and restricts HAMR analysis to specific genomic regions listed in BED)
	--paired_ends, -pe
		indicates paired-end sequencing was used 
	--filter_ends,-fe
		excludes the first and last nucleotides of a read from the analysis



Copyright:
	Copyright (c) 2015 University of Pennsylvania

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
