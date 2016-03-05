
HAMR takes a bam file that has been sorted with redundant reads removed and generates a HAMR predicted_mods.txt output

Operating instructions:
	To run HAMR, run hamr_code_sk.py  following arguments while in this directory:
		python2.7 hamr_code_sk.py [bam] [genome_fas] [min_qual] [min_cov] [seq_err] [hypothesis] [max_p] [max_fdr] [refpercent] [out_file] (optional)--paired_ends (optional)--filter_ends

	Required Arguments:
		bam 		-The Directory location and name of bam file consisting of nonredundant reads that you wish to use.
		genome_fas 	-The Directory Genome fasta file that you wish to use. WARNING: remember to index the reference using samtools faifx before running hamr_code_sk.py
		min_qual	-The minimum threshhold of a read's quality score for it to be analyzed
		min_cov		-The minimum coverage of a nucleotide for it to be analyzed
		seq_err		-The percentage of mismatches based solely on sequencing error
		hypothesis	-The hypothesis to be tested, either "H1" or "H4"
		max_p		-The maximum p-value cutoff
		max_fdr		-The maximum FDR cutoff
		refpercent	-The percentage of reads that must match the reference nucleotide
		out_file	-The directory location and filename for the  HAMR predicted_mods output
 
	Optional Arguements:
		--paired_ends, -pe	Indicates paired-end sequencing was used
		--filter_ends,-fe	Excludes the first and last nucleotides of a read from the analysis



Other Info:
	The compiled rnapileup script include in the directory is statically compiled. It can be run without the samtools library installed. In the event of bugfixes in the samtools library, however, rnapileupV2.2.1 will need to be recompiled before the bugfixes take effect.


Copyright:
	Copyright (c) 2013 University of Pennsylvania

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
