# genome-wide
python hamr.py trial.bam  genomes/Arabidopsis_thaliana.TAIR10.25.dna.genome.fa models/euk_trna_mods.Rdata HAMRtest arabidopsis  30 10 0.01 H4  1 0.05 0.05

# BED-restricted
python hamr.py trial.bam  genomes/Arabidopsis_thaliana.TAIR10.25.dna.genome.fa models/euk_trna_mods.Rdata HAMRtest arabidopsis_bed 30 10 0.01 H4  1 0.05 0.05 --normalization_bed region.bed 

