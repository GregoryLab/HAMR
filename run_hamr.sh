# genome-wide
python hamr_code_v2.2.3.pavel.py trial.bam  /mnt/data/users/pkuksa/datasets/arabidopsis/Arabidopsis_thaliana.TAIR10.25.dna.genome.fa models/euk_trna_mods.Rdata testoutnosplice testnosplice 30 10 0.01 H4  1 0.05 0.05

# BED-restricted
python hamr_code_v2.2.3.pavel.py trial.bam  /mnt/data/users/pkuksa/datasets/arabidopsis/Arabidopsis_thaliana.TAIR10.25.dna.genome.fa models/euk_trna_mods.Rdata testoutnosplice2 testnosplice 30 10 0.01 H4  1 0.05 0.05 --normalization_bed region.bed 

