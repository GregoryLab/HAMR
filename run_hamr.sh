set -e

# HAMR parameters
minQ=30
minReadCov=50
seqErr=0.05
minP=0.01
minQ=0.05
minRefPct=0.05
BAM="trial.human.bam"
GENOME="hg19.fa"
OUTDIR="HAMRtrial"
model=models/euk_trna_mods.Rdata

# genome-wide
EXPNAME="trial"
HAMR_command="python hamr.py ${BAM} ${GENOME} ${model} ${OUTDIR} ${EXPNAME} ${minQ} ${minReadCov} ${seqErr} H4  ${minP} ${minQ} ${minRefPct}" 
echo -e "\n\nHAMR (genome-wide):\n${HAMR_command}"
${HAMR_command} || echo -e "\n\n***Failed to run:\n${HAMR_command} "

# BED-restricted
EXPNAME="trial_region"
HAMR_command="python hamr.py ${BAM} ${GENOME} ${model} ${OUTDIR} ${EXPNAME}  ${minQ} ${minReadCov} ${seqErr} H4  ${minP} ${minQ} ${minRefPct} --target_bed region.human.bed"
echo -e "\n\nHAMR (BED-restricted):\n${HAMR_command}"
${HAMR_command} || echo -e "\n\n***Failed to run:\n${HAMR_command}"

