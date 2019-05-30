#!/bin/sh
#
#BATCH --verbose
#SBATCH --job-name=htseq
#SBATCH --output=/scratch/las821/reports/slurm_htseq_%j.out
#SBATCH --error=/scratch/las821/reports/slurm_htseq_%j.err
#SBATCH --time=1:00:00
#SBATCH --nodes=1
#SBATCH --mem=80GB
#SBATCH --mail-type=ALL
#SBATCH --mail-user=las821@nyu.edu

module load htseq/intel/0.6.1p1
module load samtools/intel/1.6

cd $WORKING_DIR/BAM
val=$SLURM_ARRAY_TASK_ID
file=`sed -n ${val}p files.txt`

perl /scratch/cgsb/ercan/scripts/rna/getRNASeqCounts_stranded.pl $file
exit 0;
