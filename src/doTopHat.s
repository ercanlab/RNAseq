#!/bin/sh
#
#BATCH --verbose
#SBATCH --job-name=tophat
#SBATCH --output=/scratch/las821/reports/slurm_tophat_%j.out
#SBATCH --error=/scratch/las821/reports/slurm_tophat_%j.err
#SBATCH --time=5:00:00
#SBATCH --nodes=2
#SBATCH --mem=50GB
#SBATCH --mail-type=ALL
#SBATCH --mail-user=las821@nyu.edu

module load samtools/intel/1.6
module load tophat/intel/2.1.1
module load bowtie2/intel/2.2.9

val=$SLURM_ARRAY_TASK_ID
file=`sed -n ${val}p files.txt`

perl /scratch/cgsb/ercan/scripts/rna/getAlignment_Tophat.pl $file
exit 0;
