#!/bin/sh
#
#BATCH --verbose
#SBATCH --job-name=deseq
#SBATCH --nodes=1
#SBATCH --tasks-per-node=1
#SBATCH --mem=2GB
#SBATCH --time=01:00:00

module purge
module load r/intel/3.4.2

Rscript $SBATCH_SCRIPTS/Initial_DESeq2_Analysis_v2.R $YAML_CONFIG