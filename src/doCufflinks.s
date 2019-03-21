#!/bin/sh
#
#BATCH --verbose
#SBATCH --job-name=cufflinks
#SBATCH --output=/scratch/las821/reports/slurm_cufflinks_%j.out
#SBATCH --error=/scratch/las821/reports/slurm_cufflinks_%j.err
#SBATCH --time=3:00:00
#SBATCH --nodes=1
#SBATCH --mem=80GB
#SBATCH --mail-type=ALL
#SBATCH --mail-user=las821@nyu.edu

module load tophat/intel/2.1.1
module load bowtie2/intel/2.2.9
module load cufflinks/2.2.1

#module load cufflinks/intel/2.2.0
#module load tophat/intel/2.0.12
#module load bowtie/intel/1.1.1

### REQUIRED MODULES

if [ "$SLURM_JOBTMP" == "" ]; then

    export SLURM_JOBTMP=/state/partition1/$USER

    mkdir -p $SLURM_JOBTMP

fi

export TMPDIR=$SLURM_JOBTMP

#clear tmp

rm -rf $TMPDIR



val=$SLURM_ARRAY_TASK_ID
file=`sed -n ${val}p files.txt`

perl /scratch/cgsb/ercan/scripts/rna/doCufflinks_stranded.pl $file
exit 0;
