#!/bin/sh

# get configuration variables
export WORKING_DIR=$(pwd)
export YAML_CONFIG=$WORKING_DIR/config_v1.yaml
export MAIL=$(egrep 'mail:' $YAML_CONFIG | awk '!/^#/ && $0=$2')
export NYUID=$(egrep 'nyuid:' $YAML_CONFIG | awk '!/^#/ && $0=$2')
export SBATCH_SCRIPTS=$(egrep 'sbatch_scripts:' $YAML_CONFIG | awk '!/^#/ && $0=$2')
export TITLE=$(egrep 'experiment_title:' $YAML_CONFIG | awk '!/^#/ && $0=$2')

# Parse input
if [[ "$#" -gt 0 ]]; then
  echo "Error: this script takes no arguments."
  echo "Usage: ./run_rnaseq.sh"
  exit -1;
fi

mkdir -p /scratch/$NYUID/reports

sbatch --output=/scratch/$NYUID/reports/slurm_rnaseq_%j.out\
       --mail-type=ALL\
       --mail-user=$MAIL\
       --export=TITLE,WORKING_DIR,MAIL,SBATCH_SCRIPTS,YAML_CONFIG\
       $SBATCH_SCRIPTS/run_rnaseq.s
