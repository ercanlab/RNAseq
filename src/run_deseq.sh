#!/bin/sh

WORKING_DIR=$(pwd)
export YAML_CONFIG=$WORKING_DIR/config_v1.yaml
if [[ ! -f "$YAML_CONFIG" ]]; then
  echo "Error: could not find the YAML config file: $YAML_CONFIG"
  exit -1;
fi

# get configuration variables
export MAIL=$(egrep 'mail:' $YAML_CONFIG | awk '!/^#/ && $0=$2')
export NYUID=$(egrep 'nyuid:' $YAML_CONFIG | awk '!/^#/ && $0=$2')
export SBATCH_SCRIPTS=$(egrep 'sbatch_scripts:' $YAML_CONFIG | awk '!/^#/ && $0=$2')

mkdir -p /scratch/$NYUID/reports

sbatch --output=$WORKING_DIR/reports/slurm_deseq_%j.out\
       --error=$WORKING_DIR/reports/slurm_deseq_%j.err\
       --mail-type=ALL\
       --mail-user=$MAIL\
       --export=YAML_CONFIG,SBATCH_SCRIPTS\
       $SBATCH_SCRIPTS/run_deseq.s
