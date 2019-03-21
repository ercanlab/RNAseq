#!/bin/bash
#
#BATCH --verbose
#SBATCH --job-name=rnaseq
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --mem=50GB
#SBATCH --mail-type=ALL

merge_bams(){
  rm -rf replicates
  mkdir -p replicates

  curr_id=-1
  egrep -A 3 "id: [0-9]+" $YAML_CONFIG | while read -r line ; do
    if [ -n "$(grep 'id:' <<< $line)" ]; then
      curr_id="$(cut -d ' ' -f 3 <<< $line)"
    
    elif [ -n "$(grep 'fastq:' <<< $line)" ]; then
      fastq="$(cut -d ' ' -f 2 <<< $line)"
      echo $fastq >> "replicates/reps_${curr_id}.txt"
    fi
  done

  to_remove=()
  for f in replicates/*; do
    i=0
    out=""
    while read -r line; do
      array[$i]=${line%%.fastq*}_accepted_hits.bam
      if [ $i -gt 0 ]; then
        out=${out}_
      fi 
      out=${out}${line%%.fastq*}
      i=$(($i+1))
    done < <(cat $f)
    if [ $i -ge 2 ]; then
      to_remove+=(${array[*]})
      out=${out}_accepted_hits.bam
      samtools merge $out ${array[*]}
    fi
    unset array
  done


  rm ${to_remove[*]}
  rm -rf replicates
}

##
# Waits for a sbatch job to finish.
# Input: string, the message returned by the sbatch command
#        e.g: 'Submitted batch job 4424072'
##
wait_for_job(){
  # extract only the jobid from the job output
  jobout="$1"
  jobid="${jobout##* }"

  is_running=$(squeue -j $jobid | wc -l | awk '$0=$1')
  while [ $is_running -gt 1 ]
  do
    sleep 15
    is_running=$(squeue -j $jobid | wc -l | awk '$0=$1')
  done
}

module load samtools/intel/1.6

ercan_rna="/scratch/cgsb/ercan/rna"

mkdir -p $WORKING_DIR/reports
cd $WORKING_DIR

echo "Running RNAseq for experiment '$TITLE'"

echo "Unzipping files"
gunzip *.gz

ls *fastq > files.txt
n=$(wc -l files.txt | awk '$0=$1')

echo "Running TopHat..."
job_out=$(sbatch --output=$WORKING_DIR/reports/slurm_tophat_%j.out\
                 --error=$WORKING_DIR/reports/slurm_tophat_%j.err\
                 --mail-type=ALL\
                 --mail-user=$MAIL\
                 --array=1-$n\
                 $SBATCH_SCRIPTS/doTopHat.s)

wait_for_job "$job_out"
echo "TopHat finished"

rm $WORKING_DIR/files.txt

mkdir FastQ
mv *fastq FastQ
mkdir tophat
mv *_ce10_tophat/ tophat

# Move the read alignments to archive.
mkdir -p ReadAlignments
mv *txt ReadAlignments

# Copy the bam files to the ercanlab folder. 
cp tophat/*_ce10_tophat/*accepted_hits.bam $ercan_rna/bam
# Copy the tophat folders to the ercanlab folder. 
cp -r tophat/*_ce10_tophat $ercan_rna/tophat/

# (for ease) Copy all bamfiles to a new folder. 
# Make a new folder.
mkdir bamfiles
#Copy all bamfiles into the new folder and cd into it.
cp tophat/*_ce10_tophat/*accepted_hits.bam bamfiles

cd bamfiles

merge_bams

ls *bam > files.txt
n=$(wc -l files.txt | awk '$0=$1')

echo "Running Cufflinks..."
job_out=$(sbatch --output=$WORKING_DIR/reports/slurm_cufflinks_%j.out\
                 --error=$WORKING_DIR/reports/slurm_cufflinks_%j.err\
                 --mail-type=ALL\
                 --mail-user=$MAIL\
                 --array=1-${n}\
                 $SBATCH_SCRIPTS/doCufflinks.s)

wait_for_job "$job_out"
echo "Cufflinks finished"

for subdir_path in */; do
  # remove ending forward slash
  subdir=${subdir_path%\/}
  mv $subdir/genes.fpkm_tracking $subdir/$subdir.genes.fpkm_tracking;
done

cd $WORKING_DIR
mkdir -p cufflinks
mv bamfiles/*/ cufflinks/

cp -r cufflinks/* $ercan_rna/cufflinks/

mkdir -p fpkm
cp cufflinks/*/*.genes.fpkm_tracking fpkm

cd bamfiles
ls *bam > files.txt
n=$(wc -l files.txt | awk '$0=$1')

echo "Running HtSeq-Count..."
job_out=$(sbatch --output=$WORKING_DIR/reports/slurm_htseq_%j.out\
                 --error=$WORKING_DIR/reports/slurm_htseq_%j.err\
                 --mail-type=ALL\
                 --mail-user=$MAIL\
                 --array=1-$n\
                 --export=WORKING_DIR\
                 $SBATCH_SCRIPTS/getCounts2.s)

wait_for_job "$job_out"
echo "HtSeq-Count finished"
cd $WORKING_DIR
mkdir -p counts
mv bamfiles/*counts* counts/
cp counts/*counts* $ercan_rna/counts/
