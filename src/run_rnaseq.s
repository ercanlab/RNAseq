#!/bin/bash
#
#BATCH --verbose
#SBATCH --job-name=rnaseq
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --mem=50GB
#SBATCH --mail-type=ALL

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

#Start of workflow

#Load in packages
module load samtools/intel/1.6

#Set directory where scripts will be found
ercan_rna="/scratch/cgsb/ercan/rna"

#Create directory to store generate reports
mkdir -p $WORKING_DIR/reports
cd $WORKING_DIR

echo "Running RNAseq for experiment '$TITLE'"

### MAP THE FASTQ FILES WITH TOPHAT

echo "Unzipping files"
gunzip *.gz

#This will assay how many fastq files we have to map
ls *fastq > files.txt
n=$(wc -l files.txt | awk '$0=$1')

#THis will run tophat to map the fastq files to the genome. Anything mapped
#succesfully will be found in accepted hits file in BAM directory. Other output
#files cn be found in the tophat directory i.e. unmapped reads
echo "Running TopHat..."
job_out=$(sbatch --output=$WORKING_DIR/reports/slurm_tophat_%j.out\
                 --error=$WORKING_DIR/reports/slurm_tophat_%j.err\
                 --mail-type=ALL\
                 --mail-user=$MAIL\
                 --array=1-$n\
                 $SBATCH_SCRIPTS/doTopHat.s)

wait_for_job "$job_out"
echo "TopHat finished"

#Clean up list of input files
rm $WORKING_DIR/files.txt

#Clean up working directory.
#Move Fastq files to their own directory
mkdir Fastq
mv *fastq Fastq
#Move tophat files to their own directory
# NOTE: if genome changes then this will also need to be changed.
mkdir tophat
mv *_ce10_tophat/ tophat

# Move the read alignments to archive.
mkdir -p ReadAlignments
mv *txt ReadAlignments

# (for ease) Copy all bamfiles to a new folder.
# Make a new folder.
mkdir BAM
#Copy all bamfiles into the new folder and cd into it.
cp tophat/*_ce10_tophat/*accepted_hits.bam BAM
cd BAM

#Rename files to remove accepted_hits suffixes
SUFF='_accepted_hits.bam'
suff='.bam'
for i in $(ls *$SUFF)
do
  mv -f $i ${i%$SUFF}$suff
done

# Copy the bam files to the ercanlab folder.
cp *.bam $ercan_rna/bam
# Copy the tophat folders to the ercanlab folder.
cp -r tophat/*_ce10_tophat $ercan_rna/tophat/

### MERGE ANY TECHNICAL REPLICATES

# Need to merge technical replicates from resequencing. Will parses YAML_CONFIG
# to find the id# for every library. If there are repeats because of
# resequencing samtools will be used to merge bamfiles together.

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
      array[$i]=${line%%.fastq*}.bam
      if [ $i -gt 0 ]; then
        out=${out}_
      fi
      out=${out}${line%%.fastq*}
      i=$(($i+1))
    done < <(cat $f)
    if [ $i -ge 2 ]; then
      to_remove+=(${array[*]})
      out=${out}.bam
      samtools merge $out ${array[*]}
    fi
    unset array
  done


  rm ${to_remove[*]}
  rm -rf replicates
}

#Apply merge bams function
merge_bams

###  GENERATE FPKMs BY ASSIGNING READS TO GENES
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

#This will add the library infromation to the fpkm output files from cufflinks
for subdir_path in */; do
  # remove ending forward slash
  subdir=${subdir_path%\/}
  mv $subdir/genes.fpkm_tracking $subdir/$subdir.genes.fpkm_tracking;
done

#Move cufflinks outputs to their own directory
cd $WORKING_DIR
mkdir -p cufflinks
mv bamfiles/*/ cufflinks/

#Save the outputs to ERcan lab directories
cp -r cufflinks/* $ercan_rna/cufflinks/

### MAKE THE COUNT FILES

# count how many bam files there are to count
cd BAM
ls *bam > files.txt
n=$(wc -l files.txt | awk '$0=$1')

#Run HT-seq cocunt onbam files to get counts per gene
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

#Move count data to its own directory
mkdir -p counts
mv bamfiles/*counts* counts/
cp counts/*counts* $ercan_rna/counts/
