# RNA-seq analysis

This repository contains all the scripts to complete RNA-seq analysis on the prince HPC at NYU. Methodology can be generalized to other systems, but they have been tailor made to be run by members of the Ercan lab. This is the second version of the .

The below instructions will outline the steps to run the ChIP-seq analysis pipeline.

### ErcanLab_ChIP-seq_analysis_v2_slurm

Author: Sevinc Ercan - se71@nyu.edu

Implementation:
Sarah Albritton - sarahea283@gmail.com
Diogo Mesquita - dam740@nyu.edu  
Lena Street - las821@nyu.edu   
Matt Paul - matthew.paul.2006@gmail.com  

Date: 05.2019

TO DO:
Replace Schematic
Write Methodology

#### RNA pipeline flowchart overview
![chip_pipeline_overview](https://github.com/ercanlab/ChIPseq/blob/master/specs/chip_pipeline_v5.png)

#### Method Overview
Sequenced reads were mapped to WS220/ce10 build of the C. elegans genome using Bowtie2 version 2.3.2 (Langmead and Salzberg, 2012). Standard mapping parameters were used for mapping with Bowtie 2. Duplicates above a expected cutoff were removed using MACS2 version 2.1.1 (Zhang, et al 2008). MACS2 determined the cutoff (typically above 1 or 2 copies) using a binomial distribution test with a p-value of 1e-5. Biological replicates were merged together with samtools (Li, et al 2009). Coverage was calculated over 100bp bins for Input subtracted ChIP enrichment scores/Chip to input ratio enrichment scores using deepTools version 3.02 (Ramírez, et al 2016). Peaks were called using MACS2. Peaks were defined using lenient settings for all biological replicates separately and using the strict settings for the merge of the biological replicates. Peaks were only considered to be true peaks if they were found in the merged replicate and in the majority of the separate replicates.

* Langmead B., Salzberg S. Fast gapped-read alignment with Bowtie 2. Nat Methods. 2012 Mar 4;9(4):357-9. doi: 10.1038/nmeth.1923.
* Zhang Y., Liu T., Meyer C.A., Eeckhoute J., Johnson D.S., Bernstein B.E., Nusbaum C., Myers R.M., Brown M., Li W., Liu X.S. Model-based analysis of ChIP-Seq (MACS). Genome Biol. 2008;9(9):R137. doi: 10.1186/gb-2008-9-9-r137.
* Li H., Handsaker B., Wysoker A., Fennell T., Ruan J., Homer N., Marth G., Abecasis G., Durbin R. and 1000 Genome Project Data Processing Subgroup. The Sequence alignment/map (SAM) format and SAMtools. Bioinformatics. 2009 Aug 15;25(16):2078-9. doi: 10.1093/bioinformatics/btp352.
* Ramírez F., Ryan D.P., Grüning B., Bhardwaj V., Kilpert F., Richter A.S., Heyne S., Dündar F., Manke T. deepTools2: A next Generation Web Server for Deep-Sequencing Data Analysis. Nucleic Acids Res. 2016 Jul 8;44(W1):W160-5. doi: 10.1093/nar/gkw257.


#### Note:
Analysis is done in two steps. First quantification of of transcript levels for each biological sample. Second differential gene expression analysis between experimental datasets. It is not necessary to do both at the same time (i.e. if you have a single RNA sample with no control sample to compare to, there is no need to do differential gene expression analysis as there is nothing to compare to).

#### 1. Create a new directory WD for the data (WD for working directory). This can be any name. Here we have used 'NewRNAseqData'.

```sh
mkdir -p /scratch/$NYUID/NewRNAseqData
cd /scratch/$NYUID/NewRNAseqData
```

#### Note:
NYUID is a placeholder. In commands given replace any instance of NYUID with your own NYUID:

```sh
mkdir /scratch/$NYUID/NewRNAseqData    #NYUID must be replaced in this directory
```

Example:
```sh
mkdir /scratch/mrp420/NewRNAseqData   #This is correct command as I have put my own NYUID (mrp420) into the directory
```

#### 2. Copy the data you want to analyze into your directory. These can be bam, fastq files or both. The directory of the new data will be given by Gencore once they have finshed sequencing. Otherwise if it is an older files the bam files can be found in the ercan lab scratch, or the locations of fastqs will be in the data google spreadsheets.

```sh
cp /scratch/cgsb/gencore/out/Ercan/DirectoryName/FileName* /scratch/$NYUID/NewRNAseqData/
```

#### 3. Rename the files to fit with the ercan lab naming system. Ensure there are NO SPACES (if you are using BAM files instead of Fastq files, in the following examples replace .fastq by .bam)

```sh
mv Oldfilename.fastq Newfilename.fastq
```

Naming Convention:

    RNA:
    <sequence_ID>_<strain>_<stage>_<replicate>.fastq
    Example: SEA38_N2_Emb_rep1.fastq

    If the experiment involved RNAi or other perturbation, indicate so in the Strain name
    Example: SEA38_N2DPY27RNAi_Emb_rep1.fastq

#### 4. Copy the analysis script to the WD

```sh
cp /scratch/cgsb/ercan/scripts/rna/slurm/run_rnaseq.sh .
```

#### 5. Copy the configuration file to the WD and edit it to add in the details of the RNAseq files

```sh
cp /scratch/cgsb/ercan/scripts/rna/slurm/config_v1.yaml .
```

After copying the configuration file (config_v1.yaml) you must edit it with the information about your data. The instructions on how to do this are in the header of the file itself. Note that the input files can be either Fastq or BAM files (simply make sure that they have the correct extension in their file name - i.e. either .fastq or .bam)

It is also possible to copy the file from google drive, edit it locally on your computer, then upload it to the HPC. The configuration file (config_v1.yaml) can be found on gdrive at:
https://drive.google.com/drive/u/1/folders/0B7OXyPwHskXaWWsyT2RZVEJ6X1U

Once edited locally on your computer you can transfer it up to the HPC with the following command:

```sh
scp path/to/config_v1.yaml $NYUID@prince.hpc.nyu.edu:/scratch/NYUID/NewRNAseq/
```

#### Note:
The config file must contain only the details for the samples you want to run. If there are extra entries leaving them blank or using a pound sign (#) to mask them does not work. They will be parsed and confuse the metadata files resulting in errors. Instead just delete blank entries.

#### 6. Run the RNA-seq pipeline. This command must be run from inside the WD.

```sh
source run_rnaseq.sh
```

This job will take several hours. The length of time will depend on how many files you are running.  
To check if it is still running can use the squeue command:

```sh
squeue -u $NYUID
```

In the output of this command look for the row named 'rnaseq' for current status and runtime so far. The command also gives the $JOBID of the job. These can be used to access the logs of the job (i.e. the printed messages). Logs will be stored in /scratch/$NYUID/reports/slurm_rnaseq_$JOBID.out  
As a sanity check you can look at the top of this log file to make sure the program read the configuration file correctly. You can look at this log file while the job is running:

```sh
less /scratch/$NYUID/reports/slurm_rnaseq_$JOBID.out
```

You will also receive emails that detail the start and end of each step in the pipeline. You can keep an eye on these outputs to check when the pipeline ends.

#### 7. Once the pipeline is finished, copy the read alignments from bowtie to archive. This command must be run from inside the WD.

```sh
cp ReadAlignments/* /archive/s/se71/ercanlab/ReadAlignments
```

#### 8. Check how the job performed

###### A. Check your emails:
You will get a series of emails for the start and end of the different jobs that the RNA-seq pipeline spawns. Check the emails that signal the end of a job and make sure that they have a COMPLETED status (as opposed to FAILED) and that the exit code is [0-0].  

###### B. Run the find errors command:
Some errors are not reported in the emails. And this command might catch some of these. If it finds no errors this command will output nothing. If it finds errors then it outputs the file where it found the errors and the corresponding error messages (note that the file name indicates in which part of the pipeline the error occurred).
```sh
source /scratch/cgsb/ercan/scripts/generic/find_errors.sh       
```
You will see something like this if there is an error. The filename indicates that there was an error while running bowtie (this example is from the ChIP-seq pipeline). Specifically in this case there was an error because the job ran out of memory:
![chip_pipeline_error_example](https://github.com/ercanlab/ChIPseq/blob/master/specs/chip_pipeline_error_example.png)

###### C. Check the log files for the job as a whole

```sh
cat /scratch/$NYUID/reports/slurm_rnaseq_$JOBID.out
```

###### D. Check the directory structure
Once the job had completed successfully the WD should look as below:
![rna_pipeline_directory_example](https://github.com/ercanlab/RNAseq/blob/master/rnaseq_output.png)

#### Note:
If there is an error and you successfully fix the problem, before you decide to rerun make sure you clear your working directory of any files associated with previous run. There should just be input files (fastq or bam), config file and RNA-seq script.

Running this code will reset directory if you are in the WD. Check there are no conflicts between rm command and names of Fastq files:

```sh
mv Fastq/* .
rm -rf bam* cou* cuff* fpk* Re* re* t*
```
If you started with BAM files (or if the issue is not at the level of mapping and you do not want to run bowtie again) you can substitute moving the contents of Fastq directory with moving contents of BAM directory.

#### Note:
These 4 examples show how ways to check that the pipeline has worked. Sometimes errors way still slip through. Keep an eye out to make sure you get all the desired output files.

###### E. Look at results
Most of the data you are interested in will be kept in the summary directory. If you are interested in specific files i.e. isoform information, you can check the relevant subdirectory for more detailed metadata files.

The summary directory will contain two types of files. (1) The summary files contain all types of mean normalized count information (i.e. mean TPMs, CPMs and FPKM) from all of the replicates for each condition. There are also the standard deviation. (2) All files contain the normalized count information for every individual replicate.

For analysis choose between using CPMs or TPMs. If you are comparing expression between genes under the same condition (i.e. within sample), use the TPM. This accounts for gene length and normalizes to sequence length variation within the pool (and therefore is better then fpkm). If you are comparing gene expression between different conditions or different strains use the CPMs or continuue on to differential expression analysis with DEseq.

#### Differential expression analysis
If you are not interested in completing differential analysis there is no need to carry on from this point. You will already have files containing CPM, TPM and FPKM from all samples. For differential expression analysis we use the Bioconductor tool DEseq (https://bioconductor.org/packages/release/bioc/html/DESeq.html).

#### 1. Get the DEseq script
Go to the folder where the initial RNAseq pipeline was run. Copy the DEseq script into this directory.

```sh
cd /scratch/NYUID/NewRNAseq
cp /scratch/cgsb/ercan/scripts/rna/slurm/run_deseq.sh .
```
#### 2. Run DEseq script
Run the DEseq script. There is no need to add any more information, as the original config file from the RNAseq pipeline will provide all the required information.
```sh
source run_deseq.sh
```

This job may take several hours. The length of time will depend on how many files you are running.  
To check if it is still running can use the squeue command:

```sh
squeue -u $NYUID
```

In the output of this command look for the row named 'deseq' for current status and runtime so far. The command also gives the $JOBID of the job. These can be used to access the logs of the job (i.e. the printed messages). Logs will be stored in /scratch/$NYUID/reports/slurm_deseq_$JOBID.out  
As a sanity check you can look at the top of this log file to make sure the program read the configuration file correctly. You can look at this log file while the job is running:

```sh
less /scratch/$NYUID/reports/slurm_deseq_$JOBID.out
```

You will also receive emails that detail the start and end of the script. You can keep an eye on these outputs to check when the pipeline ends.

#### 8. Check results and how the job performed

###### A. Check your emails:
You will get an email for the start and end of the DEseq script. Check the email that signals the end of the job and make sure that they have a COMPLETED status (as opposed to FAILED) and that the exit code is [0-0].  

###### B. Check the log files for the job as a whole

```sh
cat /scratch/$NYUID/reports/slurm_rnaseq_$JOBID.out
```
###### C. Check the directory structure
Once the job had completed successfully the WD should look as below:
![rna_pipeline_directory_example](https://github.com/ercanlab/RNAseq/blob/master/deseq_output.png).
The results can be found under the folder that is named after the experiment_title field you submitted in the RNA-seq configuration file. In this case, the results are in dpy-27-RNAi_test
