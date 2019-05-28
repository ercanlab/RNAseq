################################################################################
## Performs standard DEseq analysis. Also saves DEseq results                 ##
##                                                                            ##
## usage: Rscript Initial_DESeq2_Analysis.R YAML_CONFIG                       ##
##                                                                            ##
## positional arguments:                                                      ##
##  YAML_CONFIG                                                               ##
##  An example of the configuration file is on Prince at:                     ##
##  /scratch/cgsb/ercan/scripts/rna/slurm/config_deseq.yaml                   ##
##  or on gdrive:                                                             ##
##  https://drive.google.com/open?id=1HCEOuFQQsObFf5QVLvF3n0a-894Ts9Ze        ##
##                                                                            ##
################################################################################


################################################################################
## The package DESeq2 provides methods to test for differential expression by ##
## use of negative binomial generalized linear models                         ##
## Documentation for DESSeq2 is found at:                                     ##
## https://bioconductor.org/packages/3.7/bioc/vignettes/DESeq2/inst/doc/DESeq2.html ##
################################################################################


## Load up required packages
library('stringr')
suppressMessages(library('DESeq2'))
library('yaml')
suppressMessages(library('gplots'))
library('RColorBrewer')
library("openxlsx")
suppressMessages(library("dplyr"))

####################################
##### Define utility functions #####
####################################

## Extract condition and rep names to add to normalized data frame.
getColNames <- function(sampleCondition){
  curr_condition <- ""
  col_names <- vector(mode="character", length=length(sampleCondition))
  i=1
  for (cond in sampleCondition){
    if (cond != curr_condition){
      curr_condition=cond
      rep_n <- 1
    }
    col_names[i] <- paste0(curr_condition, "_r", rep_n)
    rep_n <- rep_n + 1
    i <- i + 1
  }
  col_names
}

## Make a dataframe of rpkms
make_fpkm_df <- function(dir, files, sampleCondition, to_threshold=FALSE){
  # todo
  g1 <- 'T21D12.9'
  g2 <- 'F07C6.4'
  keep <- c('gene_id', 'FPKM')
  df = NULL
  for (i in seq_along(files)){
    path <- file.path(dir, files[i])
    fpkm <- read.csv(path, sep='\t')[keep]
    fpkm <- filter(fpkm,gene_id!=g1)
    fpkm <- filter(fpkm,gene_id!=g2)
    if (to_threshold){
      fpkm <- fpkm[fpkm['FPKM'] > 1,]
    }
    if (is.null(df)){
      df <- fpkm
    } else {
      df <- merge(df, fpkm, by='gene_id', suffixes=c(files[i-1], files[i]))
    }
  }
  rownames(df) <- df$gene_id
  df <- df[, !(names(df) %in% 'gene_id')]
  names(df) <- getColNames(sampleCondition)
  df$gene_id<-rownames(df)
  return(df)
}

#This function does comparisons of gene expression between conditions using DEeq
do.DEseq<-function(tableA,tableB,conditionA,conditionB,condA.avg, condB.avg, condition_type){
  TABLE<-rbind(tableA, tableB)
  DDS<-DESeqDataSetFromHTSeqCount(sampleTable = TABLE, directory = counts_dir, design= ~ condition)
  filt.DDS <- DDS[ rowSums(counts(dds)) > 1, ]
  if (condition_type[conditionA] == 'control'){
    filt.DDS$condition <- relevel(filt.DDS$condition, ref = conditionA)
  } else if (condition_type[conditionB] == 'control'){
    filt.DDS$condition <- relevel(filt.DDS$condition, ref = conditionB)
  }
  DDS<-DESeq(filt.DDS)
  DDS.res<-results(DDS, alpha=0.05)
  return(DDS.res)
}

#Modify data frame to give appropriate headers to filenames
updatePairwiseDataFrame <- function(df, res, col.basename){
  df[paste0(col.basename,'.log2.deseq')] <- res$log2FoldChange
  df[paste0(col.basename,'.log2')] <- res$unadj.log2
  df[paste0(col.basename,'.log2.pval')] <- res$padj
  return(df)
}

## Plot correlations between chr heatmaps
plotSpearmenHeatmap <- function(df, file, title){
  c <- cor(df, method="spearman", use="complete.obs")
  pdf(file=file)
  mypalette<-brewer.pal(11,"Spectral")
  morecols<-colorRampPalette(mypalette)
  heatmap.2(c, main=title,col=rev(morecols(50)),trace="none", margins=c(17,6),denscol='black')
  dev.off()
}

# Plot pairwise counts comparisons along with Rsquared value for every pairwise comparison
plotPairwiseCountData <- function(df, file){
  df<-as.data.frame(df)
  reps.list<-colnames(df)
  n_reps = length(reps.list)
  pdf(file=filepath)
  par(mfrow=(c(3,3)))
  par(mar=c(4,4,1,1))
  for(i in 1:n_reps){
    x.data<-df[,reps.list[i]]
    for(j in 1:n_reps){
      if(j!=i){
        y.data<-df[,reps.list[j]]
        plot_limits<-c(0,signif(range(df)[2]+(range(df)[2]/10), digits = 2))
        plot((x.data),(y.data),xlab=reps.list[i],ylab=reps.list[j],xlim=plot_limits,ylim=plot_limits) #you may want to change the x and y limits
        r2<-round(summary(lm(y.data ~ x.data))$r.squared,4)
        text(plot_limits[2]*(1.5/5),plot_limits[2]*(4/5),paste('r2=',r2))
      }else{
      plot_limits<-c(0,signif(range(df)[2]+(range(df)[2]/10), digits = 2))
      plot(0,0,type='l',xlab='',ylab='',axes=F,xlim=plot_limits,ylim=plot_limits)
      text(x=plot_limits[2]/(1.9), y=plot_limits[2]/2,reps.list[i], cex=2)
      }
    }
  }
  dev.off()
}

y.data<-df[,reps.list[1]]
round(summary(lm(y.data ~ x.data))$r.squared,4)

#Helps prepare dataframe for plotting and removes mitochondrial dna
joinChromosome <- function(df, c_elegans_annots){
  df <- as.data.frame(df)
  df$gene.name <- rownames(df)
  merged <- merge(df, c_elegans_annots, by.x="gene.name", by.y="Sequence.Name.(Gene)")
  merged <- filter(merged, Chr.Name != "MtDNA")
  return(merged)
}


#Make a boxplot of log2FC for values of each chromosome
boxPlotDESeqByChromosome <- function(df, file, c_elegans_annots, title){
  merged <- joinChromosome(df, c_elegans_annots)
  pdf(file=file)
  ylow = min(merged$log2FoldChange, na.rm=TRUE)
  yhigh = max(merged$log2FoldChange, na.rm=TRUE)
  boxplot(log2FoldChange~Chr.Name,data=merged, main=title,
          xlab="Chromosome", ylab="Log2 Fold Change", ylim=c(ylow, yhigh))

  stats_by_chromosome <-
    merged %>%
    group_by(Chr.Name) %>%
    summarize(mean = mean(log2FoldChange, na.rm = TRUE), std=sd(log2FoldChange, na.rm = TRUE))

  mean <- sapply(stats_by_chromosome$mean, sprintf, fmt="%.3f")
  std <- sapply(stats_by_chromosome$std, sprintf, fmt="%.3f")
  text(c(1), y=0.9*ylow, labels=c("mean:"))
  text(seq(6), y=0.98*ylow, labels=mean)
  dev.off()
}

#Scatterplot comparing gene expression. X and A are highlighted
scatterPlotDeseq <- function(res, c_elegans_annots, file, title){
  merged <- joinChromosome(res, c_elegans_annots)
  merged$is_X <- merged$Chr.Name == "X"
  pdf(file=file)
  ylow <- min(merged$log2FoldChange, na.rm=TRUE)
  yhigh <- max(merged$log2FoldChange, na.rm=TRUE)
  plotMA(merged[c("baseMean", "log2FoldChange", "is_X")], main=title,
         xlab="Base mean", ylab="Log2 Fold Change", colSig='cyan',
         ylim=c(ylow, yhigh))

  legend('bottomright','groups',c("X genes","Autosomal genes"), pch = 16,
         col=c('cyan', 'black'),ncol=2,bty ="n")
  dev.off()
}

## This function is used to validate that the config file has right attributes
validateConfig <- function(conf){
  required <- c(
    'experiment_title',
    'infiles',
    'c_elegans_wbid_to_gene',
    'c_elegans_annots',
    'nyuid',
    'mail',
    'sbatch_scripts'
  )
  missing <- required[!(required %in% names(conf))]
  if (length(missing) > 0){
    stop("Attributes missing from configuration: ", paste(missing, collapse="; "))
  }
  invisible(lapply(conf$infiles, validateInfiles))
}

## Ensure that config file has inputs within it
validateInfiles <- function(x){
  required <- c(
    "id",
    "fastq",
    "condition",
    "type"
  )
    missing <- required[!(required %in% names(x))]
    if (length(missing) > 0){
      stop("Attributes missing from element in infiles: ", paste(missing, collapse="; "))
    }
}

#Uses gene names from gtf file to extract out wormbase ID
getGenesWbId <- function(c_elegans_annots, genes){
  gene.to.wbid<-read.table(file=c_elegans_annots,header=F,stringsAsFactors=F)
  colnames(gene.to.wbid)<-c('gene','wbid')
  relevant_genes_ix <- match(genes, gene.to.wbid$gene)
  wbid<-gene.to.wbid$wbid[relevant_genes_ix]
}

## Read in annotation file and transform to remove duplicates subset to relevant
## information
getCelegansAnnotations <- function(file){
  c_elegans_annots <- read.xlsx(file)
  relevant_cols <- c("Gene.WB.ID", "Sequence.Name.(Gene)", "Chr.Name")
  c_elegans_annots <- c_elegans_annots[,relevant_cols]
  c_elegans_annots <- c_elegans_annots[!duplicated(c_elegans_annots[,"Sequence.Name.(Gene)"]),]
  c_elegans_annots <- c_elegans_annots[complete.cases(c_elegans_annots[,relevant_cols]),]
  return(c_elegans_annots)
}

pValueLogFoldChangeByChromosome <- function(df){
  null_mean <- NULL
  for (chr in unique(df$Chr.Name)){
    if (chr %in% c("MtDNA", "X")){
      next
    }
    data <- df[df$Chr.Name == chr, "log2FoldChange"]
    m <- mean(data, na.rm=TRUE)
    if (is.null(null_mean)){
      null_mean <- c(m)
    } else {
      null_mean <- c(null_mean, m)
    }
  }
  null_mean <- mean(null_mean)

  res <- list()
  for (chr in unique(df$Chr.Name)){
    if (chr == "MtDNA"){
      next
    }
    data <- df[df$Chr.Name == chr, "log2FoldChange"]
    n <- length(data)
    t <- (mean(data, na.rm=TRUE) - null_mean) / (sd(data, na.rm=TRUE)/sqrt(n))
    p <- 2*pt(-abs(t),df=n-1)
    res[chr] <- p
  }
  return(res)
}

## Creates a list of input files and experimental metadata
readInFiles <- function(infiles){
  idx = 1
  id_to_idx = new.env()
  for (ele in infiles){
    id  <- toString(ele$id)
    if (is.null(id_to_idx[[id]])){
      id_to_idx[[id]] <- idx
      idx = idx + 1
    }
  }
  num_files <- length(id_to_idx)
  bam_suffix <- "_accepted_hits_counts.txt"
  fpkm_suffix <- "_cufflinks.genes.fpkm_tracking"
  files_by_id <- vector("list", num_files)
  conditions <- vector("list", num_files)
  types <- vector("list", num_files)

  for (ele in infiles){
    idx  <- id_to_idx[[toString(ele$id)]]
    if (is.null(files_by_id[[idx]])){
      files_by_id[[idx]] <- list(ele$fastq)
    } else {
      files_by_id[[idx]] <- list(files_by_id[[idx]], ele$fastq)
    }
    conditions[[idx]] <- ele$condition
    types[[idx]] <- ele$type
  }
  count_files <- unlist(lapply(files_by_id, function(x)
    paste0(paste0(str_replace_all(unlist(x), c(".fastq.gz"="",".fastq"="")), collapse="_"), bam_suffix)))
  fpkm_files <- unlist(lapply(files_by_id, function(x)
    paste0(paste0(str_replace_all(unlist(x), c(".fastq.gz"="",".fastq"="")), collapse="_"), fpkm_suffix)))
  return(list(
    "count_files"=count_files,
    "fpkm_files"=fpkm_files,
    "conditions"=unlist(conditions),
    "types"=unlist(types)
  ))
}

########################
##### Start script #####
########################

## Load in the arguments from command line (location of YAML file) ##
args <- commandArgs(trailingOnly = TRUE)
if (length(args)==0) {
  stop("At least one argument must be supplied (the deseq yaml config file).n", call.=FALSE)
}

## HTSeq input - read in config file and check it has the right format
conf <- yaml.load_file(args[1])
validateConfig(conf)

## Define file paths to get to processed experimental files ##
ercan_rna <- "/scratch/cgsb/ercan/rna/deseq"
scratch_dir <- file.path(ercan_rna, conf$experiment_title)
# Checks to see if the experimental comparisons have been made before
if (dir.exists(scratch_dir)){
  stop(sprintf("the experiment title '%s' already exists. Change this title or remove the %s directory.", conf$experiment_title, scratch_dir))
}

## Define directories for input/output files
working_dir <- dirname(args[1])
fpkm_dir <- file.path(working_dir, "fpkm")
counts_dir <- file.path(working_dir, "counts")
deseq_dir <- file.path(working_dir, conf$experiment_title)
out_dir <- file.path(deseq_dir, "results")
dir.create(out_dir,showWarnings = FALSE, recursive = TRUE)

## Creates a list of input files and experimental metadata
inFiles_data <- readInFiles(conf$infiles)

conditions <- levels(factor(sampleCondition))
sampleType <- inFiles_data$type
condition_type <- vector()
#for loop creates a vector of treatments
for (i in 1:length(sampleCondition)){
  cond <- sampleCondition[i]
  if (!(cond %in% names(condition_type))) {condition_type[cond] = sampleType[i]}
}
sampleFiles <- inFiles_data$count_files
fpkm_files <- inFiles_data$fpkm_files

## Read in the HTseq outputs into dds format ##

sampleTable <- data.frame(sampleName=sampleFiles, fileName=sampleFiles, condition=sampleCondition)
dds <- DESeqDataSetFromHTSeqCount(sampleTable = sampleTable, directory = counts_dir, design= ~ condition)

## Filter out low/no count genes. Need more then 1 count per gene ##
filt.dds <- dds[ rowSums(counts(dds)) > 1, ]

print(paste0('Of the ', nrow(dds), ' genes, ', nrow(filt.dds), ' have >1 reads summed across conditions. ', (nrow(filt.dds)/nrow(dds))*100, '% of genes remain'))

## Get normalized count values ##
#DEseq works to estimate variability of samples and applies -ve binomial model
filt.dds<-DESeq(filt.dds)

#Extract out the normalized count value from the DEseq analysis object
normalized.count.data<-(assays(filt.dds)[["mu"]])
colnames(normalized.count.data)<-getColNames(sampleCondition)
normalized.count.data <- as.data.frame(normalized.count.data)

#Extract out gene names, WB gene names and chr for annotated genes
c_elegans_annots <- getCelegansAnnotations(conf$c_elegans_annots)


## Make heatmap plot of FPKM values of replicates based on spearman correlation
#DATANAME <- paste0(conditions, collapse = 'vs')
#thresholded_fpkm_data <- make_fpkm_df(fpkm_dir, fpkm_files, sampleCondition, to_threshold = TRUE)
## Add in the Chr info and Wormbase ID
#merged <- merge(thresholded_fpkm_data, c_elegans_annots, by.x="gene_id", by.y="Sequence.Name.(Gene)")
## Remove mitochondrial data
#merged <- filter(merged, Chr.Name != "MtDNA")
#not_for_plot <- c('Chr.Name', 'gene_id', 'Gene.WB.ID')
#for (chr in unique(merged$Chr.Name)){
#  #Set names and titles
#  filepath <- file.path(out_dir, paste0(DATANAME,'.heatmap.spearman.thresholded.fpkm.', chr, '.pdf'))
#  title <- paste0("Correlations of FPKM of Chr ", chr, " genes \n DATANAME")
#  #pull out a single chromosome
#  df<-filter(merged, Chr.Name == chr)
#  #Drop uneeded annotation from plotting
#  df <- df[, !(names(df) %in% not_for_plot)]
#  #plot the correlations
#  plotSpearmenHeatmap(df, filepath, title)
#}
# plot all chromosomes together
#filepath <- file.path(out_dir, paste0(DATANAME,'.heatmap.spearman.thresholded.fpkm.all.pdf'))
#title <- paste0("Correlations of FPKM values \n DATANAME")
#plotSpearmenHeatmap(merged[, !(names(merged) %in% not_for_plot)], filepath, title)

## get average count data table ##
#Extract one condition at a time from normalized data, and then calculate mean
# and stdev for each gene under that condition
conditions_avg <- list()
for (cond in conditions){
  conditions_avg[[paste0(cond, '_mean')]] <- rowMeans(normalized.count.data[,sampleCondition == cond, drop=FALSE])
  #conditions_avg[[paste0(cond, '_stdev')]] <- rowSds(data.matrix(normalized.count.data[,sampleCondition == cond, drop=FALSE]))
}

#Add gene names to average count data
genes <- names(conditions_avg[[1]])
wbid <- getGenesWbId(conf$c_elegans_wbid_to_gene, genes)
avg.count.data<-data.frame(wbid=wbid, conditions_avg)

#Save the normalized counts
filepath <- file.path(out_dir,'avg.count.data.txt')
write.table(format(avg.count.data, digits=2, scientific=FALSE),file=filepath,row.names=T,col.names=T,quote=F,sep='\t')


#######################################
## plot replicates pairwise with Rsquared values ##
#######################################
#filepath <- file.path(out_dir, paste0(DATANAME, '.replicates.counts.vs.counts.Rsq.pdf'))
#plotPairwiseCountData(normalized.count.data, filepath)

########################
## do pairwise comparisons ##
########################

#Generate a list with of all count data and its condition
condTables = list()
for (cond in conditions){
  condTables[[cond]] = sampleTable[sampleTable['condition']==cond,]
}

n_conditions <- length(conditions)
pairwise_res_df <- data.frame(row.names = genes)
for (i in seq(n_conditions-1)){
  for (j in seq(i+1,n_conditions)){
  if(i!=j){
    deseq.df <- do.DEseq(condTables[[i]], condTables[[j]],
                        conditions[[i]], conditions[[j]],
                        conditions_avg[[i]], conditions_avg[[j]], condition_type)

    basename <- paste0(conditions[[i]],'vs',conditions[[j]])
    filepath <- file.path(out_dir, paste0(basename,'.deseq.txt'))
    write.table(format(as.data.frame(deseq.df), digits=2, scientific=FALSE),file=filepath,row.names=T,col.names=T,quote=F,sep='\t')
    #Modify data frame col names and minimise to just FC values and pvalues
    pairwise_res_df <- updatePairwiseDataFrame(pairwise_res_df, deseq.df, basename)

    filepath <- file.path(out_dir, paste0(basename,'.deseq.boxplot.by.chromosome.pdf'))
    boxPlotDESeqByChromosome(deseq.df, filepath, c_elegans_annots, "Log Fold Change By Chromosome")

    filepath <- file.path(out_dir, paste0(basename,'.deseq.scatterplot.pdf'))
    scatterPlotDeseq(deseq.df, c_elegans_annots, filepath, "Log Fold change vs expression")
  }
}}

#### Make a table with all of the relevant data ####
genes <- names(conditions_avg[[1]])
wbid <- getGenesWbId(conf$c_elegans_wbid_to_gene, genes)
samplemeans.df<-as.data.frame(conditions_avg)

summary.data.df<-as.data.frame(cbind(wbid,samplemeans.df, pairwise_res_df))
filepath <- file.path(out_dir, 'deseq.summaryoverview.txt')
write.table(format(summary.data.df, digits=2, scientific=FALSE),file=filepath,row.names=T,col.names=T,quote=F,sep='\t')

# todo: organize folder - should have conf in top dir then a counts folder, a fpkm folder and a output folder
file.copy(args[1], deseq_dir)
file.copy(deseq_dir, ercan_rna, recursive = TRUE)
