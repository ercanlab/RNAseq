##########################################################################################
## Performs standard DEseq analysis. Also saves DEseq results
##
## usage: Rscript Initial_DESeq2_Analysis.R YAML_CONFIG
##
## positional arguments:
##  YAML_CONFIG           the configuration file. See /scratch/cgsb/ercan/scripts/rna/slurm/config_deseq.yaml
##                        in Prince or https://drive.google.com/open?id=1HCEOuFQQsObFf5QVLvF3n0a-894Ts9Ze for an example
##                        
#########################################################################################

###################################################################################################################################
## The package DESeq2 provides methods to test for differential expression by use of negative binomial generalized linear models #################################
## documentation: https://bioconductor.org/packages/3.7/bioc/vignettes/DESeq2/inst/doc/DESeq2.html
###################################################################################################################################
suppressMessages(library('DESeq2'))
library('yaml')
suppressMessages(library('gplots'))
library('RColorBrewer')
library("openxlsx")
suppressMessages(library("dplyr"))

####################################
##### Define utility functions #####
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

do.DEseq<-function(tableA,tableB,conditionA,conditionB,condA.avg, condB.avg, condition_type){
  TABLE<-rbind(tableA, tableB)
  DDS<-DESeqDataSetFromHTSeqCount(sampleTable = TABLE, directory = counts_dir, design= ~ condition)
  if (condition_type[conditionA] == 'control'){
    DDS$condition <- relevel(DDS$condition, ref = conditionA)
  } else if (condition_type[conditionB] == 'control'){
    dds$condition <- relevel(dds$condition, ref = conditionB)
  }
  DDS<-DESeq(DDS)
  DDS.res<-results(DDS, alpha=0.05)
  DDS.res$unadj.log2<-log2(condA.avg/condB.avg)
  return(DDS.res)
}

updatePairwiseDataFrame <- function(df, res, col.basename){
  df[paste0(col.basename,'.log2.deseq')] <- res$log2FoldChange
  df[paste0(col.basename,'.log2')] <- res$unadj.log2
  df[paste0(col.basename,'.log2.pval')] <- res$padj
  return(df)
}

plotSpearmenHeatmap <- function(df, file, title){
  c <- cor(df, method="spearman", use="complete.obs")
  pdf(file=file)
  mypalette<-brewer.pal(11,"Spectral")
  morecols<-colorRampPalette(mypalette)
  heatmap.2(c, main=title,col=rev(morecols(50)),trace="none", margins=c(17,6),denscol='black')
  dev.off()
}

plotPairwiseCountData <- function(df, file){
  df<-as.data.frame(df)
  reps.list<-colnames(df)
  n_reps = length(reps.list)
  pdf(file=filepath,2*n_reps,2*n_reps)
  par(mfrow=(c(n_reps,n_reps)))
  for(i in 1:n_reps){
    x.data<-df[,reps.list[i]]
    for(j in 1:n_reps){
      if(j >=i){
        
        y.data<-df[,reps.list[j]]
        plot((x.data),(y.data),xlab=reps.list[i],ylab=reps.list[j],xlim=c(0,110000),ylim=c(0,110000)) #you may want to change the x and y limits
        r2<-round(summary(lm(y.data ~ x.data))$r.squared,4)
        text(25000,100000,paste('r2=',r2)) #you may want to change where the rqsured value is recorded
      }
      else{plot(1,1,type='n',bty='n',xlab='',ylab='',axes=F)}
    }
  }
  dev.off()
}

joinChromosome <- function(df, c_elegans_annots){
  df <- as.data.frame(df)
  df$gene.name <- rownames(df)
  merged <- merge(df, c_elegans_annots, by.x="gene.name", by.y="Sequence.Name.(Gene)")
  merged <- filter(merged, Chr.Name != "MtDNA")
  return(merged)
}

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

getGenesWbId <- function(c_elegans_annots, genes){
  gene.to.wbid<-read.table(file=c_elegans_annots,header=F,stringsAsFactors=F)
  colnames(gene.to.wbid)<-c('gene','wbid')
  relevant_genes_ix <- match(genes, gene.to.wbid$gene)
  wbid<-gene.to.wbid$wbid[relevant_genes_ix]
}

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

readInFiles <- function(infiles){
  bam_suffix <- "_accepted_hits_counts.txt"
  fpkm_suffix <- "_cufflinks.genes.fpkm_tracking"
  max_id <- max(unlist(lapply(infiles, function(x) x$id)))
  files_by_id <- vector("list", max_id)
  conditions <- vector("list", max_id)
  types <- vector("list", max_id)
  
  for (ele in infiles){
    id  <- ele$id
    if (is.null(files_by_id[[id]])){
      files_by_id[[id]] <- list(ele$fastq)
    } else {
      files_by_id[[id]] <- list(files_by_id[[id]], ele$fastq)
    }
    conditions[[id]] <- ele$condition
    types[[id]] <- ele$type
  }
  count_files <- unlist(lapply(files_by_id, function(x)
    paste0(paste0(gsub(".fastq.gz", "", unlist(x)), collapse="_"), bam_suffix)))
  fpkm_files <- unlist(lapply(files_by_id, function(x)
    paste0(paste0(gsub(".fastq.gz", "", unlist(x)), collapse="_"), fpkm_suffix)))
  return(list(
    "count_files"=count_files,
    "fpkm_files"=fpkm_files,
    "conditions"=unlist(conditions),
    "types"=unlist(types)
  ))
}
########################
##### Start script #####

args <- commandArgs(trailingOnly = TRUE)
if (length(args)==0) {
  stop("At least one argument must be supplied (the deseq yaml config file).n", call.=FALSE)
}

## HTSeq input ##
conf <- yaml.load_file(args[1])
validateConfig(conf)

ercan_rna <- "/scratch/cgsb/ercan/rna/deseq"
scratch_dir <- file.path(ercan_rna, conf$experiment_title)
if (dir.exists(scratch_dir)){
  stop(sprintf("the experiment title '%s' already exists. Change this title or remove the %s directory.", conf$experiment_title, scratch_dir))
}

working_dir <- dirname(args[1])
fpkm_dir <- file.path(working_dir, "fpkm")
counts_dir <- file.path(working_dir, "counts")
deseq_dir <- file.path(working_dir, conf$experiment_title)
out_dir <- file.path(deseq_dir, "results")
dir.create(out_dir,showWarnings = FALSE, recursive = TRUE)

inFiles_data <- readInFiles(conf$infiles)
sampleCondition <- inFiles_data$conditions
conditions <- levels(factor(sampleCondition))
sampleType <- inFiles_data$type
condition_type <- vector()
for (i in 1:length(sampleCondition)){
  cond <- sampleCondition[i]
  if (!(cond %in% names(condition_type))) {condition_type[cond] = sampleType[i]}
}
sampleFiles <- inFiles_data$count_files
fpkm_files <- inFiles_data$fpkm_files

sampleTable <- data.frame(sampleName=sampleFiles, fileName=sampleFiles, condition=sampleCondition)
dds <- DESeqDataSetFromHTSeqCount(sampleTable = sampleTable, directory = counts_dir, design= ~ condition)

## filter out low/no count genes ## 
# nrow(dds)
# [1] 20389
filt.dds <- dds[ rowSums(counts(dds)) > 1, ]
# Example
# nrow(filt.dds)
# [1] 19841
# That is, 548 genes had no reads across all 40 samples. 

## get normalized count values ##
#
dds<-DESeq(dds)
normalized.count.data<-(assays(dds)[["mu"]])
colnames(normalized.count.data)<-getColNames(sampleCondition)
normalized.count.data <- as.data.frame(normalized.count.data)

c_elegans_annots <- getCelegansAnnotations(conf$c_elegans_annots)

## make heatmap plot of FPKM values of replicates based on spearman correlation
DATANAME <- paste0(conditions, collapse = 'vs')
thresholded_fpkm_data <- make_fpkm_df(fpkm_dir, fpkm_files, sampleCondition, to_threshold = TRUE)
merged <- merge(thresholded_fpkm_data, c_elegans_annots, by.x="gene_id", by.y="Sequence.Name.(Gene)")
merged <- filter(merged, Chr.Name != "MtDNA")
not_for_plot <- c('Chr.Name', 'gene_id', 'Gene.WB.ID')
for (chr in unique(merged$Chr.Name)){
  filepath <- file.path(out_dir, paste0(DATANAME,'.heatmap.spearman.thresholded.fpkm.', chr, '.pdf'))
  title <- paste0("Correlations of FPKM of ", chr, " genes")
  
  df<-filter(merged, Chr.Name == chr)
  df <- df[, !(names(df) %in% not_for_plot)]
  plotSpearmenHeatmap(df, filepath, title)
}
# plot all chromosomes together
filepath <- file.path(out_dir, paste0(DATANAME,'.heatmap.spearman.thresholded.fpkm.all.pdf'))
title <- paste0("Correlations of FPKM values")
plotSpearmenHeatmap(merged[, !(names(merged) %in% not_for_plot)], filepath, title)

## get average count data table ##
conditions_avg <- list()
for (cond in conditions){
  conditions_avg[[cond]] <- rowMeans(normalized.count.data[,sampleCondition == cond, drop=FALSE])
}

genes <- names(conditions_avg[[1]])
wbid <- getGenesWbId(conf$c_elegans_wbid_to_gene, genes)


avg.count.data<-data.frame(wbid=wbid, conditions_avg)
filepath <- file.path(out_dir,'avg.count.data.txt')
write.table(avg.count.data,file=filepath,row.names=T,col.names=T,quote=F,sep='\t')

#######################################
## plot replicates pairwise with Rsquared values ##
#######################################
filepath <- file.path(out_dir, paste0(DATANAME, '.replicates.counts.vs.counts.Rsq.pdf'))
plotPairwiseCountData(normalized.count.data, filepath)

########################
## do pairwise comparisons ##
########################
condTables = list()
for (cond in conditions){
  condTables[[cond]] = sampleTable[sampleTable['condition']==cond,]
}

n_conditions <- length(conditions)
pairwise_res_df <- data.frame(row.names = genes)
for (i in seq(n_conditions-1)){
  for (j in seq(i+1,n_conditions)){
    deseq.df <- do.DEseq(condTables[[i]], condTables[[j]],
                        conditions[[i]], conditions[[j]], 
                        conditions_avg[[i]], conditions_avg[[j]], condition_type)
    
    basename <- paste0(conditions[[i]],'vs',conditions[[j]])
    filepath <- file.path(out_dir, paste0(basename,'.deseq.txt'))
    write.table(deseq.df,file=filepath,row.names=T,col.names=T,quote=F,sep='\t')
    pairwise_res_df <- updatePairwiseDataFrame(pairwise_res_df, deseq.df, basename)
    
    filepath <- file.path(out_dir, paste0(basename,'.deseq.boxplot.by.chromosome.pdf'))
    boxPlotDESeqByChromosome(deseq.df, filepath, c_elegans_annots, "Log Fold Change By Chromosome")
    
    filepath <- file.path(out_dir, paste0(basename,'.deseq.scatterplot.pdf'))
    scatterPlotDeseq(deseq.df, c_elegans_annots, filepath, "Log Fold change vs expression")
  }
}

#### Make a table with all of the relevant data ####
genes <- names(conditions_avg[[1]])
wbid <- getGenesWbId(conf$c_elegans_wbid_to_gene, genes)
samplemeans.df<-as.data.frame(conditions_avg)

summary.data.df<-as.data.frame(cbind(wbid,samplemeans.df, pairwise_res_df))
filepath <- file.path(out_dir, 'deseq.summaryoverview.txt')
write.table(summary.data.df,file=filepath,row.names=T,col.names=T,quote=F,sep='\t')

# todo: organize folder - should have conf in top dir then a counts folder, a fpkm folder and a output folder 
file.copy(fpkm_dir, deseq_dir, recursive = TRUE)
file.copy(counts_dir, deseq_dir, recursive = TRUE)
file.copy(args[1], deseq_dir)
file.copy(deseq_dir, ercan_rna, recursive = TRUE)