################################################################################
## Takes the count and fpkm data and then normalizes it to generate summary   ##
## csv files that contain CPM, FPKM And TPMs for every replicate and means    ##
##                                                                            ##
## usage: Rscript create_summary_outputs.R WORKING_DIR                        ##
##                                                                            ##
## positional arguments:                                                      ##
##  WORKING DIR                                                               ##
##  The working directory is essential. This will be direct from pipeline.    ##
##                                                                            ##
################################################################################


args <- commandArgs(trailingOnly = TRUE)
setwd(paste0(args[1],'/counts/'))

#Identify the count files and read them into a list
pathToFiles <- as.list(paste0(getwd(),'/',list.files()))
counts <- lapply(pathToFiles, read.table, header = F, stringsAsFactors=F)
#Extract the names of each replicate
my_file_names<-list.files()

#Function that will extract sepecifically the names and conditions
name_extractor<-function(file_names){
name_breakdown<-strsplit(file_names,'[_.]')
if(length(name_breakdown[[1]])==5){
my_rep<-name_breakdown[[1]][4]
my_strain<-name_breakdown[[1]][2]
my_stage<-name_breakdown[[1]][3]
my_library<-name_breakdown[[1]][1]
}else{
my_rep<-name_breakdown[[1]][4]
my_strain<-name_breakdown[[1]][2]
my_stage<-name_breakdown[[1]][3]
library_IDs<-seq(1,(length(name_breakdown[[1]])-1),by=4)
my_library<-paste0(name_breakdown[[1]][library_IDs], collapse="-")
}
condition<-paste(my_strain,my_stage,sep='_')
final_name<-paste(my_strain,my_stage,my_rep,my_library,sep='_')
return(list(condition,final_name))
}

#apply function to the names
my_result <- apply(data.matrix(my_file_names),1,name_extractor)
conditionNames<-sapply(my_result, "[[", 1)
sampleNames<-sapply(my_result, "[[", 2)

#Add the names
for(i in 1:length(counts)){
  colnames(counts[[i]])[2] <- sampleNames[i]
}

#What are conditions?
conditions<-unique(conditionNames)

#Flatten datasets into a single dataframe
counts <- Reduce(merge, counts)
colnames(counts)[1]<-'Feature'
counts <- counts[-(1:5),]

#Calculate counts per million
cpm<- apply(counts[,-1],2, function(x) (x/sum(x))*1000000)
cpm<-as.data.frame(cbind(counts[,1],cpm), stringsAsFactors=F)
colnames(cpm)[1]<-'Feature'

#Calculate means and SD for CPM
for (i in 1:length(conditions)){
my_ind<-which(conditionNames==conditions[i])
mean_cpm<-apply(data.matrix(cpm[,(my_ind+1)]),1,mean)
sd_cpm<-apply(data.matrix(cpm[,(my_ind+1)]),1,sd)
cpm_summary<-cbind(mean_cpm,sd_cpm)
colnames(cpm_summary)<-paste0(c('cpm_mean_','cpm_sd_'),conditions[i])
assign(paste0('cpm_summary_',conditions[i]),cpm_summary)
}

#Identify the FPKM files from cufflinks and read them into a list
setwd('../cufflinks/')
pathToFiles <- as.list(paste0(getwd(),'/',list.files('.', recursive=T, pattern='genes.fpkm_tracking')))
fpkm <- lapply(pathToFiles, read.table, header = T, stringsAsFactors=F, colClasses=c('character', rep("NULL",5 ),'character', rep("NULL",2 ),'numeric', rep("NULL",3 )))

for(i in 1:length(fpkm)){
  colnames(fpkm[[i]])[3] <- sampleNames[i]
}

fpkm <- Reduce(merge, fpkm)
colnames(fpkm)[1]<-'Feature'

#Calculate means and SD for FPKM
for (i in 1:length(conditions)){
my_ind<-which(conditionNames==conditions[i])
mean_fpkm<-apply(data.matrix(fpkm[,(my_ind+2)]),1,mean)
sd_fpkm<-apply(data.matrix(fpkm[,(my_ind+2)]),1,sd)
fpkm_summary<-cbind(mean_fpkm,sd_fpkm)
colnames(fpkm_summary)<-paste0(c('fpkm_mean_','fpkm_sd_'),conditions[i])
assign(paste0('fpkm_summary_',conditions[i]),fpkm_summary)
}

#Calculate TPM values
sample_sums<-apply(data.matrix(fpkm[,-c(1:2)]),2,sum)
tpm<-t(t(fpkm[,-c(1:2)])/sample_sums)*10^6
tpm<-as.data.frame(cbind(fpkm[,c(1:2)],tpm), stringsAsFactors=F)

#Calculate means and SD for TPM
for (i in 1:length(conditions)){
my_ind<-which(conditionNames==conditions[i])
mean_tpm<-apply(data.matrix(tpm[,(my_ind+2)]),1,mean)
sd_tpm<-apply(data.matrix(tpm[,(my_ind+2)]),1,sd)
tpm_summary<-cbind(mean_tpm,sd_tpm)
colnames(tpm_summary)<-paste0(c('tpm_mean_','tpm_sd_'),conditions[i])
assign(paste0('tpm_summary_',conditions[i]),tpm_summary)
}

#Create summary output files and write them
setwd('../')
dir.create('summary')
for (i in 1:length(conditions)){
summary_cpm_tpm_fpkm<-fpkm[,c(1:2)]
my_cpm<-get(paste0('cpm_summary_',conditions[i]))
my_cpm<-as.data.frame(cbind(cpm[,1], my_cpm), stringsAsFactors=F)
colnames(my_cpm)[1]<-'Feature'
my_tpm<-get(paste0('tpm_summary_',conditions[i]))
my_fpkm<-get(paste0('fpkm_summary_',conditions[i]))
summary_cpm_tpm_fpkm<-cbind(summary_cpm_tpm_fpkm,my_tpm,my_fpkm)
write.csv(merge(summary_cpm_tpm_fpkm,my_cpm, all=TRUE), file=paste0(args[1],'/summary/summary_cpm_tpm_fpkm',conditions[i],'.csv'))
}

#Write out the files
cpm<-merge(fpkm[,1:2],cpm)
write.csv(tpm, file=paste0(args[1],'/summary/all_tpms_',paste0(conditions,collapse='-'),'.csv'))
write.csv(fpkm, file=paste0(args[1],'/summary/all_fpkms_',paste0(conditions,collapse='-'),'.csv'))
write.csv(cpm, file=paste0(args[1],'/summary/all_cpms',paste0(conditions,collapse='-'),'.csv'))

## Make heatmap plot of CPM values of replicates based on spearman correlation

# Filter out Mitochondria genes
cpm<-cpm[-grep('MtDNA',cpm[,2]),]

chrs<-c('CHROMOSOME_I:','CHROMOSOME_II:','CHROMOSOME_III:','CHROMOSOME_IV:','CHROMOSOME_V:', 'CHROMOSOME_X:')

## Plot correlations between chr heatmaps
plotSpearmenHeatmap <- function(df, file, title){
  c <- cor(df, method="spearman", use="complete.obs")
  pdf(file=file)
  mypalette<-brewer.pal(11,"Spectral")
  morecols<-colorRampPalette(mypalette)
  heatmap.2(c, main=title,col=rev(morecols(50)),trace="none", margins=c(17,6),denscol='black')
  dev.off()
}

#Run it over each chromosome
suppressMessages(library('gplots'))
library('RColorBrewer')
for (chr in chrs){
  #Set names and titles
  filepath <- paste0(args[1],'/summary/',paste0(conditions,collapse='-'),chr,'.heatmap.spearman.thresholded.fpkm.pdf')
  title <- paste0("Correlations of CPM of", chr, "genes\n",paste0(conditions,collapse='-'))
  #pull out a single chromosome
  df<-cpm[grep(chr,cpm[,2]),]
  #plot the correlations
  plotSpearmenHeatmap(df[,-c(1:2)], filepath, title)
}
# plot all chromosomes together
filepath <- paste0(args[1],'/summary/',paste0(conditions,collapse='-'),'ALL.heatmap.spearman.thresholded.fpkm.pdf')
title <- paste0("Correlations of CPM of all genes\n",paste0(conditions,collapse='-'))
  plotSpearmenHeatmap(cpm[,-c(1:2)], filepath, title)
