

#Identify the count files and read them into a list
pathToFiles <- as.list(paste0(getwd(),'/',list.files()))
counts <- lapply(pathToFiles, read.table, header = F, stringsAsFactors=F)
#Extract the names of each replicate
my_file_names<-list.files()

#Function that will extract sepecifcally the names
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
final_name<-paste(my_strain,my_stage,my_rep,my_library,sep='_')
return(final_name)
}

#apply function to the names
sampleNames <- apply(data.matrix(my_file_names),1,name_extractor)

#Add the names
for(i in 1:length(counts)){
  colnames(counts[[i]])[2] <- sampleNames[i]
}

#Flatten datasets into a single dataframe
counts <- Reduce(merge, counts)
colnames(counts)[1]<-'Feature'
counts <- counts[-(1:5),]

#Calculate counts per million
cpm<- apply(counts[,-1],2, function(x) (x/sum(x))*1000000)
cpm<-as.data.frame(cbind(counts[,1],cpm), stringsAsFactors=F)

#Calculate counts per million
cpm<- apply(counts[,-1],2, function(x) (x/sum(x))*1000000)
cpm<-as.data.frame(cbind(counts[,1],cpm), stringsAsFactors=F)


#outputs
1)counts
2)tpm
3)rpkm
4)cpm
5)means of above and Sdev fro summary table
