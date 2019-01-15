# Title: 1_Differential_Expression
# Purpose: Identify differentially expressed lncRNAs and mRNAs from the microarray at the 4 hour IL1B time point
# Packages required: None
# Input files: "lncRNA_microarray_results" and "mRNA_microarray_results"
# Output files: "lncRNA_IL1B_differentiallyexpressed" and "mRNA_IL1B_differentiallyexpressed"

#---------------------------------------------------------------------------------------------------------------------------
# Opening up expression data for lncRNAs in the microarray file. Note only using the NS and 4hr IL1B time points, leaving out the 24hr IL1B time point for simplicity. 

lncRNA_IL1B_all <- read.csv("Path/lncRNA_microarray_results.csv")[,1:11]

# Getting rid off incomplete rows
lncRNA_IL1B_all<- lncRNA_IL1B_all[complete.cases(lncRNA_IL1B_all),]

# Removing duplicated lncRNAs
lncRNA_IL1B_all<-lncRNA_IL1B_all[!duplicated(lncRNA_IL1B_all$lncRNA),]

#---------------------------------------------------------------------------------------------------------------------------
# Converting normalized expression values out of log2 to allow for averages to be taken

lncRNA_IL1B_all[,6:11]<- mapply(function(x) 2^(as.numeric(x)),lncRNA_IL1B_all[,6:11])

#---------------------------------------------------------------------------------------------------------------------------
# Creating columns where to record average NS and 4hr IL1B treatment values, the fold change in expression, and p-values

average_NS<-vector(length=nrow(lncRNA_IL1B_all))
average_4hrIL1B<-vector(length=nrow(lncRNA_IL1B_all))
fold_change <-vector(length=nrow(lncRNA_IL1B_all))
pvalue_4tx <-vector(length=nrow(lncRNA_IL1B_all))

# Appending columns to the original dataframe

lncRNA_IL1B_all<-cbind(lncRNA_IL1B_all, average_NS, average_4hrIL1B, fold_change, pvalue_4tx)

#---------------------------------------------------------------------------------------------------------------------------
# Calculate the average values for NS and 4hr IL1B

lncRNA_IL1B_all$average_NS<-sapply(1:nrow(lncRNA_IL1B_all), function(i) mean(as.numeric(lncRNA_IL1B_all[i,6:8])))

lncRNA_IL1B_all$average_4hrIL1B<-sapply(1:nrow(lncRNA_IL1B_all), function(i) mean(as.numeric(lncRNA_IL1B_all[i,9:11])))

#---------------------------------------------------------------------------------------------------------------------------
# Calculate the fold change between NS and 4hr IL1B

lncRNA_IL1B_all$fold_change<-sapply(1:nrow(lncRNA_IL1B_all), function(i) as.numeric(lncRNA_IL1B_all$average_4hrIL1B[i]/as.numeric(lncRNA_IL1B_all$average_NS[i])))
  
#---------------------------------------------------------------------------------------------------------------------------
# Perform student t-test for the NS vs 4hr IL1B

lncRNA_IL1B_all$pvalue_4tx<-sapply(1:nrow(lncRNA_IL1B_all), function(i) t.test(as.numeric(lncRNA_IL1B_all[i,6:8]), as.numeric(lncRNA_IL1B_all[i,9:11]), paired=TRUE)$p.value )

#---------------------------------------------------------------------------------------------------------------------------
# Reverting expression values back to log2 for NS, 4hr IL1B, averages, and the fold change

lncRNA_IL1B_all[,6:14]<- mapply(function(x) log2(as.numeric(x)), (lncRNA_IL1B_all[,6:14]))

#---------------------------------------------------------------------------------------------------------------------------
# Retrieving differentially expressed lncRNAs

lncRNA_IL1B_differentiallyexpressed<-lncRNA_IL1B_all[(((lncRNA_IL1B_all$fold_change>=1)|(lncRNA_IL1B_all$fold_change<=-1))&(lncRNA_IL1B_all$pvalue_4tx<=0.05)),]

# Save "lncRNA_IL1B_differentiallyexpressed" for further analysis


#---------------------------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------------------------

# Repeating the exact same procedure for mRNA microarray expression set.
# The only difference is that the mRNA expression set has the mRNA Gene Names in the dataframe
# which the lncRNA expression set does not.

#---------------------------------------------------------------------------------------------------------------------------
# Opening up the file
mRNA_IL1B_all <- read.csv("Path/mRNA_microarray_results.csv")[,1:12]

# Getting rid off incomplete rows
mRNA_IL1B_all<- mRNA_IL1B_all[complete.cases(mRNA_IL1B_all),]

# Removing duplicated mRNAs
mRNA_IL1B_all<-mRNA_IL1B_all[!duplicated(mRNA_IL1B_all$mRNA),]

#---------------------------------------------------------------------------------------------------------------------------
# Converting normalized expression values out of log2 to allow for averages to be taken

mRNA_IL1B_all[,7:12]<- mapply(function(x) 2^(as.numeric(x)),mRNA_IL1B_all[,7:12])

#---------------------------------------------------------------------------------------------------------------------------
# Creating columns where to record average NS and 4hr IL1B treatment values, the fold change in expression, and p-values

average_NS<-vector(length=nrow(mRNA_IL1B_all))
average_4hrIL1B<-vector(length=nrow(mRNA_IL1B_all))
fold_change <-vector(length=nrow(mRNA_IL1B_all))
pvalue_4tx <-vector(length=nrow(mRNA_IL1B_all))

# Appending columns to the original dataframe

mRNA_IL1B_all<-cbind(mRNA_IL1B_all, average_NS, average_4hrIL1B, fold_change, pvalue_4tx)

#---------------------------------------------------------------------------------------------------------------------------
# Calculate the average values for NS and 4hr IL1B

mRNA_IL1B_all$average_NS<-sapply(1:nrow(mRNA_IL1B_all), function(i) mean(as.numeric(mRNA_IL1B_all[i,7:9])))

mRNA_IL1B_all$average_4hrIL1B<-sapply(1:nrow(mRNA_IL1B_all), function(i) mean(as.numeric(mRNA_IL1B_all[i,10:12])))

#---------------------------------------------------------------------------------------------------------------------------
# Calculate the fold change between NS and 4hr IL1B

mRNA_IL1B_all$fold_change<-sapply(1:nrow(mRNA_IL1B_all), function(i) as.numeric(mRNA_IL1B_all$average_4hrIL1B[i]/as.numeric(mRNA_IL1B_all$average_NS[i])))

#---------------------------------------------------------------------------------------------------------------------------
# Perform student t-test for the 4hr IL1B vs NS 

mRNA_IL1B_all$pvalue_4tx<-sapply(1:nrow(mRNA_IL1B_all), function(i) t.test(as.numeric(mRNA_IL1B_all[i,7:9]), as.numeric(mRNA_IL1B_all[i,10:12]), paired=TRUE)$p.value )

#---------------------------------------------------------------------------------------------------------------------------
# Reverting expression values back to log2 for NS, 4hr IL1B, averages, and the fold change

mRNA_IL1B_all[,7:15]<- mapply(function(x) log2(as.numeric(x)), (mRNA_IL1B_all[,7:15]))

#---------------------------------------------------------------------------------------------------------------------------
# Retrieving differentially expressed mRNAs

mRNA_IL1B_differentiallyexpressed<-mRNA_IL1B_all[(((mRNA_IL1B_all$fold_change>=1)|(mRNA_IL1B_all$fold_change<=-1))&(mRNA_IL1B_all$pvalue_4tx<=0.05)),]

# Save "mRNA_IL1B_differentiallyexpressed" for further analysis