# Title: 1_Differential_Expression
# Purpose: Identify differentially expressed lncRNAs and mRNAs from the microarray at the 4 hour IL1B time point
# Packages required: None
# Input files: "lncRNA_microarray_results" and "mRNA_microarray_results"
# Output files: "lncRNA_IL1B_differentiallyexpressed" and "mRNA_IL1B_differentiallyexpressed"
# "lncRNA_IL1B_differentiallyexpressed_24tx" and "mRNA_IL1B_differentiallyexpressed_24tx"
#---------------------------------------------------------------------------------------------------------------------------
# Opening up expression data for lncRNAs in the microarray file. 

lncRNA_IL1B_all <- read.csv("Path/Table_S1_lncRNA_microarray_results.csv")[,1:14]

# Getting rid off incomplete rows
lncRNA_IL1B_all<- lncRNA_IL1B_all[complete.cases(lncRNA_IL1B_all),]

# Removing duplicated lncRNAs
lncRNA_IL1B_all<-lncRNA_IL1B_all[!duplicated(lncRNA_IL1B_all$lncRNA),]

#---------------------------------------------------------------------------------------------------------------------------
# Converting normalized expression values out of log2 to allow for averages to be taken

lncRNA_IL1B_all[,6:14]<- mapply(function(x) 2^(as.numeric(x)),lncRNA_IL1B_all[,6:14])

#---------------------------------------------------------------------------------------------------------------------------
# Creating columns where to record average NS, 4hr and 24hr IL1B treatment values, the fold change in expression, and p-values

average_NS<-vector(length=nrow(lncRNA_IL1B_all))
average_4tx<-vector(length=nrow(lncRNA_IL1B_all))
average_24tx<-vector(length=nrow(lncRNA_IL1B_all))

fold_change <-vector(length=nrow(lncRNA_IL1B_all))
fold_change_24tx <-vector(length=nrow(lncRNA_IL1B_all))

pvalue_4tx <-vector(length=nrow(lncRNA_IL1B_all))
pvalue_24tx <-vector(length=nrow(lncRNA_IL1B_all))

# Appending columns to the original dataframe

lncRNA_IL1B_all_24tx<-cbind(lncRNA_IL1B_all[,c(1:8,12:14)], average_NS, average_24tx, fold_change_24tx, pvalue_24tx)

lncRNA_IL1B_all<-cbind(lncRNA_IL1B_all[,1:11], average_NS, average_4tx, fold_change, pvalue_4tx)

#---------------------------------------------------------------------------------------------------------------------------
# Calculate the average values for NS, 4hr IL1B, and 24hr IL1B

lncRNA_IL1B_all$average_NS<-sapply(1:nrow(lncRNA_IL1B_all), function(i) mean(as.numeric(lncRNA_IL1B_all[i,6:8])))

lncRNA_IL1B_all$average_4tx<-sapply(1:nrow(lncRNA_IL1B_all), function(i) mean(as.numeric(lncRNA_IL1B_all[i,9:11])))

lncRNA_IL1B_all_24tx$average_NS<-sapply(1:nrow(lncRNA_IL1B_all_24tx), function(i) mean(as.numeric(lncRNA_IL1B_all_24tx[i,6:8])))

lncRNA_IL1B_all_24tx$average_24tx<-sapply(1:nrow(lncRNA_IL1B_all_24tx), function(i) mean(as.numeric(lncRNA_IL1B_all_24tx[i,9:11])))

#---------------------------------------------------------------------------------------------------------------------------
# Calculate the fold change between NS and 4hr IL1B or NS and 24hr IL1B

lncRNA_IL1B_all$fold_change<-sapply(1:nrow(lncRNA_IL1B_all), function(i) as.numeric(lncRNA_IL1B_all$average_4tx[i]/as.numeric(lncRNA_IL1B_all$average_NS[i])))

lncRNA_IL1B_all_24tx$fold_change_24tx<-sapply(1:nrow(lncRNA_IL1B_all_24tx), function(i) as.numeric(lncRNA_IL1B_all_24tx$average_24tx[i]/as.numeric(lncRNA_IL1B_all_24tx$average_NS[i])))

#---------------------------------------------------------------------------------------------------------------------------
# Perform student t-test for the NS vs 4hr IL1B or NS vs 24hr IL1B

lncRNA_IL1B_all$pvalue_4tx<-sapply(1:nrow(lncRNA_IL1B_all), function(i) t.test(as.numeric(lncRNA_IL1B_all[i,6:8]), as.numeric(lncRNA_IL1B_all[i,9:11]), paired=TRUE)$p.value )

lncRNA_IL1B_all_24tx$pvalue_24tx<-sapply(1:nrow(lncRNA_IL1B_all_24tx), function(i) t.test(as.numeric(lncRNA_IL1B_all_24tx[i,6:8]), as.numeric(lncRNA_IL1B_all_24tx[i,9:11]), paired=TRUE)$p.value )

#---------------------------------------------------------------------------------------------------------------------------
# Reverting expression values back to log2 for NS, 4hr IL1B, 24hr IL1B, averages, and the fold change

lncRNA_IL1B_all[,6:14]<- mapply(function(x) log2(as.numeric(x)), (lncRNA_IL1B_all[,6:14]))


lncRNA_IL1B_all_24tx[,6:14]<- mapply(function(x) log2(as.numeric(x)), (lncRNA_IL1B_all_24tx[,6:14]))
------------------------------------------------------------------------------------------------------------------------
# Retrieving differentially expressed lncRNAs

  
lncRNA_IL1B_differentiallyexpressed<-lncRNA_IL1B_all[(((lncRNA_IL1B_all$fold_change>=1)|(lncRNA_IL1B_all$fold_change<=-1))&(lncRNA_IL1B_all$pvalue_4tx<=0.05)),]


lncRNA_IL1B_differentiallyexpressed_24tx<-lncRNA_IL1B_all_24tx[(((lncRNA_IL1B_all_24tx$fold_change_24tx>=1)|(lncRNA_IL1B_all_24tx$fold_change_24tx<=-1))&(lncRNA_IL1B_all_24tx$pvalue_24tx<=0.05)),]

# Save "lncRNA_IL1B_differentiallyexpressed" for further analysis
# Save "lncRNA_IL1B_differentiallyexpressed_24tx" for reference

#---------------------------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------------------------

# Repeating the exact same procedure for mRNA microarray expression set.
# The only difference is that the mRNA expression set has the mRNA Gene Names in the dataframe
# which the lncRNA expression set does not.

#---------------------------------------------------------------------------------------------------------------------------
# Opening up the file
mRNA_IL1B_all <- read.csv("Path/Table_S2_mRNA_microarray_results.csv")[,1:15]

# Getting rid off incomplete rows
mRNA_IL1B_all<- mRNA_IL1B_all[complete.cases(mRNA_IL1B_all),]

# Removing duplicated mRNAs
mRNA_IL1B_all<-mRNA_IL1B_all[!duplicated(mRNA_IL1B_all$mRNA),]

#---------------------------------------------------------------------------------------------------------------------------
# Converting normalized expression values out of log2 to allow for averages to be taken

mRNA_IL1B_all[,7:15]<- mapply(function(x) 2^(as.numeric(x)),mRNA_IL1B_all[,7:15])

#---------------------------------------------------------------------------------------------------------------------------
# Creating columns where to record average NS, 4hr, 24hr IL1B treatment values, the fold change in expression, and p-values

average_NS<-vector(length=nrow(mRNA_IL1B_all))
average_4tx<-vector(length=nrow(mRNA_IL1B_all))
average_24tx<-vector(length=nrow(mRNA_IL1B_all))

fold_change <-vector(length=nrow(mRNA_IL1B_all))
fold_change_24tx <-vector(length=nrow(mRNA_IL1B_all))

pvalue_4tx <-vector(length=nrow(mRNA_IL1B_all))
pvalue_24tx <-vector(length=nrow(mRNA_IL1B_all))

# Appending columns to the original dataframe

mRNA_IL1B_all_24tx<-cbind(mRNA_IL1B_all[,c(1:9,13:15)], average_NS, average_24tx, fold_change_24tx, pvalue_24tx)

mRNA_IL1B_all<-cbind(mRNA_IL1B_all[,1:12], average_NS, average_4tx, fold_change, pvalue_4tx)

#---------------------------------------------------------------------------------------------------------------------------
# Calculate the average values for NS and 4hr IL1B/ NS and 24hr IL1B

mRNA_IL1B_all$average_NS<-sapply(1:nrow(mRNA_IL1B_all), function(i) mean(as.numeric(mRNA_IL1B_all[i,7:9])))
mRNA_IL1B_all$average_4tx<-sapply(1:nrow(mRNA_IL1B_all), function(i) mean(as.numeric(mRNA_IL1B_all[i,10:12])))


mRNA_IL1B_all_24tx$average_NS<-sapply(1:nrow(mRNA_IL1B_all_24tx), function(i) mean(as.numeric(mRNA_IL1B_all_24tx[i,7:9])))
mRNA_IL1B_all_24tx$average_24tx<-sapply(1:nrow(mRNA_IL1B_all_24tx), function(i) mean(as.numeric(mRNA_IL1B_all_24tx[i,10:12])))


#---------------------------------------------------------------------------------------------------------------------------
# Calculate the fold change between NS and 4hr IL1B/ NS and 24hr IL1B

mRNA_IL1B_all$fold_change<-sapply(1:nrow(mRNA_IL1B_all), function(i) as.numeric(mRNA_IL1B_all$average_4tx[i]/as.numeric(mRNA_IL1B_all$average_NS[i])))

mRNA_IL1B_all_24tx$fold_change_24tx<-sapply(1:nrow(mRNA_IL1B_all_24tx), function(i) as.numeric(mRNA_IL1B_all_24tx$average_24tx[i]/as.numeric(mRNA_IL1B_all_24tx$average_NS[i])))

#---------------------------------------------------------------------------------------------------------------------------
# Perform student t-test for the 4hr IL1B vs NS / 24hr IL1B vs NS

mRNA_IL1B_all$pvalue_4tx<-sapply(1:nrow(mRNA_IL1B_all), function(i) t.test(as.numeric(mRNA_IL1B_all[i,7:9]), as.numeric(mRNA_IL1B_all[i,10:12]), paired=TRUE)$p.value )

mRNA_IL1B_all_24tx$pvalue_24tx<-sapply(1:nrow(mRNA_IL1B_all_24tx), function(i) t.test(as.numeric(mRNA_IL1B_all_24tx[i,7:9]), as.numeric(mRNA_IL1B_all_24tx[i,10:12]), paired=TRUE)$p.value )

#---------------------------------------------------------------------------------------------------------------------------
# Reverting expression values back to log2 for NS, 4hr IL1B, 24hr IL1B, averages, and the fold change

mRNA_IL1B_all[,7:15]<- mapply(function(x) log2(as.numeric(x)), (mRNA_IL1B_all[,7:15]))

mRNA_IL1B_all_24tx[,7:15]<- mapply(function(x) log2(as.numeric(x)), (mRNA_IL1B_all_24tx[,7:15]))

#---------------------------------------------------------------------------------------------------------------------------
# Retrieving differentially expressed mRNAs

mRNA_IL1B_differentiallyexpressed<-mRNA_IL1B_all[(((mRNA_IL1B_all$fold_change>=1)|(mRNA_IL1B_all$fold_change<=-1))&(mRNA_IL1B_all$pvalue_4tx<=0.05)),]

mRNA_IL1B_differentiallyexpressed_24tx<-mRNA_IL1B_all_24tx[(((mRNA_IL1B_all_24tx$fold_change_24tx>=1)|(mRNA_IL1B_all_24tx$fold_change_24tx<=-1))&(mRNA_IL1B_all_24tx$pvalue_24tx<=0.05)),]

# Save "mRNA_IL1B_differentiallyexpressed" for further analysis
# Save "mRNA_IL1B_differentiallyexpressed_24tx" for future reference
