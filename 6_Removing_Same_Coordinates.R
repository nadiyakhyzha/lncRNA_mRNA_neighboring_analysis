# Title: 6_Removing_Same_Coordinates
# Purpose: To remove mRNAde-lncRNA and mRNAde-lncRNAde pairs where the two transcripts
# have the same coordinates. Afterwards tabulate the number of Transcript Neighbors per mRNAde
# and number of mRNAde that have a neighboring transcript
# Packages required: None
# Input files: "original_mRNAde_lncRNA_[distance]", mRNAde_mRNAde_[distance], mRNAde_mRNA_[distance], "lncRNA_IL1B_differentiallyexpressed", "mRNA_IL1B_differentiallyexpressed",
# "mRNAde_lncRNA_transcriptional_classification", and "mRNAde_mRNA_transcriptional_classification"
# Output files: "mRNAde_lncRNAde_[distance]", "mRNAde_lncRNA_[distance]", "mRNA_lncRNA_tabulation"
# Output graphs: "mRNAlncRNA_mRNAmRNA_tabulation", "mRNA_permRNA_perlncRNA"

# -----------------------------------------------------------------------------------------------------------------------
# Loading the mRNA-lncRNA neighboring pair csv files 
# Start by setting up a list of genomic distances for which to open up the files
distance<-as.integer(c(seq(5000, 100000, 5000),125000,150000,175000,200000,
                       225000, 250000,275000, 300000, 400000, 500000))

# Opening up the mRNA-lncRNA neighboring pairs for the selected genomic distances and combining all these dataframes in one list
mRNA_lncRNA_list<-lapply(seq_along(distance), function(i)   assign(paste("mRNAde_lncRNA_",distance[i], sep=""), read.csv(paste("Path/original_mRNAde_lncRNA_", distance[i], ".csv", sep=""), stringsAsFactors=FALSE)[,2:3]))

# -----------------------------------------------------------------------------------------------------------------------
# Opening up the differentially expressed list of lncRNAs and mRNAs
lncRNA_IL1B <- read.csv("Path/lncRNA_IL1B_differentiallyexpressed.csv", stringsAsFactors=FALSE)
mRNA_IL1B<- read.csv("Path/mRNA_IL1B_differentiallyexpressed.csv", stringsAsFactors=FALSE)

# -----------------------------------------------------------------------------------------------------------------------
# Next searching through the list of differentially expressed lncRNAs and mRNAs and retrieving their expression 
# fold changes in response to IL1B. All lncRNAs that are not differentially expressed are given the value of "0"

  
mRNA_foldchange<-sapply(seq_along(mRNA_lncRNA_list), function(i) sapply(1:nrow(mRNA_lncRNA_list[[i]]), function(j) mRNA_IL1B$fold_change[mRNA_lncRNA_list[[i]]$mRNA_Symbol[j]==mRNA_IL1B$mRNA]))
lncRNA_foldchange<-sapply(seq_along(mRNA_lncRNA_list), function(i) sapply(1:nrow(mRNA_lncRNA_list[[i]]), function(j) ifelse(sum(mRNA_lncRNA_list[[i]]$lncRNA_Symbol[j]==lncRNA_IL1B$lncRNA)>=1, lncRNA_IL1B$fold_change[mRNA_lncRNA_list[[i]]$lncRNA_Symbol[j]==lncRNA_IL1B$lncRNA], 0 )))

mRNA_lncRNA_list<-lapply(seq_along(mRNA_lncRNA_list), function(i) cbind(mRNA_lncRNA_list[[i]],mRNA_foldchange[[i]], lncRNA_foldchange[[i]]))

# -----------------------------------------------------------------------------------------------------------------------
# Separating the list into mRNAs with regulated lncRNAs or not regulated lncRNAs
# To start two new lists of dataframes are created as a copy of the original list of dataframes
regulated_mRNA_lncRNA_list<-mRNA_lncRNA_list
unregulated_mRNA_lncRNA_list<-mRNA_lncRNA_list

# All of the mRNA-lncRNA pairs where the lncRNA is differentially expressed will be kept in the regulated_dataframe_list
regulated_mRNA_lncRNA_list<-lapply(1:30, function(i) regulated_mRNA_lncRNA_list[[i]][((regulated_mRNA_lncRNA_list[[i]]$`lncRNA_foldchange[[i]]`>=1)|(regulated_mRNA_lncRNA_list[[i]]$`lncRNA_foldchange[[i]]`<=-1)),])

# All of the mRNA-lncRNA pairs where the lncRNA is not differentially expressed will be kept in the unregulated_dataframe_list
unregulated_mRNA_lncRNA_list<-lapply(1:30, function(i) unregulated_mRNA_lncRNA_list[[i]][!((unregulated_mRNA_lncRNA_list[[i]]$`lncRNA_foldchange[[i]]`>=1)|(unregulated_mRNA_lncRNA_list[[i]]$`lncRNA_foldchange[[i]]`<=-1)),])

# -----------------------------------------------------------------------------------------------------------------------
# The next section is based on previously performed analysis where certain mRNA-lncRNA pairs were found
# to have the same coordinates. FOr the purposes of continuing analysis these mRNA-lncRNA pairs will
# be removed from the list
# To do so the classification lists will be opened up

mRNAde_lncRNA_classification_list<-lapply(seq_along(distance), function(i)   assign(paste("mRNAde_lncRNA_",distance[i], "classification", sep=""), read.csv(paste("Path/mRNAde_lncRNA_transcriptional_classification", distance[i], ".csv", sep=""), stringsAsFactors=FALSE)[,-1]))
mRNAde_lncRNAde_classification_list<-lapply(seq_along(distance), function(i)   assign(paste("mRNAde_lncRNAde_",distance[i], "classification", sep=""), read.csv(paste("Path/mRNAde_lncRNAde_transcriptional_classification", distance[i], ".csv", sep=""), stringsAsFactors=FALSE)[,-1]))

# Now filtering out all of the mRNA-lncRNA pairs where they have the same coordinates
mRNAde_lncRNA_list<-lapply(seq_along(distance), function(i) unregulated_mRNA_lncRNA_list[[i]][mRNAde_lncRNA_classification_list[[i]]$V7!="same coordinates",])
mRNAde_lncRNAde_list<-lapply(seq_along(distance), function(i) regulated_mRNA_lncRNA_list[[i]][mRNAde_lncRNAde_classification_list[[i]]$V7!="same coordinates",])

# -----------------------------------------------------------------------------------------------------------------------
# Saving the regulated_mRNA_lncRNA_list and unregulated_mRNA_lncRNA_list files as csv files

lapply(seq_along(mRNAde_lncRNA_list),
       function(i)(write.csv(mRNAde_lncRNA_list[[i]],
                             file =paste0('mRNAde_lncRNA_',distance[i],'.csv'))))

lapply(seq_along(mRNAde_lncRNAde_list),
       function(i)(write.csv(mRNAde_lncRNAde_list[[i]],
                             file =paste0('mRNAde_lncRNAde_',distance[i],'.csv'))))

# -----------------------------------------------------------------------------------------------------------------------
# Performing tabulations of mRNA-lncRNA/mRNA-mRNA pairs
# To begin open the mRNAde-mRNA and mRNAde-mRNAde pairs

mRNAde_mRNA_list<-lapply(seq_along(distance), function(i)   assign(paste("mRNAde_mRNA_",distance[i], sep=""), read.csv(paste("Path/mRNAde_mRNA_", distance[i], ".csv", sep=""), stringsAsFactors=FALSE)[,2:3]))
mRNAde_mRNAde_list<-lapply(seq_along(distance), function(i)   assign(paste("mRNAde_mRNAde_",distance[i], sep=""), read.csv(paste("Path/mRNAde_mRNAde_", distance[i], ".csv", sep=""), stringsAsFactors=FALSE)[,2:3]))

# Tabulation of unique # of mRNAde that have a neighboring transcript for a specified genomic distance

mRNAde_lncRNA_tab<-sapply(seq_along(distance), function(i) length(unique(mRNAde_lncRNA_list[[i]]$mRNA_Symbol)))
mRNAde_lncRNAde_tab<-sapply(seq_along(distance), function(i) length(unique(mRNAde_lncRNAde_list[[i]]$mRNA_Symbol)))
mRNAde_mRNA_tab<-sapply(seq_along(distance), function(i) length(unique(mRNAde_mRNA_list[[i]]$mRNA_Symbol)))
mRNAde_mRNAde_tab<-sapply(seq_along(distance), function(i) length(unique(mRNAde_mRNAde_list[[i]]$mRNA_Symbol)))

# Combind the previous tabulations into one dataframe
mRNA_lncRNA_tabulation<-cbind(mRNAde_lncRNA_tab, mRNAde_lncRNAde_tab, mRNAde_mRNA_tab, mRNAde_mRNAde_tab)
rownames(mRNA_lncRNA_tabulation)<-distance
# Save this file

# -----------------------------------------------------------------------------------------------------------------------
# Grooming the dataframe for plotting
tabulation<-c(mRNA_lncRNA_tabulation[10,1], mRNA_lncRNA_tabulation[10,2], mRNA_lncRNA_tabulation[10,3], mRNA_lncRNA_tabulation[10,4],
              mRNA_lncRNA_tabulation[22,1], mRNA_lncRNA_tabulation[22,2], mRNA_lncRNA_tabulation[22,3], mRNA_lncRNA_tabulation[22,4],
              mRNA_lncRNA_tabulation[28,1], mRNA_lncRNA_tabulation[28,2], mRNA_lncRNA_tabulation[28,3], mRNA_lncRNA_tabulation[28,4])

# Plot
pdf("mRNAlncRNA_mRNAmRNA_tabulation.pdf", height = 5, width = 7, paper="special")
barplot(unlist(tabulation), col=c("white", "red", "light blue", "blue", 
                                  "white", "red", "light blue", "blue", 
                                  "white", "red", "light blue", "blue"), 
        ylab="# mRNA with neighboring transcript", space=c(0.25,0.25,0.25,0.25,1.5,0.25,0.25,0.25,1.5,0.25,0.25,0.25))
dev.off()

# -----------------------------------------------------------------------------------------------------------------------
# Another way of tabulating the neighboring pairs.
# Here will calculate the number of lncRNAs or mRNAs that each mRNAde has neighboring to it
# for a specified genomic distance

# Do this by creating a table where the number of times an mRNAde is present is tabulated. Afterwards calculate
# the median number of transcript pairs that each mRNAde has.

mRNAde_lncRNA_50000<-data.frame(table(mRNAde_lncRNA_list[[10]]$mRNA_Symbol))
median(mRNAde_lncRNA_50000[,2])
#2
mRNAde_lncRNA_150000<-data.frame(table(mRNAde_lncRNA_list[[22]]$mRNA_Symbol))
median(mRNAde_lncRNA_150000[,2])
#3
mRNAde_lncRNA_300000<-data.frame(table(mRNAde_lncRNA_list[[28]]$mRNA_Symbol))
median(mRNAde_lncRNA_300000[,2])
#6

mRNAde_lncRNAde_50000<-data.frame(table(mRNAde_lncRNAde_list[[10]]$mRNA_Symbol))
median(mRNAde_lncRNAde_50000[,2])
#1

mRNAde_lncRNAde_150000<-data.frame(table(mRNAde_lncRNAde_list[[22]]$mRNA_Symbol))
median(mRNAde_lncRNAde_150000[,2])
#1

mRNAde_lncRNAde_300000<-data.frame(table(mRNAde_lncRNAde_list[[28]]$mRNA_Symbol))
median(mRNAde_lncRNAde_300000[,2])
#1




mRNAde_mRNA_50000<-data.frame(table(mRNAde_mRNA_list[[10]]$mRNA_Symbol))
median(mRNAde_mRNA_50000[,2])
#2
mRNAde_mRNA_150000<-data.frame(table(mRNAde_mRNA_list[[22]]$mRNA_Symbol))
median(mRNAde_mRNA_150000[,2])
#4
mRNAde_mRNA_300000<-data.frame(table(mRNAde_mRNA_list[[28]]$mRNA_Symbol))
median(mRNAde_mRNA_300000[,2])
#7

mRNAde_mRNAde_50000<-data.frame(table(mRNAde_mRNAde_list[[10]]$mRNA_Symbol))
median(mRNAde_mRNAde_50000[,2])
#1

mRNAde_mRNAde_150000<-data.frame(table(mRNAde_mRNAde_list[[22]]$mRNA_Symbol))
median(mRNAde_mRNAde_150000[,2])
#1

mRNAde_mRNAde_300000<-data.frame(table(mRNAde_mRNAde_list[[28]]$mRNA_Symbol))
median(mRNAde_mRNAde_300000[,2])
#1

# -----------------------------------------------------------------------------------------------------------------------
# Now plotting the number of neighboring transcripts per mRNAde

pdf("mRNA_permRNA_perlncRNA.pdf", height = 5, width = 7, paper="special")
boxplot(mRNAde_lncRNA_50000[,2], mRNAde_lncRNAde_50000[,2], mRNAde_mRNA_50000[,2], mRNAde_mRNAde_50000[,2], 
        mRNAde_lncRNA_150000[,2], mRNAde_lncRNAde_150000[,2], mRNAde_mRNA_150000[,2], mRNAde_mRNAde_150000[,2], 
        mRNAde_lncRNA_300000[,2], mRNAde_lncRNAde_300000[,2], mRNAde_mRNA_300000[,2], mRNAde_mRNAde_300000[,2], 
        col=c("white", "red", "light blue", "blue", 
              "white", "red", "light blue", "blue", 
              "white", "red", "light blue", "blue"), 
        ylab="Transcripts per mRNA", notch=TRUE, outline=FALSE,
        at=c(1,2,3,4,5.5,6.5,7.5,8.5,10,11,12,13))
dev.off()

