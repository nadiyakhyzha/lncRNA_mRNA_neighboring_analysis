# Title: 4_Neighboring_Pair_Subset_de
# Purpose: To retrieving the expression for all lncRNAs and mRNAs and then subsetting into those mRNAs with regulated or not regulated lncRNAs.
# Packages required: None
# Input files: "original_mRNAde_lncRNA_[distance]", "original_mRNAde_mRNA_[distance]", "lncRNA_IL1B_differentiallyexpressed", "mRNA_IL1B_differentiallyexpressed"
# Output files: "mRNAde_lncRNAde_[distance]", "mRNAde_lncRNA_[distance]", mRNAde_mRNAde_[distance], mRNAde_mRNA_[distance] 

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
# Saving the regulated_mRNA_lncRNA_list and unregulated_mRNA_lncRNA_list files as csv files

lapply(seq_along(regulated_mRNA_lncRNA_list),
       function(i)(write.csv(regulated_mRNA_lncRNA_list[[i]],
                             file =paste0('mRNAde_lncRNAde_',distance[i],'.csv'))))

lapply(seq_along(unregulated_mRNA_lncRNA_list),
       function(i)(write.csv(unregulated_mRNA_lncRNA_list[[i]],
                             file =paste0('mRNAde_lncRNA_',distance[i],'.csv'))))




#-----------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------
# Performing the same thing but for the mRNA-mRNA datasets

# Opening up the mRNA-lncRNA neighboring pairs for the selected genomic distances and combining all these dataframes in one list
mRNA_mRNA_list<-lapply(seq_along(distance), function(i)   assign(paste("mRNAde_mRNA_",distance[i], sep=""), read.csv(paste("Path/original_mRNAde_mRNA_", distance[i], ".csv", sep=""), stringsAsFactors=FALSE)[,2:3]))

# -----------------------------------------------------------------------------------------------------------------------
# Next searching through the list of differentially expressed mRNAs and retrieving their expression 
# fold changes in response to IL1B. All mRNAs that are not differentially expressed are given the value of "0"


mRNA_foldchange<-sapply(seq_along(mRNA_mRNA_list), function(i) sapply(1:nrow(mRNA_mRNA_list[[i]]), function(j) mRNA_IL1B$fold_change[mRNA_mRNA_list[[i]]$mRNA_Symbol[j]==mRNA_IL1B$mRNA]))
mRNA2_foldchange<-sapply(seq_along(mRNA_mRNA_list), function(i) sapply(1:nrow(mRNA_mRNA_list[[i]]), function(j) ifelse(sum(mRNA_mRNA_list[[i]]$mRNA2_Symbol[j]==mRNA_IL1B$mRNA)>=1, mRNA_IL1B$fold_change[mRNA_mRNA_list[[i]]$mRNA2_Symbol[j]==mRNA_IL1B$mRNA], 0 )))

mRNA_mRNA_list<-lapply(seq_along(mRNA_mRNA_list), function(i) cbind(mRNA_mRNA_list[[i]],mRNA_foldchange[[i]], mRNA2_foldchange[[i]]))

# -----------------------------------------------------------------------------------------------------------------------
# Separating the list into mRNAs with regulated mRNA2 or not regulated mRNA2
regulated_mRNA_mRNA_list<-mRNA_mRNA_list
unregulated_mRNA_mRNA_list<-mRNA_mRNA_list

# All of the mRNA-mRNA pairs where the mRNA2 is differentially expressed will be kept in the regulated_dataframe_list
regulated_mRNA_mRNA_list<-lapply(1:30, function(i) regulated_mRNA_mRNA_list[[i]][((regulated_mRNA_mRNA_list[[i]]$`mRNA2_foldchange[[i]]`>=1)|(regulated_mRNA_mRNA_list[[i]]$`mRNA2_foldchange[[i]]`<=-1)),])

# All of the mRNA-mRNA pairs where the mRNA2 is not differentially expressed will be kept in the unregulated_dataframe_list
unregulated_mRNA_mRNA_list<-lapply(1:30, function(i) unregulated_mRNA_mRNA_list[[i]][!((unregulated_mRNA_mRNA_list[[i]]$`mRNA2_foldchange[[i]]`>=1)|(unregulated_mRNA_mRNA_list[[i]]$`mRNA2_foldchange[[i]]`<=-1)),])

# -----------------------------------------------------------------------------------------------------------------------
# Saving the filtered files as csv's

lapply(seq_along(regulated_mRNA_mRNA_list),
       function(i)(write.csv(regulated_mRNA_mRNA_list[[i]],
                             file =paste0('mRNAde_mRNAde_',distance[i],'.csv'))))

lapply(seq_along(unregulated_mRNA_mRNA_list),
       function(i)(write.csv(unregulated_mRNA_mRNA_list[[i]],
                             file =paste0('mRNAde_mRNA_',distance[i],'.csv'))))
