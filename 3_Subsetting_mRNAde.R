# Title: 3_Subsetting_mRNAde
# Purpose: Take output files from FindingNeighboringlncRNA/ FindingNeighboringmRNA functions and filter out so that only
# neighbors to differentially expressed mRNAs (mRNAde) are selected
# Packages: None
# Input files: "mRNA_microarray_results", "mRNA_lncRNA_[distance]" and "mRNA_mRNA_[distance]"
# Output files: "original_mRNAde_lncRNA_[distance]", "original_mRNAde_mRNA_[distance]"
# -----------------------------------------------------------------------------------------------------------------------
# Loading the mRNA-lncRNA neighboring pair csv files 
# Start by setting up a list of genomic distances for which to open up the files
distance<-as.integer(c(seq(5000, 100000, 5000),125000,150000,175000,200000,
                       225000, 250000,275000, 300000, 400000, 500000))

# Opening up the mRNA-lncRNA neighboring pairs for the selected genomic distances and combining all these dataframes in one list
mRNA_lncRNA_list<-lapply(seq_along(distance), function(i)   assign(paste("mRNA_lncRNA_",distance[i], sep=""), read.csv(paste("Path/mRNA_lncRNA_", distance[i], ".csv", sep=""), stringsAsFactors=FALSE)[,2:3]))

# -----------------------------------------------------------------------------------------------------------------------
# Opening up the differentially expressed list of mRNAs
mRNA_IL1B<- read.csv("Path/mRNA_microarray_results.csv", stringsAsFactors=FALSE)

# -----------------------------------------------------------------------------------------------------------------------
# Here the mRNA-lncRNA pairs will be sifted to find which pairs have a mRNA that is differentially expressed.
# If the mRNA is not differentially expressed, it is marked as NA and than removed using the complete.cases function

test_mRNA_lncRNA<-sapply(1:30, function(i) sapply(1:nrow(mRNA_lncRNA_list[[i]]), function(j) {ifelse((sum(mRNA_lncRNA_list[[i]]$mRNA_Symbol[j]==mRNA_IL1B$mRNA)==0), NA, mRNA_lncRNA_list[[i]]$mRNA_Symbol[j])}))

mRNA_lncRNA_list<-lapply(1:30, function(i) mRNA_lncRNA_list[[i]][complete.cases(test_mRNA_lncRNA[[i]]),])

# -----------------------------------------------------------------------------------------------------------------------
# Saving the mRNAde_lncRNA neighboring pairs as csv's

lapply(seq_along(mRNA_lncRNA_list),
       function(i)(write.csv(mRNA_lncRNA_list[[i]],
                             file =paste0("original_mRNAde_lncRNA_",distance[i],'.csv'))))



# -----------------------------------------------------------------------------------------------------------------------
# -----------------------------------------------------------------------------------------------------------------------
# -----------------------------------------------------------------------------------------------------------------------

# Repeating the same thing but for the mRNA-mRNA list

# Opening up the mRNA-mRNA neighboring pairs for the selected genomic distances and combining all these dataframes in one list

mRNA_mRNA_list<-lapply(seq_along(distance), function(i) assign(paste("mRNA_mRNA_",distance[i], sep=""), read.csv(paste("Path/mRNA_mRNA_", distance[i], ".csv", sep=""), stringsAsFactors=FALSE)[,2:3]))

# -----------------------------------------------------------------------------------------------------------------------
# Here the mRNA-lncRNA pairs will be sifted to find which pairs have a mRNA that is differentially expressed.
# If the mRNA is not differentially expressed, it is marked as NA and than removed using the complete.cases function

test_mRNA_mRNA<-sapply(1:30, function(i) sapply(1:nrow(mRNA_mRNA_list[[i]]), function(j) {ifelse((sum(mRNA_mRNA_list[[i]]$mRNA_Symbol[j]==mRNA_IL1B$mRNA)==0), NA, mRNA_mRNA_list[[i]]$mRNA_Symbol[j])}))

mRNA_mRNA_list<-lapply(1:30, function(i) mRNA_mRNA_list[[i]][complete.cases(test_mRNA_mRNA[[i]]),])

# -----------------------------------------------------------------------------------------------------------------------
# Saving the filtered files as csv's

lapply(seq_along(mRNA_mRNA_list),
       function(i)(write.csv(mRNA_mRNA_list[[i]],
                             file =paste0("original_mRNAde_mRNA_",distance[i],'.csv'))))
