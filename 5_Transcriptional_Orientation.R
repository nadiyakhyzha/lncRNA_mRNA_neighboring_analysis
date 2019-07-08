# Title: 5_Transcriptional_Orientation
# Purpose: To assign the type of transcriptional orientation between the mRNA-lncRNA/ mRNA-mRNA pairs
# Packages: "GenomicRanges"
require(GenomicRanges)
# Input files: "mRNAde_lncRNAde_[distance]", "mRNAde_lncRNA_[distance]", mRNAde_mRNAde_[distance], mRNAde_mRNA_[distance],
# "lncRNA_IL1B_differentiallyexpressed", "mRNA_IL1B_differentiallyexpressed"
# Output files: "mRNAde_lncRNA_transcriptional_classification_[distance]", "mRNAde_lncRNAde_transcriptional_classification_[distance]",
# "mRNAde_mRNA_transcriptional_classification_[distance]", "mRNAde_mRNA2_transcriptional_classification_[distance]",
# "mRNAdelncRNAde_TranscriptionalOrientation_[distance]", and "mRNAdemRNAde_TranscriptionalOrientation_[distance]"

# -----------------------------------------------------------------------------------------------------------------------
# Loading the mRNA-lncRNA neighboring pair csv files 
# Start by setting up a list of genomic distances for which to open up the files
distance<-as.integer(c(seq(5000, 100000, 5000),125000,150000,175000,200000,
                       225000, 250000,275000, 300000, 400000, 500000))

# Opening up the mRNA-lncRNA neighboring pairs for the selected genomic distances and combining all these dataframes in one list
mRNAde_lncRNA_list<-lapply(seq_along(distance), function(i)   assign(paste("mRNAde_lncRNA_",distance[i], sep=""), read.csv(paste("Path/mRNAde_lncRNA_",distance[i],".csv", sep=""), stringsAsFactors=FALSE)[,c(2:3)]))
mRNAde_lncRNAde_list<-lapply(seq_along(distance), function(i)   assign(paste("mRNAde_lncRNAde_",distance[i], sep=""), read.csv(paste("Path/mRNAde_lncRNAde_",distance[i],".csv", sep=""), stringsAsFactors=FALSE)[,c(2:3)]))

# -----------------------------------------------------------------------------------------------------------------------------------------------
# Now opening up the files with the coordinates of lncRNAs and mRNAs

lncRNA_coordinates <- read.csv("Path/lncRNA_IL1B_differentiallyexpressed.csv")
lncRNA_coordinates<-lncRNA_coordinates[complete.cases(lncRNA_coordinates),]
lncRNA_coordinates<-lncRNA_coordinates[!duplicated(lncRNA_coordinates$lncRNA),]

mRNA_coordinates <- read.csv("Path/mRNA_IL1B_differentiallyexpressed.csv")
mRNA_coordinates<-mRNA_coordinates[complete.cases(mRNA_coordinates),]
mRNA_coordinates<-mRNA_coordinates[!duplicated(mRNA_coordinates$mRNA),]

#-----------------------------------------------------------------------------------------------------------
# Here use the names of mRNA/lncRNA of interest to retrieve its coordinates from the
# lncRNA_coordinates/ mRNA_coordinates dataframe.
# The coordinates data is than combined

# Perform first for lncRNAs from mRNAde_lncRNA pairs
lncRNA_chr<-sapply(seq_along(distance), function(i) sapply(1:nrow(mRNAde_lncRNA_list[[i]]), function(j) as.character(lncRNA_coordinates[as.character(mRNAde_lncRNA_list[[i]]$lncRNA_Symbol[j])==as.character(lncRNA_coordinates$lncRNA),2])))
lncRNA_start<-sapply(seq_along(distance), function(i) sapply(1:nrow(mRNAde_lncRNA_list[[i]]), function(j) as.character(lncRNA_coordinates[as.character(mRNAde_lncRNA_list[[i]]$lncRNA_Symbol[j])==as.character(lncRNA_coordinates$lncRNA),3])))
lncRNA_end<-sapply(seq_along(distance), function(i) sapply(1:nrow(mRNAde_lncRNA_list[[i]]), function(j) as.character(lncRNA_coordinates[as.character(mRNAde_lncRNA_list[[i]]$lncRNA_Symbol[j])==as.character(lncRNA_coordinates$lncRNA),4])))
lncRNA_strand<-sapply(seq_along(distance), function(i) sapply(1:nrow(mRNAde_lncRNA_list[[i]]), function(j) as.character(lncRNA_coordinates[as.character(mRNAde_lncRNA_list[[i]]$lncRNA_Symbol[j])==as.character(lncRNA_coordinates$lncRNA),5])))

lncRNA_frommRNAdelncRNA_list<-lapply(seq_along(distance), function(i) as.data.frame(cbind(mRNAde_lncRNA_list[[i]]$lncRNA_Symbol, lncRNA_chr[[i]], lncRNA_start[[i]], lncRNA_end[[i]], lncRNA_strand[[i]]), stringsAsFactors=FALSE))

# Next perform for mRNAs from mRNAde_lncRNA pairs
mRNA_chr<-sapply(seq_along(distance), function(i) sapply(1:nrow(mRNAde_lncRNA_list[[i]]), function(j) as.character(mRNA_coordinates[as.character(mRNAde_lncRNA_list[[i]]$mRNA_Symbol[j])==as.character(mRNA_coordinates$mRNA),3])))
mRNA_start<-sapply(seq_along(distance), function(i) sapply(1:nrow(mRNAde_lncRNA_list[[i]]), function(j) as.character(mRNA_coordinates[as.character(mRNAde_lncRNA_list[[i]]$mRNA_Symbol[j])==as.character(mRNA_coordinates$mRNA),4])))
mRNA_end<-sapply(seq_along(distance), function(i) sapply(1:nrow(mRNAde_lncRNA_list[[i]]), function(j) as.character(mRNA_coordinates[as.character(mRNAde_lncRNA_list[[i]]$mRNA_Symbol[j])==as.character(mRNA_coordinates$mRNA),5])))
mRNA_strand<-sapply(seq_along(distance), function(i) sapply(1:nrow(mRNAde_lncRNA_list[[i]]), function(j) as.character(mRNA_coordinates[as.character(mRNAde_lncRNA_list[[i]]$mRNA_Symbol[j])==as.character(mRNA_coordinates$mRNA),6])))

mRNA_frommRNAdelncRNA_list<-lapply(seq_along(distance), function(i) as.data.frame(cbind(mRNAde_lncRNA_list[[i]]$mRNA_Symbol, mRNA_chr[[i]], mRNA_start[[i]], mRNA_end[[i]], mRNA_strand[[i]]), stringsAsFactors=FALSE))

# Next perform for lncRNAs from mRNAde_lncRNAde pairs
lncRNA_chr<-sapply(seq_along(distance), function(i) sapply(1:nrow(mRNAde_lncRNAde_list[[i]]), function(j) as.character(lncRNA_coordinates[as.character(mRNAde_lncRNAde_list[[i]]$lncRNA_Symbol[j])==as.character(lncRNA_coordinates$lncRNA),2])))
lncRNA_start<-sapply(seq_along(distance), function(i) sapply(1:nrow(mRNAde_lncRNAde_list[[i]]), function(j) as.character(lncRNA_coordinates[as.character(mRNAde_lncRNAde_list[[i]]$lncRNA_Symbol[j])==as.character(lncRNA_coordinates$lncRNA),3])))
lncRNA_end<-sapply(seq_along(distance), function(i) sapply(1:nrow(mRNAde_lncRNAde_list[[i]]), function(j) as.character(lncRNA_coordinates[as.character(mRNAde_lncRNAde_list[[i]]$lncRNA_Symbol[j])==as.character(lncRNA_coordinates$lncRNA),4])))
lncRNA_strand<-sapply(seq_along(distance), function(i) sapply(1:nrow(mRNAde_lncRNAde_list[[i]]), function(j) as.character(lncRNA_coordinates[as.character(mRNAde_lncRNAde_list[[i]]$lncRNA_Symbol[j])==as.character(lncRNA_coordinates$lncRNA),5])))

lncRNA_frommRNAdelncRNAde_list<-lapply(seq_along(distance), function(i) as.data.frame(cbind(mRNAde_lncRNAde_list[[i]]$lncRNA_Symbol, lncRNA_chr[[i]], lncRNA_start[[i]], lncRNA_end[[i]], lncRNA_strand[[i]]), stringsAsFactors=FALSE))

# Next perform for mRNAs from mRNAde_lncRNAde pairs
mRNA_chr<-sapply(seq_along(distance), function(i) sapply(1:nrow(mRNAde_lncRNAde_list[[i]]), function(j) as.character(mRNA_coordinates[as.character(mRNAde_lncRNAde_list[[i]]$mRNA_Symbol[j])==as.character(mRNA_coordinates$mRNA),3])))
mRNA_start<-sapply(seq_along(distance), function(i) sapply(1:nrow(mRNAde_lncRNAde_list[[i]]), function(j) as.numeric(mRNA_coordinates[as.character(mRNAde_lncRNAde_list[[i]]$mRNA_Symbol[j])==as.character(mRNA_coordinates$mRNA),4])))
mRNA_end<-sapply(seq_along(distance), function(i) sapply(1:nrow(mRNAde_lncRNAde_list[[i]]), function(j) as.numeric(mRNA_coordinates[as.character(mRNAde_lncRNAde_list[[i]]$mRNA_Symbol[j])==as.character(mRNA_coordinates$mRNA),5])))
mRNA_strand<-sapply(seq_along(distance), function(i) sapply(1:nrow(mRNAde_lncRNAde_list[[i]]), function(j) as.character(mRNA_coordinates[as.character(mRNAde_lncRNAde_list[[i]]$mRNA_Symbol[j])==as.character(mRNA_coordinates$mRNA),6])))

mRNA_frommRNAdelncRNAde_list<-lapply(seq_along(distance), function(i) as.data.frame(cbind(mRNAde_lncRNAde_list[[i]]$mRNA_Symbol, mRNA_chr[[i]], mRNA_start[[i]], mRNA_end[[i]], mRNA_strand[[i]]), stringsAsFactors=FALSE))

# -----------------------------------------------------------------------------------------------------------------------------------------------
# Afterwards the function "transcriptional classifier" is used on the mRNA/lncRNA coordinates
# to define them as either being forward codirectional, reverse codirectional, convergent, divergent,
# antisense or intersecting

transcriptional_classifier_lncRNA<-function(lncRNA, mRNA){
 
  # The inputs for the function is a list of lncRNAs with their coordinates and a list of
  # mRNAs with their coordinates
  
  # Start with first making sure that all of the column labels are appropriate and that
  # the correct classes of data are used 
  
  lncRNA<-lncRNA
  colnames(lncRNA)<-c("lncRNA", "lncRNA_chr", "lncRNA_start", "lncRNA_end", "lncRNA_strand")
  lncRNA$lncRNA_start<-as.numeric(lncRNA$lncRNA_start)
  lncRNA$lncRNA_end<-as.numeric(lncRNA$lncRNA_end)
  lncRNA$lncRNA_strand<-as.character(lncRNA$lncRNA_strand)
  
  mRNA<-mRNA
  colnames(mRNA)<-c("mRNA", "mRNA_chr", "mRNA_start", "mRNA_end", "mRNA_strand")
  mRNA$mRNA_start<-as.numeric(mRNA$mRNA_start)
  mRNA$mRNA_end<-as.numeric(mRNA$mRNA_end)
  mRNA$mRNA_strand<-as.character(mRNA$mRNA_strand)

  lncRNA[((lncRNA$lncRNA_start>mRNA$mRNA_start)&(lncRNA$lncRNA_end>mRNA$mRNA_end)&(lncRNA$lncRNA_start<mRNA$mRNA_end)),6]<-as.character("antisense")
  lncRNA[((lncRNA$lncRNA_start<mRNA$mRNA_start)&(lncRNA$lncRNA_end<mRNA$mRNA_end)&(lncRNA$lncRNA_end>mRNA$mRNA_start)),6]<-as.character("antisense")
  
  lncRNA[((lncRNA$lncRNA_start>mRNA$mRNA_start)&(lncRNA$lncRNA_end>mRNA$mRNA_end)&(lncRNA$lncRNA_start>mRNA$mRNA_end)),6]<-as.character("downstream")
  lncRNA[((lncRNA$lncRNA_start<mRNA$mRNA_start)&(lncRNA$lncRNA_end<mRNA$mRNA_end)&(lncRNA$lncRNA_end<mRNA$mRNA_start)),6]<-as.character("upstream")
  
  lncRNA[((lncRNA$lncRNA_start<=mRNA$mRNA_start)&(lncRNA$lncRNA_end>=mRNA$mRNA_end))|((lncRNA$lncRNA_start>=mRNA$mRNA_start)&(lncRNA$lncRNA_end<=mRNA$mRNA_end)),6]<-as.character("intersecting")
  lncRNA[((lncRNA$lncRNA_start==mRNA$mRNA_start)&(lncRNA$lncRNA_end==mRNA$mRNA_end)),6]<-as.character("same coordinates")
  
  
  lncRNA[(lncRNA$lncRNA_strand=="+")&(mRNA$mRNA_strand=="+")&(lncRNA$V6=="upstream"),7]<-as.character("forward codirectional")
  lncRNA[(lncRNA$lncRNA_strand=="-")&(mRNA$mRNA_strand=="-")&(lncRNA$V6=="upstream"),7]<-as.character("reverse codirectional")
  lncRNA[(lncRNA$lncRNA_strand=="+")&(mRNA$mRNA_strand=="-")&(lncRNA$V6=="upstream"),7]<-as.character("convergent")
  lncRNA[(lncRNA$lncRNA_strand=="-")&(mRNA$mRNA_strand=="+")&(lncRNA$V6=="upstream"),7]<-as.character("divergent")
  
  lncRNA[(lncRNA$lncRNA_strand=="+")&(mRNA$mRNA_strand=="+")&(lncRNA$V6=="downstream"),7]<-as.character("forward codirectional")
  lncRNA[(lncRNA$lncRNA_strand=="-")&(mRNA$mRNA_strand=="-")&(lncRNA$V6=="downstream"),7]<-as.character("reverse codirectional")
  lncRNA[(lncRNA$lncRNA_strand=="+")&(mRNA$mRNA_strand=="-")&(lncRNA$V6=="downstream"),7]<-as.character("divergent")
  lncRNA[(lncRNA$lncRNA_strand=="-")&(mRNA$mRNA_strand=="+")&(lncRNA$V6=="downstream"),7]<-as.character("convergent")

  lncRNA[(lncRNA$V6=="same coordinates"),7]<-as.character("same coordinates")
  lncRNA[(lncRNA$V6=="antisense"),7]<-as.character("antisense")
  lncRNA[(lncRNA$V6=="intersecting"),7]<-as.character("intersecting")
  
    return(lncRNA)
}

# -----------------------------------------------------------------------------------------------------------------------------------------------
# Running the transcription_classifier function for the mRNAde_lncRNA list and
# the mRNAde_lncRNAde list

mRNAde_lncRNA_list<-lapply(seq_along(distance), function(i) transcriptional_classifier_lncRNA(lncRNA_frommRNAdelncRNA_list[[i]], mRNA_frommRNAdelncRNA_list[[i]]))

mRNAde_lncRNAde_list<-lapply(seq_along(distance), function(i) transcriptional_classifier_lncRNA(lncRNA_frommRNAdelncRNAde_list[[i]], mRNA_frommRNAdelncRNAde_list[[i]]))

# -----------------------------------------------------------------------------------------------------------------------------------------------
# Saving the transcriptional orientation of mRNA-lncRNA pairs

lapply(seq_along(mRNAde_lncRNA_list),
       function(i)(write.csv(as.data.frame(mRNAde_lncRNA_list[[i]]),
                             file =paste0("mRNAde_lncRNA_transcriptional_classification",distance[i],'.csv'))))

lapply(seq_along(mRNAde_lncRNAde_list),
       function(i)(write.csv(as.data.frame(mRNAde_lncRNAde_list[[i]]),
                             file =paste0("mRNAde_lncRNAde_transcriptional_classification",distance[i],'.csv'))))



# -----------------------------------------------------------------------------------------------------------------------------------------------
# -----------------------------------------------------------------------------------------------------------------------------------------------
# -----------------------------------------------------------------------------------------------------------------------------------------------
# Repeating for the mRNA_mRNA datasets
# Loading the mRNA-lncRNA neighboring pair csv files 
# Start by setting up a list of genomic distances for which to open up the files
distance<-as.integer(c(seq(5000, 100000, 5000),125000,150000,175000,200000,
                       225000, 250000,275000, 300000, 400000, 500000))

# Opening up the mRNA-mRNA neighboring pairs for the selected genomic distances and combining all these dataframes in one list
mRNAde_mRNA_list<-lapply(seq_along(distance), function(i)   assign(paste("mRNAde_mRNA_",distance[i], "_up_downstream", sep=""), read.csv(paste("C:/Users/User/Desktop/Test/mRNAde_mRNA/mRNAde_mRNA_",distance[i],".csv", sep=""), stringsAsFactors=FALSE)[,c(2:3)]))
mRNAde_mRNAde_list<-lapply(seq_along(distance), function(i)   assign(paste("mRNAde_mRNAde_",distance[i], "_up_downstream", sep=""), read.csv(paste("C:/Users/User/Desktop/Test/mRNAde_mRNAde/mRNAde_mRNAde_",distance[i],".csv", sep=""), stringsAsFactors=FALSE)[,c(2:3)]))

# -----------------------------------------------------------------------------------------------------------------------------------------------
# Now opening up the file with the coordinates of mRNAs

mRNA_coordinates <- read.csv("Path/mRNA_IL1B_differentiallyexpressed.csv")
mRNA_coordinates<-mRNA_coordinates[complete.cases(mRNA_coordinates),]
mRNA_coordinates<-mRNA_coordinates[!duplicated(mRNA_coordinates$mRNA),]

#-----------------------------------------------------------------------------------------------------------
# Here use the names of mRNAs of interest to retrieve its coordinates from the
# mRNA_coordinates dataframe.
# The coordinates data is than combined

# Perform first for mRNAs from mRNAde_mRNA pairs
mRNA_chr<-sapply(seq_along(distance), function(i) sapply(1:nrow(mRNAde_mRNA_list[[i]]), function(j) as.character(mRNA_coordinates[as.character(mRNAde_mRNA_list[[i]]$mRNA_Symbol[j])==as.character(mRNA_coordinates$mRNA),3])))
mRNA_start<-sapply(seq_along(distance), function(i) sapply(1:nrow(mRNAde_mRNA_list[[i]]), function(j) as.character(mRNA_coordinates[as.character(mRNAde_mRNA_list[[i]]$mRNA_Symbol[j])==as.character(mRNA_coordinates$mRNA),4])))
mRNA_end<-sapply(seq_along(distance), function(i) sapply(1:nrow(mRNAde_mRNA_list[[i]]), function(j) as.character(mRNA_coordinates[as.character(mRNAde_mRNA_list[[i]]$mRNA_Symbol[j])==as.character(mRNA_coordinates$mRNA),5])))
mRNA_strand<-sapply(seq_along(distance), function(i) sapply(1:nrow(mRNAde_mRNA_list[[i]]), function(j) as.character(mRNA_coordinates[as.character(mRNAde_mRNA_list[[i]]$mRNA_Symbol[j])==as.character(mRNA_coordinates$mRNA),6])))

mRNA_frommRNAdemRNA_list<-lapply(seq_along(distance), function(i) as.data.frame(cbind(mRNAde_mRNA_list[[i]]$mRNA_Symbol, mRNA_chr[[i]], mRNA_start[[i]], mRNA_end[[i]], mRNA_strand[[i]]), stringsAsFactors=FALSE))

# Next perform for mRNA2s from mRNAde_mRNA pairs
mRNA2_chr<-sapply(seq_along(distance), function(i) sapply(1:nrow(mRNAde_mRNA_list[[i]]), function(j) as.character(mRNA_coordinates[as.character(mRNAde_mRNA_list[[i]]$mRNA2_Symbol[j])==as.character(mRNA_coordinates$mRNA),3])))
mRNA2_start<-sapply(seq_along(distance), function(i) sapply(1:nrow(mRNAde_mRNA_list[[i]]), function(j) as.character(mRNA_coordinates[as.character(mRNAde_mRNA_list[[i]]$mRNA2_Symbol[j])==as.character(mRNA_coordinates$mRNA),4])))
mRNA2_end<-sapply(seq_along(distance), function(i) sapply(1:nrow(mRNAde_mRNA_list[[i]]), function(j) as.character(mRNA_coordinates[as.character(mRNAde_mRNA_list[[i]]$mRNA2_Symbol[j])==as.character(mRNA_coordinates$mRNA),5])))
mRNA2_strand<-sapply(seq_along(distance), function(i) sapply(1:nrow(mRNAde_mRNA_list[[i]]), function(j) as.character(mRNA_coordinates[as.character(mRNAde_mRNA_list[[i]]$mRNA2_Symbol[j])==as.character(mRNA_coordinates$mRNA),6])))

mRNA2_frommRNAdemRNA_list<-lapply(seq_along(distance), function(i) as.data.frame(cbind(mRNAde_mRNA_list[[i]]$mRNA2_Symbol, mRNA2_chr[[i]], mRNA2_start[[i]], mRNA2_end[[i]], mRNA2_strand[[i]]), stringsAsFactors=FALSE))

# Next perform for mRNAs from mRNAde_mRNAde pairs
mRNA_chr<-sapply(seq_along(distance), function(i) sapply(1:nrow(mRNAde_mRNAde_list[[i]]), function(j) as.character(mRNA_coordinates[as.character(mRNAde_mRNAde_list[[i]]$mRNA_Symbol[j])==as.character(mRNA_coordinates$mRNA),3])))
mRNA_start<-sapply(seq_along(distance), function(i) sapply(1:nrow(mRNAde_mRNAde_list[[i]]), function(j) as.character(mRNA_coordinates[as.character(mRNAde_mRNAde_list[[i]]$mRNA_Symbol[j])==as.character(mRNA_coordinates$mRNA),4])))
mRNA_end<-sapply(seq_along(distance), function(i) sapply(1:nrow(mRNAde_mRNAde_list[[i]]), function(j) as.character(mRNA_coordinates[as.character(mRNAde_mRNAde_list[[i]]$mRNA_Symbol[j])==as.character(mRNA_coordinates$mRNA),5])))
mRNA_strand<-sapply(seq_along(distance), function(i) sapply(1:nrow(mRNAde_mRNAde_list[[i]]), function(j) as.character(mRNA_coordinates[as.character(mRNAde_mRNAde_list[[i]]$mRNA_Symbol[j])==as.character(mRNA_coordinates$mRNA),6])))

mRNA_frommRNAdemRNAde_list<-lapply(seq_along(distance), function(i) as.data.frame(cbind(mRNAde_mRNAde_list[[i]]$mRNA_Symbol, mRNA_chr[[i]], mRNA_start[[i]], mRNA_end[[i]], mRNA_strand[[i]]), stringsAsFactors=FALSE))

# Next perform for mRNA2s from mRNAde_mRNAde pairs
mRNA2_chr<-sapply(seq_along(distance), function(i) sapply(1:nrow(mRNAde_mRNAde_list[[i]]), function(j) as.character(mRNA_coordinates[as.character(mRNAde_mRNAde_list[[i]]$mRNA2_Symbol[j])==as.character(mRNA_coordinates$mRNA),3])))
mRNA2_start<-sapply(seq_along(distance), function(i) sapply(1:nrow(mRNAde_mRNAde_list[[i]]), function(j) as.numeric(mRNA_coordinates[as.character(mRNAde_mRNAde_list[[i]]$mRNA2_Symbol[j])==as.character(mRNA_coordinates$mRNA),4])))
mRNA2_end<-sapply(seq_along(distance), function(i) sapply(1:nrow(mRNAde_mRNAde_list[[i]]), function(j) as.numeric(mRNA_coordinates[as.character(mRNAde_mRNAde_list[[i]]$mRNA2_Symbol[j])==as.character(mRNA_coordinates$mRNA),5])))
mRNA2_strand<-sapply(seq_along(distance), function(i) sapply(1:nrow(mRNAde_mRNAde_list[[i]]), function(j) as.character(mRNA_coordinates[as.character(mRNAde_mRNAde_list[[i]]$mRNA2_Symbol[j])==as.character(mRNA_coordinates$mRNA),6])))

mRNA2_frommRNAdemRNAde_list<-lapply(seq_along(distance), function(i) as.data.frame(cbind(mRNAde_mRNAde_list[[i]]$mRNA2_Symbol, mRNA2_chr[[i]], mRNA2_start[[i]], mRNA2_end[[i]], mRNA2_strand[[i]]), stringsAsFactors=FALSE))


# -----------------------------------------------------------------------------------------------------------------------------------------------
#This part here test whether the coordinates of the mRNA and mRNA2 are in comparison to each other in genomic space

transcriptional_classifier_mRNA<-function(mRNA2, mRNA){
  mRNA2<-mRNA2
  mRNA<-mRNA
  

  colnames(mRNA2)<-c("mRNA2", "mRNA2_chr", "mRNA2_start", "mRNA2_end", "mRNA2_strand")
  mRNA2$mRNA2_start<-as.numeric(mRNA2$mRNA2_start)
  mRNA2$mRNA2_end<-as.numeric(mRNA2$mRNA2_end)
  mRNA2$mRNA2_strand<-as.character(mRNA2$mRNA2_strand)
  
  colnames(mRNA)<-c("mRNA", "mRNA_chr", "mRNA_start", "mRNA_end", "mRNA_strand")
  mRNA$mRNA_start<-as.numeric(mRNA$mRNA_start)
  mRNA$mRNA_end<-as.numeric(mRNA$mRNA_end)
  mRNA$mRNA_strand<-as.character(mRNA$mRNA_strand)
  
  mRNA2[((mRNA2$mRNA2_start>mRNA$mRNA_start)&(mRNA2$mRNA2_end>mRNA$mRNA_end)&(mRNA2$mRNA2_start<mRNA$mRNA_end)),6]<-as.character("antisense")
  mRNA2[((mRNA2$mRNA2_start<mRNA$mRNA_start)&(mRNA2$mRNA2_end<mRNA$mRNA_end)&(mRNA2$mRNA2_end>mRNA$mRNA_start)),6]<-as.character("antisense")
  
  mRNA2[((mRNA2$mRNA2_start>mRNA$mRNA_start)&(mRNA2$mRNA2_end>mRNA$mRNA_end)&(mRNA2$mRNA2_start>mRNA$mRNA_end)),6]<-as.character("downstream")
  mRNA2[((mRNA2$mRNA2_start<mRNA$mRNA_start)&(mRNA2$mRNA2_end<mRNA$mRNA_end)&(mRNA2$mRNA2_end<mRNA$mRNA_start)),6]<-as.character("upstream")
  mRNA2[((mRNA2$mRNA2_start<=mRNA$mRNA_start)&(mRNA2$mRNA2_end>=mRNA$mRNA_end))|((mRNA2$mRNA2_start>=mRNA$mRNA_start)&(mRNA2$mRNA2_end<=mRNA$mRNA_end)),6]<-as.character("intersecting")
  mRNA2[((mRNA2$mRNA2_start==mRNA$mRNA_start)&(mRNA2$mRNA2_end==mRNA$mRNA_end)),6]<-as.character("same coordinates")
  
  mRNA2[(mRNA2$mRNA2_strand=="+")&(mRNA$mRNA_strand=="+")&(mRNA2$V6=="upstream"),7]<-as.character("forward codirectional")
  mRNA2[(mRNA2$mRNA2_strand=="-")&(mRNA$mRNA_strand=="-")&(mRNA2$V6=="upstream"),7]<-as.character("reverse codirectional")
  mRNA2[(mRNA2$mRNA2_strand=="+")&(mRNA$mRNA_strand=="-")&(mRNA2$V6=="upstream"),7]<-as.character("covergent")
  mRNA2[(mRNA2$mRNA2_strand=="-")&(mRNA$mRNA_strand=="+")&(mRNA2$V6=="upstream"),7]<-as.character("divergent")
  
  mRNA2[(mRNA2$mRNA2_strand=="+")&(mRNA$mRNA_strand=="+")&(mRNA2$V6=="downstream"),7]<-as.character("forward codirectional")
  mRNA2[(mRNA2$mRNA2_strand=="-")&(mRNA$mRNA_strand=="-")&(mRNA2$V6=="downstream"),7]<-as.character("reverse codirectional")
  mRNA2[(mRNA2$mRNA2_strand=="+")&(mRNA$mRNA_strand=="-")&(mRNA2$V6=="downstream"),7]<-as.character("divergent")
  mRNA2[(mRNA2$mRNA2_strand=="-")&(mRNA$mRNA_strand=="+")&(mRNA2$V6=="downstream"),7]<-as.character("convergent")
  
  mRNA2[(mRNA2$V6=="same coordinates"),7]<-as.character("same coordinates")
  mRNA2[(mRNA2$V6=="antisense"),7]<-as.character("antisense")
  mRNA2[(mRNA2$V6=="intersecting"),7]<-as.character("intersecting")
  
  return(mRNA2)
}

# -----------------------------------------------------------------------------------------------------------------------------------------------
#Now running the transcriptional classifier function
mRNAde_mRNA_list<- lapply(seq_along(distance), function(i) transcriptional_classifier_mRNA(mRNA2_frommRNAdemRNA_list[[i]], mRNA_frommRNAdemRNA_list[[i]]))

mRNAde_mRNAde_list<- lapply(seq_along(distance), function(i) transcriptional_classifier_mRNA(mRNA2_frommRNAdemRNAde_list[[i]], mRNA_frommRNAdemRNAde_list[[i]]))


# -----------------------------------------------------------------------------------------------------------------------------------------------
#Saving the mRNA-mRNA transcriptional orientation classification

lapply(seq_along(mRNAde_mRNA_list),
       function(i)(write.csv(as.data.frame(mRNAde_mRNA_list[[i]]),
                             file =paste0("mRNAde_mRNA_transcriptional_classification_",distance[i],'.csv'))))

lapply(seq_along(mRNAde_mRNAde_list),
       function(i)(write.csv(as.data.frame(mRNAde_mRNAde_list[[i]]),
                             file =paste0("mRNAde_mRNAde_transcriptional_classification_",distance[i],'.csv'))))



# -----------------------------------------------------------------------------------------------------------------------------------------------
# -----------------------------------------------------------------------------------------------------------------------------------------------
# -----------------------------------------------------------------------------------------------------------------------------------------------
# Tabulating the number of different neighboring pair transcriptional orientations for
# mRNA-lncRNA and mRNA-mRNA pairs

# classifier_tabulation function will could the number of different types of transcriptional orientations
# there are

classifier_tabulator<-function(file){
  file<-file
  
  antisense<-vector(length=length(file))
  forward_codirectional<-vector(length=length(file))
  reverse_codirectional<-vector(length=length(file))
  convergent<-vector(length=length(file))
  divergent<-vector(length=length(file))
  intersecting<-vector(length=length(file))
  same_coordinates<-vector(length=length(file))
  
  output<-data.frame(forward_codirectional, reverse_codirectional, convergent, divergent,
                     intersecting, same_coordinates, antisense)
  
  rownames(output)<-names(file)
  
  for (i in 1:length(file)){
    output[i,1]<-(sum(file[[i]]$V7=="forward codirectional")/nrow(file[[i]]))*100
    output[i,2]<-(sum(file[[i]]$V7=="reverse codirectional")/nrow(file[[i]]))*100
    output[i,3]<-(sum(file[[i]]$V7=="convergent")/nrow(file[[i]]))*100
    output[i,4]<-(sum(file[[i]]$V7=="divergent")/nrow(file[[i]]))*100
    output[i,5]<-(sum(file[[i]]$V7=="intersecting")/nrow(file[[i]]))*100
    output[i,6]<-(sum(file[[i]]$V7=="same coordinates")/nrow(file[[i]]))*100
    output[i,7]<-(sum(file[[i]]$V7=="antisense")/nrow(file[[i]]))*100
      }
  return(output)
}

# Now running the classifier_tabulation function and saving the file

mRNAde_lncRNA_classification_tabulation<- classifier_tabulator(mRNAde_lncRNA_list)
mRNAde_lncRNAde_classification_tabulation<- classifier_tabulator(mRNAde_lncRNAde_list)
mRNAde_mRNA_classification_tabulation<- classifier_tabulator(mRNAde_mRNA_list)
mRNAde_mRNAde_classification_tabulation<- classifier_tabulator(mRNAde_mRNAde_list)

# -----------------------------------------------------------------------------------------------------------------------------------------------
# Note that there has been some missannotations for lncRNAs as some lncRNAs have shown up as
# having idential coordinates as their neighboring mRNAs. 
# Will remove all of the neighboring pairs where the mRNA and lncRNA have identical coordinates.
