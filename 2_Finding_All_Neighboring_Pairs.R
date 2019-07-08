# Title: 2_Finding_Neighboring_Pairs
# Purpose: To identify all neighboring lncRNAs within a specific distance of any mRNA. 
# Packages required: "GenomicRanges"
require(GenomicRanges)
# Input files: "lncRNA_microarray_results", "mRNA_microarray_results" 
# Output files: "mRNA_lncRNA_[distance]" and "mRNA_mRNA_[distance]"

# -----------------------------------------------------------------------------------------------------------------------------------------------
# Opening up all lncRNAs in the microarray file  
lncRNA <- read.csv("Path/Table_S1_lncRNA_microarray_results.csv", stringsAsFactors=FALSE)[,1:11]

# Getting rid off incomplete rows
lncRNA<- lncRNA[complete.cases(lncRNA),]

# Removing duplicated lncRNAs
lncRNA<-lncRNA[!duplicated(lncRNA$lncRNA),]

# Turning lncRNA list into GRanges
lncRNA_GRanges<-makeGRangesFromDataFrame(lncRNA, keep.extra.columns=TRUE)


# Opening up all mRNAs in the microarray file
mRNA<- read.csv("Path/Table_S2_mRNA_microarray_results.csv", stringsAsFactors=FALSE)[,1:12]

# Getting rid off incomplete rows
mRNA<- mRNA[complete.cases(mRNA),]

# Removing duplicated lncRNAs
mRNA<-mRNA[!duplicated(mRNA$mRNA),]

# Turning mRNA list into GRanges
mRNA_GRanges<-makeGRangesFromDataFrame(mRNA, keep.extra.columns = TRUE)

# -----------------------------------------------------------------------------------------------------------------------------------------------
# The Finding Neighboring lncRNA function. It is based off of the overlaps function in the GRanges package

FindingNeighboringlncRNA<-function(mRNA, mRNA_GRanges, lncRNA, lncRNA_GRanges, upstream, downstream, output){

# The function will take in as input both mRNA/lncRNA list as well as the corresponding mRNA/lncRNA GRanges.
# Upstream and downstream are to specify the distance upstream or downstream from the start of the mRNA
# to search for presence of a lncRNA.
# Note that output must be a csv file
  
  mRNA<-mRNA
  mRNA_GRanges<-mRNA_GRanges
  lncRNA<-lncRNA
  lncRNA_GRanges<-lncRNA_GRanges
  upstream<-as.numeric(upstream)
  downstream<-as.numeric(downstream)
  output<-output
  
# This function will search for lncRNAs within a specified distances away from the start site of
# the mRNA. The pomoters function from the GRanges package is used to specify this distance.
  
  mRNA_GRanges<-promoters(mRNA_GRanges, upstream, downstream)
  
# The findOverlaps function is used to see whether a lncRNA falls within the specified distance
# upstream& downstream of an mRNA. Note that +/- strand specificity is ignored.
  
  overlaps<-findOverlaps(mRNA_GRanges, lncRNA_GRanges, type="any", ignore.strand=TRUE)
  
# The output of the findOverlaps function is numbers corresponding to the row number of mRNA/lncRNA
# from the input. Next step is to convert these numbers to the corresponding mRNA/lncRNA names.
# Start by creating columns where to record the mRNA/lncRNA names
  
  neighbors<-as.data.frame(overlaps)
  lncRNA_Symbol<-vector(length=nrow(neighbors))
  mRNA_Symbol<-vector(length=nrow(neighbors))
  neighbors<-cbind(neighbors, mRNA_Symbol, lncRNA_Symbol)
  
# Now retrieving the corresponding mRNA/lncRNA names
  
  neighbors$mRNA_Symbol<-sapply(1:nrow(neighbors), function(i) as.character(mRNA$mRNA[neighbors$queryHits[i]]))
  neighbors$lncRNA_Symbol<-sapply(1:nrow(neighbors), function(i) as.character(lncRNA$lncRNA[neighbors$subjectHits[i]]))

# Removing all of the extra unnecessary columns
  neighbors<-neighbors[,-(c(1:2))]
  
# Next step is to save the output file
  
  write.csv(neighbors, output)
}

# -----------------------------------------------------------------------------------------------------------------------------------------------
# Next the FindingNeighboringlncRNA function will be run across specified genomic distances

# Creating a vector of genomic distances to input into the FindingNeighboringlncRNA function
distance<-c(
  "5000", "10000", "15000", "20000", "25000", "30000", "35000", "40000", "45000", "50000", "55000", "60000", 
  "65000", "70000","75000", "80000","85000", "90000", "95000", 
  "100000", "125000", "150000", "175000", "200000", "225000", "250000", "275000",
  "300000", "400000", "500000")

# Running the FindingNeighboringlncRNA function for the specified genomic distances

sapply(1:length(distance), function(i) FindingNeighboringlncRNA(mRNA, mRNA_GRanges, lncRNA, lncRNA_GRanges, distance[i], distance[i], 
                                                                paste("mRNA_lncRNA_", distance[i], ".csv", sep="")))

# ------------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------------
# Repeating the same procedure but now looking for mRNA-mRNA pairs

# ------------------------------------------------------------------------------------------------------
# Creating the FindingNeighboringmRNA function.
# This function follows the same principles as the FindingNeighboringlncRNA function
# However, this function will filter out mRNAs that match up with their own transcripts and with
# their own splice variants

FindingNeighboringmRNA<-function(mRNA, mRNA_GRanges, mRNA2_GRanges, upstream, downstream, output){
  mRNA<-mRNA
  mRNA_GRanges<-mRNA_GRanges
  mRNA2_GRanges<-mRNA2_GRanges
  upstream<-as.numeric(upstream)
  downstream<-as.numeric(downstream)
  output<-output
  
  mRNA_GRanges<-promoters(mRNA_GRanges, upstream, downstream)
  overlaps<-findOverlaps(mRNA_GRanges, mRNA2_GRanges, type="any", ignore.strand=TRUE)
  
  neighbors<-as.data.frame(overlaps)
  mRNA_Symbol<-vector(length=nrow(neighbors))
  mRNA2_Symbol<-vector(length=nrow(neighbors))
  mRNA_Name<-vector(length=nrow(neighbors))
  mRNA2_Name<-vector(length=nrow(neighbors))
  neighbors<-cbind(neighbors, mRNA_Symbol, mRNA2_Symbol, mRNA_Name, mRNA2_Name)
  
  neighbors$mRNA_Symbol<-sapply(1:nrow(neighbors), function(i) as.character(mRNA$mRNA[neighbors$queryHits[i]]))
  neighbors$mRNA2_Symbol<-sapply(1:nrow(neighbors), function(i)as.character(mRNA$mRNA[neighbors$subjectHits[i]]))

# This addition here will removing mRNA neighbors that match up with themselves. It does so
# by removing mRNA-mRNA neighbors that have identical mRNA ID
  neighbors<-neighbors[neighbors$mRNA_Symbol!=neighbors$mRNA2_Symbol,]
  
# This addition will rid-off mRNAs that match up with their own splice variants
# First retrieve the gene names for each of the mRNA-mRNA neighbors
  neighbors$mRNA_Name<-sapply(1:nrow(neighbors), function(i) as.character(mRNA$mRNA.Symbol[neighbors$queryHits[i]]))
  neighbors$mRNA2_Name<-sapply(1:nrow(neighbors), function(i) as.character(mRNA$mRNA.Symbol[neighbors$subjectHits[i]]))
  
# Remove mRNA-mRNA pairs where both transcripts correspond to the same gene name
  neighbors<-neighbors[neighbors$mRNA_Name!=neighbors$mRNA2_Name,]
  
# Removing all of the extra unnecessary columns
  neighbors<-neighbors[,-(c(1:2,5:6))]
  
# Saving the output file 
  write.csv(neighbors, output)
}

# -----------------------------------------------------------------------------------------------------------------------------------------------
# Next the FindingNeighboringmRNA function will be run across specified genomic distances
# The vector with genomic distances has been created in the previous section

# Running the FindingNeighboringlncRNA function for the specified genomic distances

sapply(1:length(distance), function(i) FindingNeighboringmRNA(mRNA, mRNA_GRanges, mRNA_GRanges, distance[i], distance[i], 
                                                                paste("mRNA_mRNA_", distance[i], ".csv", sep="")))
