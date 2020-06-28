##################################
# 01_Get.samples.R
#################################
setwd("../shape-screen")
  
  # Read and fix labels
    metadata <- read.csv2("http://hgdownload.soe.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeAwgTfbsUniform/files.txt", header=FALSE, stringsAsFactors = FALSE)
  colnames(metadata) <- c("path", "lab", "composite", "dataType", "view", "cell", "treatment", "TF", "control", "dataVersion", "dccAccession", "controlId", "quality", "tableName", "type", "md5sum", "size")
  for(i in 1:ncol(metadata)){
    metadata[,i] <- stringr::str_split_fixed(metadata[,i], "=", 2)[,2]
  }
    
  #FIlter for good quality, K562 celltype and no treatment
  metadata.good <- metadata[metadata$quality=="good" & metadata$treatment=="None" 
                      & (metadata$cell=="K562" | metadata$cell=="GM12878"),]
  metadata.good$url <- paste("http://hgdownload.soe.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeAwgTfbsUniform/", metadata.good$tableName, ".narrowPeak.gz", sep="")
  
  #Discard HDACs
  nokeep <- grepl("HDAC", metadata.good$TF, ignore.case = FALSE)
  metadata.good<- metadata.good[!nokeep,]
  
  #Save metadata
  write.table(metadata.good, file = "ChIP-seq_data/filtered_datasets.tsv", 
              sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)
  
  #Save manifest and download ChIPseq datasets
  write.table(metadata.good$url, file = "ChIP-seq_data/manifest.tsv", sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)
  x <- paste("wget -i ", getwd(), "/ChIP-seq_data/manifest.tsv -P ", getwd(), "/ChIP-seq_data/datasets/", sep="")
  cat(x)
  system(x)