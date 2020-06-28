  ##################################
  # 06_Get.fasta.sequences.R
  #################################
  
  list.bed <- list.files("motif_detector/", pattern=".bed")
  
  getfasta <- function(list.bed){
    name <- gsub(".bed", ".fa", list.bed)
    temp <- read.table(paste(getwd(), "/motif_detector/", list.bed, sep =""), header = FALSE)
    temp$V5 <- rep(1, nrow(temp))
    temp$V6 <- temp$V4
    temp$V4 <- seq(from=1, to=nrow(temp), by=1) 
    
    write.table(temp, file = paste(getwd(), "/temp.bed", sep=""), sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)
 #Download and place hg19 genome in genome folder    
    genome <- paste(getwd(), '/genomes/hg19_Genome.fa', sep="")
    input <- paste(getwd(), "/temp.bed", sep="")
    output <- paste(getwd(), "/fasta/", name,  sep="")
    x <- paste("bedtools getfasta -fi ", genome, " -bed ", input, " -fo ", output,  sep="")
    system(x)
  }
  lapply(list.bed, getfasta)
