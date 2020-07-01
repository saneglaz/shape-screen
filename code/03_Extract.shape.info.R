  ##################################
  # 03_Extract.shape.info.R
  #################################
  # Download the DNA shape (ProT feature) from rohslab.usc.edu: 
  download.file("ftp://rohslab.usc.edu/hg19/hg19.ProT.wig.bw", "hg19.ProT.wig.bw")
  # The screening can be done with any of the others DNA shape features. In this case, only ProT data was employed.
  # extract shape data using bwtool
  list.bed <- list.files("peaklist/", pattern=".bed")
  bwtool <- function(list.bed){
    name <- gsub(".bed", "_profile.dat", list.bed)
    temp <- read.table(paste(getwd(), "/peaklist/", list.bed, sep =""), header = FALSE)
    temp$V2 <- temp$V2-1
    write.table(temp, file = paste(getwd(), "/temp.bed", sep=""), sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)
      #CHANGE FOR YOUR OWN PATH TO BWTOOL
    x <- paste("/home/frigopie/anaconda3/envs/myproject/bin/bwtool extract bed ", 
               getwd(), "/temp.bed ", 
               getwd(), "/genome/", "hg19.ProT.wig.bw ",
               getwd(), "/Shape-profile/", name, sep="")
    system(x)
  }
  lapply(list.bed, bwtool)
  
