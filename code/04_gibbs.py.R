  ##################################
  # 04_gibbs.py.R
  #################################
  
  # prepare input files for gibbs.py module
  list.bed <- list.files("peaklist/", pattern=".bed")
  list.shape <- list.files("Shape-profile/", pattern=".dat")
  gibbs_input <- function(i){
    name <- gsub("profile.dat", "shape.dat", list.shape[i])
    peaks <- read.table(paste(getwd(),"/peaklist/", list.bed[i], sep=""), header=FALSE)  
    shape <- read.table(paste(getwd(),"/Shape-profile/", list.shape[i], sep=""), header=FALSE)
    shape.info.df <- data.frame(seqnames=peaks$V1,
                               starts=peaks$V2,
                               ends=peaks$V3,
                               len=shape$V4,
                               shape=shape$V5)
    
    write.table(shape.info.df, file = paste(getwd(), "/gibbs_input/", name, sep=""), sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)
  }  
  
  lapply(seq(from=1, to=length(list.bed)), gibbs_input)
  
  # Write bash scripts for gibbs.py module
  list.dat <- list.files("gibbs_input/", pattern=".dat")
  text <- NULL
  for(i in 1:length(list.dat)){
    gibbs.py <- "python2 ../shape-motif/gibbs.py"
    input <- paste("../gibbs_input/", list.dat[i], sep="") 
    parameters <- "1000 10 5" 
    output <- "../gibbs_out"
    x <- paste(gibbs.py, input, parameters, output)
    text <- c(text, x)
  }
  bash <- "#!/bin/bash"
  writeLines(c("#!/bin/bash", text[1:90]), "gibbs1.txt")
  writeLines(c("#!/bin/bash", text[91:180]), "gibbs2.txt")
  writeLines(c("#!/bin/bash", text[181:270]), "gibbs3.txt")
  writeLines(c("#!/bin/bash", text[271:360]), "gibbs4.txt")
  
# gibbs scripts should be run in a powerful workstation, since this step is quite demanding