  ##################################
  # 10_DNA.shape.prediction.R
  #################################
  
library(MASS)
library(ggplot2)
library(dplyr)

  summary.merge <- read.table(paste("results/", "summary_overlap.anno.bed", sep=""), header=TRUE)
  summary.merge.crop <- summary.merge[,3:6]
  summary.merge.crop$start <- summary.merge.crop$start - 50
  summary.merge.crop$end <- summary.merge.crop$end + 50
  
  # Get reference sequences for selected regions with structural motifs: 
  write.table(summary.merge.crop, file = paste(getwd(), "/temp.bed", sep=""), sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)
  genome <- '/home/frigopie/BIOINFO_FORMACION/Genomes/hg19_Genome.fa'
  input <- paste(getwd(), "/temp.bed", sep="")
  output <- paste(getwd(), "/overlap/", "overlap.fa",  sep="")
  x <- paste("bedtools getfasta -fi ", genome, " -bed ", input, " -fo ", output,  sep="")
  system(x)
  dat.in <- readLines(paste("overlap/", "overlap.fa", sep=""))
  seqs <- dat.in[seq_along(dat.in) %%2 == 0]
  headers <- dat.in[seq_along(dat.in) %%2 != 0]
  
  # Insert variants to generate alternative sequences: 
  seqs_alt <- seqs
  for (i in 1:length(seqs)){
    if(grepl("ins", summary.merge$X_id[i])==TRUE){
      letter <- unlist(strsplit(summary.merge$X_id[i], "ins"))[2]
      pos_rs <- unlist(strsplit(summary.merge$X_id[i], "ins"))[1]
      pos_rs <- unlist(strsplit(pos_rs, "[.]"))[2]
      pos_rs_start <- unlist(strsplit(pos_rs, "_"))[1]
      pos_rs_end <- unlist(strsplit(pos_rs, "_"))[2]
      pos <- as.numeric(pos_rs_start) - summary.merge.crop$start[i] +1
      mut <- c(substr(seqs_alt[i], 1, pos), letter, substr(seqs_alt[i], pos+1, nchar(seqs_alt[i])))
      seqs_alt[i] <- paste(mut, collapse="", sep="") 
    }
    else{
      letter <- substr(summary.merge$X_id[i], nchar(summary.merge$X_id[i]), nchar(summary.merge$X_id[i]))
      pos_rs <- unlist(strsplit(summary.merge$X_id[i], "[.]"))[2]
      pos_relative <- substr(pos_rs, 1, nchar(pos_rs)-3)
      pos <- as.numeric(pos_relative) - summary.merge.crop$start[i]
      substr(seqs_alt[i], pos, pos) <- letter 
    }}
    
  
  fasta_alt <- c(rbind(headers, seqs_alt))
  fasta_ref <- dat.in
  write.table(fasta_alt, paste("overlap/", "fasta_ref.fa", sep=""), sep = "\t", row.names = FALSE, col.names = FALSE, quote=FALSE)
  write.table(fasta_alt, paste("overlap/", "fasta_alt.fa", sep=""), sep = "\t", row.names = FALSE, col.names = FALSE, quote=FALSE)
  
  
  # Predict DNA shape for reference and alternative sequences:
  library(DNAshapeR)
  index <- seq(from=1, to=length(seqs)*2, by=2)
  shape.ref <- list()
  shape.alt <- list()
  for (i in 1:length(seqs)){
    fasta.ref.temp <- fasta_ref[index[i]:(index[i]+1)]
    write.table(fasta.ref.temp, paste("overlap/", "fasta.ref.temp.fa", sep=""), sep = "\t", row.names = FALSE, col.names = FALSE, quote=FALSE)
    pred <- getShape("overlap/fasta.ref.temp.fa", shapeType="ProT", parse=TRUE)
    shape.ref[[i]] <- pred$ProT[1,]
    fasta.alt.temp <- fasta_alt[index[i]:(index[i]+1)]
    write.table(fasta.alt.temp, paste("overlap/", "fasta.alt.temp.fa", sep=""), sep = "\t", row.names = FALSE, col.names = FALSE, quote=FALSE)
    pred <- getShape("overlap/fasta.alt.temp.fa", shapeType="ProT", parse=TRUE)
    shape.alt[[i]] <- pred$ProT[1,]
  }
  
  # Predict DNA shape for all regions with structural motifs in order to represent the standard deviation:
  system.time(
  shape.matrix <- lapply(seq(from=1, to=nrow(summary.merge)), function(k) {
      (getShape(paste("fasta/", summary.merge$chip[k], "_", summary.merge$motif[k], ".extended.fa", sep=""), shapeType="ProT", parse=TRUE))$ProT %>% 
      write.matrix(paste("overlap/", summary.merge$chip[k], "_", summary.merge$motif[k], ".shape", sep=""))
  }))
  
  # Plot comparing DNA shapes for reference and alternative sequences
  # Imbalance measure is and indicator stimating how much the variant is affected the shape of the structural domain
  
  shape.plot <- function(i) {
    shape <- as.matrix(read.table(paste("overlap/", summary.merge$chip[i], "_", summary.merge$motif[i], ".shape", sep=""), as.is = TRUE))
    shape <- shape[,51:(ncol(shape)-50)]
    error <- NULL
    mean <- NULL
    for(k in 1:ncol(shape)){
      error <- c(error, sd(shape[,k]))
      mean <- c(mean, mean(shape[,k]))
    }
    ref.df <- data.frame(index=seq(from=1,to=(length(shape.ref[[i]])-100), by=1), shape=unlist(shape.ref[i])[51:(length(shape.ref[[i]])-50)], set=rep("ref",length(shape.ref[[i]])-100), mean=mean, error=error)
    alt.df <- data.frame(index=seq(from=1,to=(length(shape.ref[[i]])-100), by=1), shape=unlist(shape.alt[i])[51:(length(shape.ref[[i]])-50)], set=rep("alt",length(shape.ref[[i]])-100), mean=mean, error=error)
    shape.df <- rbind(ref.df, alt.df)
    imbalance <- sum(((abs(alt.df$shape - ref.df$shape))/abs(mean))/error)
    print(imbalance)
    pdf(file=paste("results/", summary.merge$chip[i], "_", summary.merge$motif[i], "_shape.pdf", sep=""))  
    print(ggplot(shape.df, aes(x=index, y=shape, color=set)) + geom_line() +
      labs(title=paste("Structural motif for", summary.merge$TF[i], "overlaping with", summary.merge$var[i], sep=" "), y="ProT shape", x="position", colour = paste("Imbalance score:", imbalance, sep=" ")) +
      geom_ribbon(aes(ymax=mean+error, ymin=mean-error), alpha = 0.2, fill="orange", colour = NA) +
      #coord_cartesian(ylim=c(30, 40)) + 
      scale_x_continuous(breaks=seq(0, length(shape.ref[[i]]), 2)) +
      theme_set(theme_classic()))
    dev.off()
    imbalance
  }
  imbalance <- NULL
  imbalance <- unlist(lapply(seq(1, nrow(summary.merge), by=1), shape.plot))
  
  summary.merge %>% rename("TF.x"="TF") %>% dplyr::select(!TF.y) %>% dplyr::mutate(imbalance) %>% dplyr::arrange(desc(imbalance)) -> summary.results
  write.table(summary.results, paste("results/", "summary_results.txt", sep=""), row.names=FALSE, col.names=TRUE)
