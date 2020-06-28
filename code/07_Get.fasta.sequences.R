  ##################################
  # 07_Get_variants_and_check_overlaping_with_shape_motif.R
  #################################
  
library(dplyr)
library(GenomicRanges)
library(ggplot2)
library(ggseqlogo)
library(seqLogo)
library(Biostrings)
library(gridExtra)

  # Download annotated variants from DisgeNet. Apply "curated" and "disease related" filters. Available in the link below:
  # https://www.disgenet.org/browser/0/0/5/0/0/25/source__CURATED/type__0__disease/  

  #Read Variants.bed and generate gr object:
  #variants <- read.table(paste("Variants/", "variants.bed", sep =""), header = FALSE)
  variants <- read.table(paste("Variants/", "variants.all.bed", sep =""), header = FALSE)
  variants$V5 <- variants$V4
  variants$V4 <- noquote(rep(".", nrow(variants)))
  colnames(variants) <- c("chrom", "start", "end", "strand", "id")
  variants <- variants[complete.cases(variants), ]
  variants.gr <- makeGRangesFromDataFrame(variants)
  
  #Check overlaping between shape motifs and variants
  list.train.chip <- list.files("gibbs_input/", pattern="positives_training_shape.dat")
  samples <- unlist(lapply(list.train.chip, function(x) {
    stringr::str_split_fixed(x, "[_]",2)[1,1]}))
  
  overlap.all <- NULL
  for(j in 1:length(samples)){
    #Check if sample files exist
    files <-list.files("motif_detector/", pattern=samples[j])
    if (length(files) == 0){
      next
    }
    else{
      #Read summary and filter good motifs
      info = file.info(paste("motif_detector/", samples[j], ".p_val_len_summary", sep=""))
      if (info$size == 0){
        next
      }
      else{
        results <- read.table(paste("motif_detector/", samples[j], ".p_val_len_summary", sep=""), header=FALSE)
        colnames(results) <- c("hyp.pval_train", "hyp.pval_val", "F1/3_train", "F1/3_val", "motif_length", "SD", "TPR_train", "FPR_train", "TPR_val", "FPR_val")
        good <- results[results$FPR_val < 0.2,]
        if (nrow(good)==0){
          next
        }
        else{
          overlap <- NULL
          #Check overlaping
          for(i in 1:nrow(good)){
            motif <- read.table(paste("motif_detector/", samples[j], "_", i, ".bed", sep=""), header=FALSE)
            colnames(motif) <- c("chrom", "start", "end", "strand")
            motif.gr <- makeGRangesFromDataFrame(motif)
            ov2 <- findOverlaps(motif.gr, variants.gr)
            motif.over <- motif[queryHits(ov2),]
            motif.over$var <- variants[subjectHits(ov2),]$id
            motif.over$motif <- noquote(rep(i, nrow(motif.over)))
            motif.over$chip <- noquote(rep(samples[j], nrow(motif.over)))
            overlap <- rbind(overlap, motif.over)
          }
          overlap.all <- rbind(overlap.all, overlap)
        }
  
        if(nrow(overlap>0)){
          #Save overlaping motifs
          write.table(overlap, paste("overlap/", samples[j], "_over.bed", sep=""), row.names=FALSE, col.names=TRUE)
          for(i in 1:nrow(overlap)) {
            
            #Generate shape plot
            shape <- read.table(paste("motif_detector/", samples[j], "_", overlap$motif[i], ".shape", sep=""), header=FALSE)
            error <- NULL
            mean <- NULL
            for(k in 1:ncol(shape)){
              error <- c(error, sd(shape[,k]))
              mean <- c(mean, mean(shape[,k]))
            }
            pdf(file=paste("overlap/", samples[j], "_", overlap$motif[i], "_shape.pdf", sep=""))  
            data <- data.frame(index=seq(1,length(mean),by=1), shape=mean, error=error)
            print(ggplot(data, aes(x=index, y=shape)) +
                    geom_line(color="red") +
                    geom_ribbon(aes(ymax=mean+error, ymin=mean-error), alpha = 0.2, fill="orange") +
                    labs(title="TF", subtitle="motif", y="Values of ProT", x="Position") +
                    theme_classic())
            dev.off()
            
            #Generate seqlogo
            pdf(file=paste("overlap/", samples[j], "_", overlap$motif[i], "_logo.pdf", sep=""))  
            dat.in <- readLines(paste("fasta/",samples[j], "_", overlap$motif[i], ".fa", sep=""))
            seqs <- dat.in[seq_along(dat.in) %%2 == 0]
            letter_counts <- consensusMatrix(seqs)
            probs <- prop.table(letter_counts[1:4,], 2)
            seqLogo(probs, ic.scale=FALSE)
            dev.off()
          }
            
  }
      }}}
      
      
  metadata <- read.table("ChIP-seq_data/filtered_datasets.tsv")
  metadata %>% dplyr::select("V17", "V8") -> metadata.crop
  colnames(metadata.crop) <- c("chip","TF")
  overlap.all <- merge(overlap.all, metadata.crop, by="chip")
  write.table(overlap.all, paste("overlap/", "summary_overlap.bed", sep=""), row.names=FALSE, col.names=TRUE)
  overlap.all <- read.table(paste("overlap/", "summary_overlap.bed", sep=""), header=TRUE)
  
