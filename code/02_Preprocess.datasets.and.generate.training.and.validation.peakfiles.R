  ##################################
  # 02_Preprocess.datasets.and.generate.training.and.validation.peakfiles.R
  #################################
  
library(diffloop)


  #Load ChIP-seq datasets
  list.datasets <- list.files("ChIP-seq_data/datasets/")
  create_shape.dat <- function(dataset){
    #load dataset
    coldata <- c("chrom", "start", "end", "name", "score", "strand", "signalValue", "pvalue", "qvalue", "peak")
    peaks <- read.table(paste("ChIP-seq_data/datasets/", dataset, sep =""), header = FALSE)
    colnames(peaks) <- coldata
    peaks$strand <- noquote(gsub(".", "+", peaks$strand))
    name <- stringr::str_split_fixed(dataset, "[.]",2)[1,1]
    
    #Filter HOT regions
    hot <- diffloop::bedToGRanges('/home/frigopie/BIOINFO_FORMACION/my_project/HOT/maphot_hs_selection_reg_k5_simP05_hot.bed')
    peaks.gr <- makeGRangesFromDataFrame(peaks, keep.extra.columns = TRUE)
    ov <- findOverlaps(peaks.gr, hot)
    peaks.gr.filtered <- peaks.gr[-subjectHits(ov),]
    
    #Create positive and negative subsets
    peaks.df.filtered <- data.frame(peaks.gr.filtered)
    peaks.df.filtered$length <- peaks.df.filtered$end - peaks.df.filtered$start
    peaks.df.filtered$center <- peaks.df.filtered$start + ((peaks.df.filtered$length)/2)
        #positive peaks are located +/-50bp around the peak summits
    peaks.positive <- peaks.df.filtered
    peaks.positive$start <- peaks.df.filtered$center - 50
    peaks.positive$end <- peaks.df.filtered$center + 50
    peaks.positive.gr <- makeGRangesFromDataFrame(peaks.positive, keep.extra.columns = TRUE)
         # each negative peak is located 100 bp upstream its corresponding positive peak
    peaks.negative <- peaks.df.filtered
    peaks.negative$start <- peaks.df.filtered$center - 250
    peaks.negative$end <- peaks.df.filtered$center - 150
    peaks.negative.gr <- makeGRangesFromDataFrame(peaks.negative, keep.extra.columns = TRUE)
    
    #Remove overlaping positive/negative peaks
    ov2 <- findOverlaps(peaks.positive.gr, peaks.negative.gr)
    peaks.positive.gr.over <- peaks.positive.gr[-queryHits(ov2),]
    peaks.negative.gr.over <- peaks.negative.gr[-subjectHits(ov2),]
    
    #Filter Top 1000 by signalValue
    if(length(peaks.positive.gr.over) >1000) {
      peaks.positive.gr.over.filtered <-sort(peaks.positive.gr.over, by= ~signalValue, decreasing=TRUE)[1:1000]
      peaks.negative.gr.over.filtered <-sort(peaks.negative.gr.over, by= ~signalValue, decreasing=TRUE)[1:1000]
    
      peaks.positive.df.over.filtered <- data.frame(seqnames=seqnames(peaks.positive.gr.over.filtered),
                                                  starts=start(peaks.positive.gr.over.filtered),
                                                  ends=end(peaks.positive.gr.over.filtered))
    
      peaks.negative.df.over.filtered <- data.frame(seqnames=seqnames(peaks.negative.gr.over.filtered),
                                                  starts=start(peaks.negative.gr.over.filtered),
                                                  ends=end(peaks.negative.gr.over.filtered))
      
      set.seed(42)
      rows <- sample(nrow(peaks.positive.df.over.filtered))
      peaks.positive.df.over.filtered.shuffled <- peaks.positive.df.over.filtered[rows,]
      peaks.positive.training <- peaks.positive.df.over.filtered.shuffled [1:500,]
      peaks.positive.validation <- peaks.positive.df.over.filtered.shuffled [501:1000,]
      
      rows <- sample(nrow(peaks.negative.df.over.filtered))
      peaks.negative.df.over.filtered.shuffled <- peaks.negative.df.over.filtered[rows,]
      peaks.negative.training <- peaks.negative.df.over.filtered.shuffled [1:500,]
      peaks.negative.validation <- peaks.negative.df.over.filtered.shuffled [501:1000,]
      
      
      #Save peak lists
      write.table(peaks.positive.training, file = paste("peaklist/", name, "_", "positives_training.bed", sep=""), sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)
      write.table(peaks.positive.validation, file = paste("peaklist/", name, "_", "positives_validation.bed", sep=""), sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)
      
      write.table(peaks.negative.training, file = paste("peaklist/", name, "_", "negatives_training.bed", sep=""), sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)
      write.table(peaks.negative.validation, file = paste("peaklist/", name, "_", "negatives_validation.bed", sep=""), sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)
     }
  }  
  lapply(list.datasets, create_shape.dat)
