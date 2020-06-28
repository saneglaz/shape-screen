  ##################################
  # 10_Peak.annotation.R
  #################################
  
library(ChIPseeker)
library(ggimage)
library(ggupset)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)


  #Annotating peaks to genes.
  summary.merge <- read.table(paste("results/", "summary_overlap.anno.bed", sep=""), header=TRUE)
  summary.merge.crop <- summary.merge[,3:6]
  summary.gr <- makeGRangesFromDataFrame(summary.merge.crop)
  filteredAnno <- annotatePeak(summary.gr, TxDb = TxDb.Hsapiens.UCSC.hg19.knownGene)
  
  
  pdf(file=paste("results/", "Peak_anno.pdf", sep=""))  
  #Displaying annotation distribution
  plotAnnoPie(filteredAnno)
  #Visualising overlaps between genomic regions.
  upsetplot(filteredAnno, vennpie=TRUE)
  #Visualize distribution of TF-binding loci relative to TSS
  plotDistToTSS(filteredAnno, title="Distribution of transcription factor-binding loci\nrelative to TSS")
  dev.off() 