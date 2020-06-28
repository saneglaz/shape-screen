  ##################################
  # 11_Functional_GO_Analysis.R 
  #################################

library(rGREAT)
library(GenomicRanges)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)

  variants <- read.table(paste("Variants/", "variants.all.bed", sep =""), header = FALSE)
  colnames(variants) <- c("chr", "start", "end", "var") 
  variants <- variants[complete.cases(variants),]
  variants <- variants[variants$chr!="chrMT",]
  variants.gr <- makeGRangesFromDataFrame(variants)
  bgAnno <- annotatePeak(variants.gr, TxDb = TxDb.Hsapiens.UCSC.hg19.knownGene)
  
  summary.merge <- read.table(paste("results/", "summary_overlap.anno.bed", sep=""), header=TRUE)
  summary.merge.crop <- summary.merge[,3:6]
  summary.gr <- makeGRangesFromDataFrame(summary.merge.crop)
  seqlevelsStyle(summary.gr) <- "UCSC"
  seqlevelsStyle(variants.gr) <- "UCSC"
  great_Job <- submitGreatJob(gr=summary.gr, bg=variants.gr, species = "hg19")
  availableCategories(great_Job)
  great_ResultTable = getEnrichmentTables(great_Job, category = "GO")
  save(great_ResultTable, file = "results/Great_Results.RData")
  
  