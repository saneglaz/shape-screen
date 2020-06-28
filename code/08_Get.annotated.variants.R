  ##################################
  # 09_Get.annotated.variants.R
  #################################

library(myvariant)
library(rsnps)
library(dplyr)
library(disgenet2r)


  overlap.all <- read.table(paste("overlap/", "summary_overlap.bed", sep=""), header=TRUE)
  variants <- read.table(paste("Variants/", "variants.all.bed", sep =""), header = FALSE)
  variants.anno <- getVariants(unique(overlap.all$var[which(!duplicated(overlap.all$var))]), fields="dbsnp")
  dput(variants.anno, file=paste("results/", "variants.anno.R", sep=""))
  variants.anno <- dget(file=paste("results/", "variants.anno.R", sep=""))
  
  # Integrate annotation info to the variants overlaping with structural motifs:
  variants.dis <- variant2disease(variant= overlap.all$var, database = "CURATED")
  variants.dis <- extract(variants.dis)
  variants.dis %>% dplyr::select("variantid", "disease_name", "diseaseid", "gene_symbol", "variant_consequence_type", "score") -> variants.dis.crop
  

  summary.merge <- merge(overlap.all, variants.dis.crop,by.x="var", by.y="variantid")
  summary.merge <- merge(summary.merge, variants.anno[,1:2], by.x="var", by.y="query")
  write.table(summary.merge, paste("results/", "summary_overlap.anno.bed", sep=""), row.names=FALSE, col.names=TRUE)
  
  