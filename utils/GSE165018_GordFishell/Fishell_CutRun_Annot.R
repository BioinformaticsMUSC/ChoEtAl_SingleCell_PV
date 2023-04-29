suppressPackageStartupMessages({
library(tidyverse)
library(ggrepel)
library(DiffBind)
library(BiocParallel)
library(ggpubr)
library(magrittr)
library(broom)
library(data.table)
library(cowplot)
library(annotatr)
library(rtracklayer)
 library(regioneR)
})



annots = c('mm10_basicgenes', 'mm10_genes_intergenic')
annotations = build_annotations(genome = 'mm10', annotations = annots)


tmp <- toGRanges("Mef2c_pv_CutTag.bed")

annotated <- annotate_regions(
                                  regions = tmp,
                                  annotations = annotations,
                                  ignore.strand = TRUE,
                                  quiet = FALSE) %>%
                                  as.data.frame() %>%
                                dplyr::select(annot.symbol, annot.type) %>%
                                distinct() 

annotated <- na.omit(annotated)
colnames(annotated) <- c("Gene","Class")


annotated$Class <- gsub("mm10_genes_","",annotated$Class)

GeneSets <- split(annotated,annotated$Class)
save(GeneSets, file = "GeneSets_Mef2c_PVcutrun.RData")



tmp <- toGRanges("Mef2c_sst_CutTag.bed")

annotated <- annotate_regions(
                                  regions = tmp,
                                  annotations = annotations,
                                  ignore.strand = TRUE,
                                  quiet = FALSE) %>%
                                  as.data.frame() %>%
                                dplyr::select(annot.symbol, annot.type) %>%
                                distinct() 

annotated <- na.omit(annotated)
colnames(annotated) <- c("Gene","Class")


annotated$Class <- gsub("mm10_genes_","",annotated$Class)

GeneSets <- split(annotated,annotated$Class)
save(GeneSets, file = "GeneSets_Mef2c_SSTcutrun.RData")