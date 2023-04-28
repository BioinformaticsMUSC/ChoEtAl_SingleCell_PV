suppressPackageStartupMessages({
library(tidyverse)
library(SingleCellExperiment)
library(scater)
library(BiocParallel)
library(ggpubr)
library(speckle)
library(magrittr)
library(broom)
library(muscat)
library(Seurat)
library(data.table)
library(BiocSingular)
library(SCopeLoomR)
})
source("utils/Utils.R")

dir.create("output_SCENIC")

load("output_relabel/JenniferCho_SeuratObj_Final_NoDoublet_Slimmed_Relabel.RData")

wt <- subset(x = seuObject_slim_nodoub_slim, subset = Genotype == "WT_F") 
het <- subset(x = seuObject_slim_nodoub_slim, subset = Genotype == "HET_F")

exprMat_wt <- wt@assays$RNA@data
cellInfo_wt <- wt@meta.data


exprMat_het <- het@assays$RNA@data
cellInfo_het <- het@meta.data

# Loom seurat
add_cell_annotation <- function(loom, cellAnnotation)
{
  cellAnnotation <- data.frame(cellAnnotation)
  if(any(c("nGene", "nUMI") %in% colnames(cellAnnotation)))
  {
    warning("Columns 'nGene' and 'nUMI' will not be added as annotations to the loom file.")
    cellAnnotation <- cellAnnotation[,colnames(cellAnnotation) != "nGene", drop=FALSE]
    cellAnnotation <- cellAnnotation[,colnames(cellAnnotation) != "nUMI", drop=FALSE]
  }
  
  if(ncol(cellAnnotation)<=0) stop("The cell annotation contains no columns")
  if(!all(get_cell_ids(loom) %in% rownames(cellAnnotation))) stop("Cell IDs are missing in the annotation")
  
  cellAnnotation <- cellAnnotation[get_cell_ids(loom),,drop=FALSE]
  # Add annotation
  for(cn in colnames(cellAnnotation))
  {
    add_col_attr(loom=loom, key=cn, value=cellAnnotation[,cn])
  }
  
  invisible(loom)
}

loom_het <- build_loom("output_SCENIC/sce_het.loom", dgem=exprMat_het)
loom_het <- add_cell_annotation(loom_het, cellInfo_het)
close_loom(loom_het)


loom_wt <- build_loom("output_SCENIC/sce_wt.loom", dgem=exprMat_wt)
loom_wt <- add_cell_annotation(loom_wt, cellInfo_wt)
close_loom(loom_wt)







