# Covariate plot
rm(list=ls())
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(RcisTarget))
suppressPackageStartupMessages(library(igraph))

# Create Output directories
dir.create("output_relabel/DGE_MAST/motif_enrichment/")
dge <- read.table("output_relabel/DGE_MAST/DGE_Sign_ForEnrich_MouseID.txt",header=T)

names <- as.character(unique(dge$Cell))

# Load the motid database
# Download the database from "https://resources.aertslab.org/cistarget/"
motifRankings <- importRankings("utils/mm9-tss-centered-10kb-10species.mc9nr.feather")
data(motifAnnotations_mgi)

# Loop across all the cool modules
input <- list()
geneLists <- list()
hits <- list()
enrich <- list()
for(i in 1:length(names))
	{
		input[[i]] <- dge %>% filter(Cell==names[[i]])
		geneLists[[i]] <- as.character(input[[i]]$Gene)
		enrich[[i]] <- cisTarget(geneLists[[i]], motifRankings,motifAnnot=motifAnnotations_mgi,nCores = 12) %>% addLogo()
		write.table(enrich[[i]],paste0("output_relabel/DGE_MAST/motif_enrichment/",names[[i]],"_Motif_Enrichment.txt"), sep="\t",quote=F)
	}
