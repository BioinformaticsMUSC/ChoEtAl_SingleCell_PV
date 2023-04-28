suppressPackageStartupMessages({
library(tidyverse)
library(ggrepel)
library(emmeans)
library(SingleCellExperiment)
library(scater)
library(BiocParallel)
library(ggpubr)
library(speckle)
library(magrittr)
library(broom)
library(muscat)
library(Seurat)
library(clustree)
library(leiden)
library(data.table)
library(cowplot)
library(scDblFinder)
library(BiocSingular)
library(RColorBrewer)
library(MetBrewer)
library(scales)
library(Nebulosa)
library(presto)
})

load("output/JenniferCho_SeuratObj_Final_NoDoublet_Slimmed.RData")


# Differential Expression
markers<- presto::wilcoxauc(seuObject_slim_nodoub_slim, 'seurat_clusters', assay = 'data')

markers$group <- paste("Cluster",markers$group,sep="_")

markers_sign <- markers %>%
                filter(padj < 0.05, pct_in > 20, pct_out < 20) %>% 
                dplyr::rename(Gene = feature, Cluster = group)

openxlsx::write.xlsx(markers_sign, 
                        file = "output/Markers_Significant.xlsx", 
                        colNames = TRUE, 
                        borders = "columns",
                        overwrite=T)

write.table(markers_sign,"output/Markers_Significant.txt",sep="\t",quote=F,row.names=F)


df <- markers_sign %>% select(Gene,Cluster)
write.table(df,"output/Markers_Significant_DGE.txt",sep="\t",quote=F,row.names=F)

# Human ID
human = biomaRt::useMart(host="https://dec2021.archive.ensembl.org", "ENSEMBL_MART_ENSEMBL", dataset="hsapiens_gene_ensembl")
mouse = biomaRt::useMart(host="https://dec2021.archive.ensembl.org", "ENSEMBL_MART_ENSEMBL", dataset="mmusculus_gene_ensembl")

MGI = biomaRt::getLDS(attributes = c("mgi_symbol"), filters = "mgi_symbol", values = df$Gene,mart = mouse, attributesL = c("hgnc_symbol"), martL = human, uniqueRows=T)

comb <- merge(df,MGI,by.x="Gene",by.y="MGI.symbol",all=F)

tmp <- comb %>%
       dplyr::select(HGNC.symbol,Cluster) %>%
       dplyr::rename(Gene = HGNC.symbol) %>%
       arrange(Cluster)

write.table(tmp,"output/Markers_Significant_DGE_HumanID.txt",sep="\t",quote=F,row.names=F)

# Run Enrichment
system("R CMD BATCH --vanilla Fisher_enrich_Markers_AllenPFC.R")
system("R CMD BATCH --vanilla Fisher_enrich_Markers_AllenPFC_Human.R")


sce <- as.SingleCellExperiment(seuObject_slim_nodoub_slim)
saveRDS(sce,"output/JenniferCho_SeuratObj_Final_NoDoublet_Slimmed.rds")

# Nebulosa 
dir.create("output/Nebulosa")

gene_list <- c("Mef2c","Slc17a7","Cux2","Calb1","Ntng1","Il1rapl2","Foxp2","Fezf2","Elavl2","Pvalb","Sst","Vip","Cck","Adarb2","Mobp","Cxcl14","C1qc","Flt1","Pdgfra","C1qa","Fn1","Vtn")

doPlot <- function(sel_name) 
{
PLOT <- plot_density(sce, sel_name, reduction="UMAP")
print(PLOT)
ggsave(sprintf("output/Nebulosa/%s.pdf", sel_name))
}
lapply(unique(gene_list), doPlot)

