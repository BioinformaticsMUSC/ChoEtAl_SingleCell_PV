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
library(Libra)
library(openxlsx)
})
source("utils/Utils.R")

#plan("multiprocess", workers = 10)
options(future.globals.maxSize = 6000 * 1024^2)

dir.create("output_relabel/DGE_MAST")

load("output_relabel/JenniferCho_SeuratObj_Final_NoDoublet_Relabel.RData")

tmp <- seuObject_slim_nodoub

tmp@meta.data$cell_type <- tmp@meta.data$Cell
tmp@meta.data$replicate <- gsub("_F","",tmp@meta.data$Genotype)
tmp@meta.data$label <- gsub("_F","",tmp@meta.data$Genotype)

dge <- run_de(tmp, de_family = 'singlecell', de_method = 'MAST',n_threads = 20)

save(dge, file = "output_relabel/DGE_MAST/JennCho_DGE.RData")

# Save as multiple table facet

list <- split(dge,dge$cell_type)
wb <- createWorkbook()
Map(function(data, nameofsheet){     
    addWorksheet(wb, nameofsheet)
    writeData(wb, nameofsheet, data)
}, list, names(list))
saveWorkbook(wb, file = "output_relabel/DGE_MAST/JennCho_Dge_FullStat.xlsx", overwrite = TRUE)


# Filter Sign
dge_sign <- dge %>% filter(p_val_adj < 0.05 & abs(avg_logFC) > 0.2)


openxlsx::write.xlsx(dge_sign, 
                     file = "output_relabel/DGE_MAST/DGE_Sign.xlsx", 
                     colNames = TRUE,
                     rowNames = FALSE, 
                     borders = "columns",
                     sheetName="Stats",
                     overwrite=T)
write.table(dge_sign,"output_relabel/DGE_MAST/DGE_Sign.txt",sep="\t",quote=F)

# Convert with human ID
# Human ID
human = biomaRt::useMart(host="https://dec2021.archive.ensembl.org", "ENSEMBL_MART_ENSEMBL", dataset="hsapiens_gene_ensembl")
mouse = biomaRt::useMart(host="https://dec2021.archive.ensembl.org", "ENSEMBL_MART_ENSEMBL", dataset="mmusculus_gene_ensembl")

MGI = biomaRt::getLDS(attributes = c("mgi_symbol"), filters = "mgi_symbol", values = dge_sign$gene,mart = mouse, attributesL = c("hgnc_symbol"), martL = human, uniqueRows=T)

comb <- merge(dge_sign,MGI,by.x="gene",by.y="MGI.symbol",all=F)

tmp <- comb %>%
       dplyr::select(HGNC.symbol,cell_type) %>%
       dplyr::rename(Gene = HGNC.symbol, Cell = cell_type) %>%
       arrange(Cell)

write.table(tmp,"output_relabel/DGE_MAST/DGE_Sign_ForEnrich_HumanID.txt",sep="\t",quote=F,row.names=F)

tmp <- dge_sign %>%
       dplyr::select(gene,cell_type) %>%
       dplyr::rename(Gene = gene, Cell = cell_type) %>%
       arrange(Cell)

write.table(tmp,"output_relabel/DGE_MAST/DGE_Sign_ForEnrich_MouseID.txt",sep="\t",quote=F,row.names=F)

# Run enrichment for disease
system("sh 06_Enrichments.sh")

# END