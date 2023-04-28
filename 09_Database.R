# Covariate plot
rm(list=ls())
suppressPackageStartupMessages({
library(here)
library(tidyverse)
library(data.table)
})

dir.create("supp_tables")

# Load chip data
chip <- load(here("output_dge","DGE_MultiComp","CocaineCPPvsCocaineNPAS4_chipseqs_enrichment","ENRICHMENT_Npas4_Cell","Npas4_Cell_GeneSets.RData")) %>%
            get()
new_names <- names(chip)
chip <- map2(chip, new_names, ~setnames(.x, 'Class', .y))

# Load DGE data
load("output_dge/DGE_MultiComp/DGE_Mast_Analysis.RData")
load("output_dge/DGE_CocaineNPAS4_Specific/CocaineNPAS4_DEG.RData")


# CocaineCPPvsCocaineNPAS4
tmp <- DE_cocaine %>% 
        group_by(cell_type) %>%
        mutate(LOG = -log10(p_val_adj), ABS = abs(avg_logFC)) %>% 
        mutate(Threshold = if_else(p_val_adj < 0.05 & ABS > 0.2, "TRUE","FALSE")) %>%
        mutate(Direction = case_when(avg_logFC > 0.2 & p_val_adj < 0.05 ~ "UpReg", avg_logFC < -0.2 & p_val_adj < 0.05 ~ "DownReg")) %>%
        filter(Threshold == "TRUE") %>% 
        dplyr::select(gene,cell_type,avg_logFC,p_val_adj, Direction) %>%
        dplyr::rename(Gene = gene)

tmp <- list(dge = tmp)


l <- c(tmp,chip)
res <- Reduce(function(x, y) {
    merge(x, y, all=TRUE, by="Gene")
}, l)

database=res[!(is.na(res$cell_type)),]
database <- database[!(duplicated(database)),]

database1  <- database %>% arrange(cell_type)



# SalineCPPvsCocaineCPP
tmp <- DE_salcoc %>% 
        group_by(cell_type) %>%
        mutate(LOG = -log10(p_val_adj), ABS = abs(avg_logFC)) %>% 
        mutate(Threshold = if_else(p_val_adj < 0.05 & ABS > 0.2, "TRUE","FALSE")) %>%
        mutate(Direction = case_when(avg_logFC > 0.2 & p_val_adj < 0.05 ~ "UpReg", avg_logFC < -0.2 & p_val_adj < 0.05 ~ "DownReg")) %>%
        filter(Threshold == "TRUE") %>% 
        dplyr::select(gene,cell_type,avg_logFC,p_val_adj, Direction) %>%
        dplyr::rename(Gene = gene)

tmp <- list(dge = tmp)


l <- c(tmp,chip)
res <- Reduce(function(x, y) {
    merge(x, y, all=TRUE, by="Gene")
}, l)

database=res[!(is.na(res$cell_type)),]
database <- database[!(duplicated(database)),]

database2  <- database %>% arrange(cell_type)


# SalineCPPvsSalineNPAS4
tmp <- DE_saline %>% 
        group_by(cell_type) %>%
        mutate(LOG = -log10(p_val_adj), ABS = abs(avg_logFC)) %>% 
        mutate(Threshold = if_else(p_val_adj < 0.05 & ABS > 0.2, "TRUE","FALSE")) %>%
        mutate(Direction = case_when(avg_logFC > 0.2 & p_val_adj < 0.05 ~ "UpReg", avg_logFC < -0.2 & p_val_adj < 0.05 ~ "DownReg")) %>%
        filter(Threshold == "TRUE") %>% 
        dplyr::select(gene,cell_type,avg_logFC,p_val_adj, Direction) %>%
        dplyr::rename(Gene = gene)

tmp <- list(dge = tmp)


l <- c(tmp,chip)
res <- Reduce(function(x, y) {
    merge(x, y, all=TRUE, by="Gene")
}, l)

database=res[!(is.na(res$cell_type)),]
database <- database[!(duplicated(database)),]

database3  <- database %>% arrange(cell_type)

# SalineNPAS4vsCocaineNPAS4
tmp <- DE_salcoc_npas4 %>% 
        group_by(cell_type) %>%
        mutate(LOG = -log10(p_val_adj), ABS = abs(avg_logFC)) %>% 
        mutate(Threshold = if_else(p_val_adj < 0.05 & ABS > 0.2, "TRUE","FALSE")) %>%
        mutate(Direction = case_when(avg_logFC > 0.2 & p_val_adj < 0.05 ~ "UpReg", avg_logFC < -0.2 & p_val_adj < 0.05 ~ "DownReg")) %>%
        filter(Threshold == "TRUE") %>% 
        dplyr::select(gene,cell_type,avg_logFC,p_val_adj, Direction) %>%
        dplyr::rename(Gene = gene)
        
tmp <- list(dge = tmp)


l <- c(tmp,chip)
res <- Reduce(function(x, y) {
    merge(x, y, all=TRUE, by="Gene")
}, l)

database=res[!(is.na(res$cell_type)),]
database <- database[!(duplicated(database)),]

database4  <- database %>% arrange(cell_type)

# CocaineNPAS4 Specific
tmp <- CocaineNPAS4_Specific_Cleaned %>% 
        dplyr::select(gene,cell_type, Direction) %>%
        dplyr::rename(Gene = gene)
        
tmp <- list(dge = tmp)


l <- c(tmp,chip)
res <- Reduce(function(x, y) {
    merge(x, y, all=TRUE, by="Gene")
}, l)

database=res[!(is.na(res$cell_type)),]
database <- database[!(duplicated(database)),]

database5  <- database %>% arrange(cell_type)

# Supp table S5
list_of_datasets <- list("CocaineCPPvsCocaineNPAS4" = database1, "SalineCPPvsCocaineCPP" = database2, "SalineCPPvsSalineNPAS4" = database3, "SalineNPAS4vsCocaineNPAS4" = database4, "CocaineNPAS4_Specific" = database5)
openxlsx::write.xlsx(list_of_datasets, file = "supp_tables/Table_S5.xlsx",colNames = TRUE, borders = "columns",overwrite=TRUE)



