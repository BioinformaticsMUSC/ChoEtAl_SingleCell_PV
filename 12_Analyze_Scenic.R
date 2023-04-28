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

load("output_relabel/JenniferCho_SeuratObj_Final_NoDoublet_Relabel.RData")

wt <- subset(x = seuObject_slim_nodoub, subset = Genotype == "WT_F") 
het <- subset(x = seuObject_slim_nodoub, subset = Genotype == "HET_F")

regulome_auc_wt <- read.csv("output_SCENIC/AUC_wt.csv")
colnames(regulome_auc_wt) <- gsub("\\.","",colnames(regulome_auc_wt))

regulome_auc_het <- read.csv("output_SCENIC/AUC_het.csv")
colnames(regulome_auc_het) <- gsub("\\.","",colnames(regulome_auc_het))


#breaksList = seq(0, 1, by = 0.5)


# PVALB
tmp_wt <- regulome_auc_wt %>% 
                    mutate(Cell = wt@meta.data$Cell) %>%
                    pivot_longer(!Cell, names_to = "Gene", values_to = "AUC_wt") %>%
                    group_by(Cell,Gene) %>%
                    filter(AUC_wt == max(AUC_wt)) %>%
                    filter(Cell == "Inh_Pvalb", AUC_wt > 0) %>%
                    arrange(desc(AUC_wt))


tmp_het <- regulome_auc_het %>% 
                    mutate(Cell = het@meta.data$Cell) %>%
                    pivot_longer(!Cell, names_to = "Gene", values_to = "AUC_het") %>%
                    group_by(Cell,Gene) %>%
                    filter(AUC_het == max(AUC_het)) %>%
                    filter(Cell == "Inh_Pvalb", AUC_het > 0) %>%
                    arrange(desc(AUC_het))


tmp <- right_join(tmp_wt,tmp_het,by=c("Cell","Gene")) %>%
                  mutate(AUC_wt = replace_na(AUC_wt, 0), AUC_het = replace_na(AUC_het, 0)) %>%
                  arrange(desc(AUC_het)) %>%
                  as_tibble() %>%
                  top_n(n = 30) %>%
                  melt() %>%
                  distinct() %>%
                  select(Gene,variable,value) %>% 
                  pivot_wider(names_from = variable, values_from = value) %>%
                  top_n(n = 5, wt = AUC_het) %>%
                  column_to_rownames("Gene") %>%
                  as.matrix() 

pdf("output_SCENIC/Heatmap_Regulome_Pvalb.pdf", width = 3, height = 3)
pheatmap(tmp,scale="none",show_rownames = T, color = viridis::viridis(n = 30, alpha = 1, begin = 0, end = 1, option = "inferno"))
dev.off()

# SST
tmp_wt <- regulome_auc_wt %>% 
                    mutate(Cell = wt@meta.data$Cell) %>%
                    pivot_longer(!Cell, names_to = "Gene", values_to = "AUC_wt") %>%
                    group_by(Cell,Gene) %>%
                    filter(AUC_wt == max(AUC_wt)) %>%
                    filter(Cell == "Inh_Sst", AUC_wt > 0) %>%
                    arrange(desc(AUC_wt))

regulome_auc_het <- read.csv("output_SCENIC/AUC_het.csv")
colnames(regulome_auc_het) <- gsub("\\.","",colnames(regulome_auc_het))

tmp_het <- regulome_auc_het %>% 
                    mutate(Cell = het@meta.data$Cell) %>%
                    pivot_longer(!Cell, names_to = "Gene", values_to = "AUC_het") %>%
                    group_by(Cell,Gene) %>%
                    filter(AUC_het == max(AUC_het)) %>%
                    filter(Cell == "Inh_Sst", AUC_het > 0) %>%
                    arrange(desc(AUC_het))


tmp <- right_join(tmp_wt,tmp_het,by=c("Cell","Gene")) %>%
                  mutate(AUC_wt = replace_na(AUC_wt, 0), AUC_het = replace_na(AUC_het, 0)) %>%
                  arrange(desc(AUC_het)) %>%
                  as_tibble() %>%
                  top_n(n = 30) %>%
                  melt() %>%
                  distinct() %>%
                  select(Gene,variable,value) %>% 
                  pivot_wider(names_from = variable, values_from = value) %>%
                  top_n(n = 5, wt = AUC_het) %>%
                  column_to_rownames("Gene") %>%
                  as.matrix() 

pdf("output_SCENIC/Heatmap_Regulome_SST.pdf", width = 3, height = 3)
pheatmap::pheatmap(tmp,scale="none",show_rownames = T, color = viridis::viridis(n = 30, alpha = 1, begin = 0, end = 1, option = "inferno"))
dev.off()


# Vip
tmp_wt <- regulome_auc_wt %>% 
                    mutate(Cell = wt@meta.data$Cell) %>%
                    pivot_longer(!Cell, names_to = "Gene", values_to = "AUC_wt") %>%
                    group_by(Cell,Gene) %>%
                    filter(AUC_wt == max(AUC_wt)) %>%
                    filter(Cell == "Inh_Vip", AUC_wt > 0) %>%
                    arrange(desc(AUC_wt))

tmp_het <- regulome_auc_het %>% 
                    mutate(Cell = het@meta.data$Cell) %>%
                    pivot_longer(!Cell, names_to = "Gene", values_to = "AUC_het") %>%
                    group_by(Cell,Gene) %>%
                    filter(AUC_het == max(AUC_het)) %>%
                    filter(Cell == "Inh_Vip", AUC_het > 0) %>%
                    arrange(desc(AUC_het))


tmp <- right_join(tmp_wt,tmp_het,by=c("Cell","Gene")) %>%
                  mutate(AUC_wt = replace_na(AUC_wt, 0), AUC_het = replace_na(AUC_het, 0)) %>%
                  arrange(desc(AUC_het)) %>%
                  as_tibble() %>%
                  top_n(n = 30) %>%
                  melt() %>%
                  distinct() %>%
                  select(Gene,variable,value) %>% 
                  pivot_wider(names_from = variable, values_from = value) %>%
                  top_n(n = 5, wt = AUC_het) %>%
                  column_to_rownames("Gene") %>%
                  as.matrix() 

pdf("output_SCENIC/Heatmap_Regulome_VIP.pdf", width = 3, height = 3)
pheatmap::pheatmap(tmp,scale="none",show_rownames = T, color = viridis::viridis(n = 30, alpha = 1, begin = 0, end = 1, option = "inferno"))
dev.off()


# Lamp5
tmp_wt <- regulome_auc_wt %>% 
                    mutate(Cell = wt@meta.data$Cell) %>%
                    pivot_longer(!Cell, names_to = "Gene", values_to = "AUC_wt") %>%
                    group_by(Cell,Gene) %>%
                    filter(AUC_wt == max(AUC_wt)) %>%
                    filter(Cell == "Inh_Lamp5", AUC_wt > 0) %>%
                    arrange(desc(AUC_wt))

tmp_het <- regulome_auc_het %>% 
                    mutate(Cell = het@meta.data$Cell) %>%
                    pivot_longer(!Cell, names_to = "Gene", values_to = "AUC_het") %>%
                    group_by(Cell,Gene) %>%
                    filter(AUC_het == max(AUC_het)) %>%
                    filter(Cell == "Inh_Lamp5", AUC_het > 0) %>%
                    arrange(desc(AUC_het))


tmp <- right_join(tmp_wt,tmp_het,by=c("Cell","Gene")) %>%
                  mutate(AUC_wt = replace_na(AUC_wt, 0), AUC_het = replace_na(AUC_het, 0)) %>%
                  arrange(desc(AUC_het)) %>%
                  as_tibble() %>%
                  top_n(n = 30) %>%
                  melt() %>%
                  distinct() %>%
                  select(Gene,variable,value) %>% 
                  pivot_wider(names_from = variable, values_from = value) %>%
                  top_n(n = 5, wt = AUC_het) %>%
                  column_to_rownames("Gene") %>%
                  as.matrix() 

pdf("output_SCENIC/Heatmap_Regulome_LAMP5.pdf", width = 3, height = 3)
pheatmap::pheatmap(tmp,scale="none",show_rownames = T, color = viridis::viridis(n = 30, alpha = 1, begin = 0, end = 1, option = "inferno"))
dev.off()



