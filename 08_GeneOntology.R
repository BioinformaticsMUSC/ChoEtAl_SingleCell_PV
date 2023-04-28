
suppressPackageStartupMessages({
library(tidyverse)
library(ggrepel)
library(BiocParallel)
library(ggpubr)
library(magrittr)
library(broom)
library(data.table)
library(cowplot)
library(BiocSingular)
library(clusterProfiler)
library(enrichR)
})

dir.create("output_relabel/DGE_MAST/functional_enrichment/")
dge <- read.table("output_relabel/DGE_MAST/DGE_Sign_ForEnrich_HumanID.txt",header=T)
l <- split(dge, dge$Cell)
l <- l[sapply(l, nrow)>0] #remove objects with less than one gene

GOI <- list()
GeneOnto <- list()

for(i in 1:length(l)){
GOI[[i]] <- bitr(as.character(l[[i]]$Gene),  
                        fromType = "SYMBOL", 
                        toType = c("ENSEMBL", "ENTREZID"), 
                        OrgDb = org.Hs.eg.db::org.Hs.eg.db)

GeneOnto[[i]] <- enrichGO(gene = unique(GOI[[i]]$ENTREZID), 
                     keyType = "ENTREZID", 
                     OrgDb = org.Hs.eg.db::org.Hs.eg.db, 
                     ont = "BP", 
                     pAdjustMethod = "BH", 
                     pvalueCutoff  = 0.2, 
                     qvalueCutoff = 0.2, 
                     readable = TRUE)

openxlsx::write.xlsx(as.data.frame(GeneOnto[[i]]), 
                     file = sprintf("output_relabel/DGE_MAST/functional_enrichment/%s.GO.xlsx", names(l)[[i]]), 
                     colNames = TRUE,
                     rowNames = TRUE, 
                     borders = "columns",
                     sheetName="GeneOnto", 
                     overwrite = TRUE)


PLOT <- dotplot(GeneOnto[[i]])
print(PLOT)
ggsave(sprintf("output_relabel/DGE_MAST/functional_enrichment/%s_dotplot.pdf", names(l)[[i]]))


PLOT_2 <- barplot(GeneOnto[[i]])
print(PLOT_2)
ggsave(sprintf("output_relabel/DGE_MAST/functional_enrichment/%s_dotplot.pdf", names(l)[[i]]))

}


# ENRICHR
dge <- read.table("output_relabel/DGE_MAST/DGE_Sign_ForEnrich_HumanID.txt",header=T)

dbs <- c('GO_Biological_Process_2018','GO_Cellular_Component_2018',
         'GO_Molecular_Function_2018')


collapsed_output <- data.frame()
for(cur in as.character(unique(dge$Cell))){
  print(cur)
  # select genes
  cur_genes <- dge %>%
    subset(Cell == cur) %>% .$Gene

  # run enrichR on different gene sets:
  cur_result <- enrichr(cur_genes, dbs)

  # collapse results into one dataframe
  for(db in dbs){
    cur_result[[db]]$cluster <- cur
    cur_result[[db]]$db <- db
    cur_result[[db]]$Diagnosis <- 'HET'
    collapsed_output <- rbind(collapsed_output, cur_result[[db]])
  }
}

collapsed_output %>%
  write.csv(file='output_relabel/DGE_MAST/functional_enrichment/ENRICHR_DEGs_GO_terms.csv')


input_bub <- collapsed_output %>% 
    filter(cluster %in% c("Inh_Pvalb","Inh_Sst","Inh_Vip","Inh_Lamp5"),db == "GO_Biological_Process_2018") %>% 
    mutate(log = -log10(P.value)) %>%
    group_by(cluster) %>%
    top_n(3,log) %>%
    mutate(Term2 = gsub("\\s*(\\([^()]*(?:(?1)[^()]*)*\\))", "", Term, perl=TRUE)) %>%    
    mutate(Term2 = as.factor(Term2), cluster = as.factor(cluster))


colors <- c("#0D0887FF", "#6A00A8FF", "#B12A90FF","#E16462FF", "#FCA636FF", "#F0F921FF")

pdf("output_relabel/DGE_MAST/functional_enrichment/ENRICHR_DEGs_GO_bubblechart.pdf", width = 4, height = 5)
ggballoonplot(input_bub, x = "cluster", y = "Term2",
              size = "log", fill = "log") +
   scale_fill_gradientn(colors = colors) +
     scale_y_discrete(labels = function(x) str_wrap(x, width = 35)) +
  guides(size = FALSE)
dev.off()