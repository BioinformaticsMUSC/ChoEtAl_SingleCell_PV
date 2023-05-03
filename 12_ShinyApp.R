library(rsconnect)
library(ShinyCell)
library(Seurat)
library(tidyverse)


load("output_relabel/JenniferCho_SeuratObj_Final_NoDoublet_Relabel.RData")


seu_filt <- seuObject_slim_nodoub_slim

seu_filt@meta.data <- seu_filt@meta.data %>% 
                            select(Genotype, pMito, nCount_RNA, nFeature_RNA, seurat_clusters, Cell)


scConf = createConfig(seu_filt)
makeShinyApp(seu_filt, scConf, gene.mapping = TRUE,
             shiny.title = "Cho et al, 2023 - PFC, Mef2c HET") 

rsconnect::setAccountInfo(name='bioinformatics-musc',
              token='838A2925138D0715F2D093E909823204',
              secret='suQLpnyt0hEXRa+LXT9IQxaYNCSIITnH2vBymGPJ')

options(repos = BiocManager::repositories())

rsconnect::deployApp("shinyApp/",
    appName = 'Cho_PFC_Mef2cHet', 
    account = 'bioinformatics-musc', 
    server = 'shinyapps.io')
