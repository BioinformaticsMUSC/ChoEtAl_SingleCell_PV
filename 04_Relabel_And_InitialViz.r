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
library(scbp)
library(scPred)
library(ComplexHeatmap)
})


# create dir for initial figure 1
dir.create("output_relabel")

#plan("multiprocess", workers = 10)
options(future.globals.maxSize = 6000 * 1024^2)

load("output/JenniferCho_SeuratObj_Final_NoDoublet.RData")

labels <- read.table("output/Labels_Clusters.txt",header=T,sep="\t")

new <- labels %>% 
                separate(variable, c("Cell", "Class","Definition"),"_",remove = FALSE) %>%
                separate(Rows, c("ID", "Cluster"),"_",remove = FALSE) %>%
                unite(Ident, c("Cell","Definition"),sep = "_",remove = FALSE,na.rm = TRUE) %>%
                unite(Ident2, c("Cell","Cluster","Definition"),sep = "_",remove = FALSE,na.rm = TRUE)


new <- new[order(as.numeric(as.character(new$Ident2))), ]

current.cluster.ids <- new$Cluster
new.cluster.ids <- as.character(new$variable)

seuObject_slim_nodoub@active.ident <- plyr::mapvalues(x = seuObject_slim_nodoub@active.ident, 
                                                     from = current.cluster.ids, 
                                                     to = new.cluster.ids)

# UMAP plot with new labels
seuObject_slim_nodoub@meta.data$Cell <- seuObject_slim_nodoub@active.ident

save(seuObject_slim_nodoub,file = "output_relabel/JenniferCho_SeuratObj_Final_NoDoublet_Relabel.RData")


umap <- as.data.frame(Embeddings(seuObject_slim_nodoub, reduction = "umap"))
meta <- as.data.frame(seuObject_slim_nodoub@meta.data)

df <- cbind(umap,meta)%>% 
  group_by(Cell) %>% 
  mutate(N = n()) %>%
    ungroup() %>% 
    arrange(Genotype) %>%
    as.data.frame()  


label <- data.frame(Cell=levels(df$Cell),label=levels(df$Cell))

label_2 <- df %>% 
  group_by(Cell) %>% 
  summarize(UMAP_1 = median(UMAP_1), UMAP_2 = median(UMAP_2),N = n()) %>% 
  left_join(label) %>%
  as.data.frame() 

col <- randomcoloR::distinctColorPalette(22)
write.table(col, "output_relabel/Colors_Used_Initial_Data.txt", sep="\t",quote=F)
#col <- data.table::fread("output_sct_reclust/Figure1/Colors_Used_Initial_Data.txt")
#col <- as.character(col$Colors)

pdf("output_relabel/Data_Integrated_UMAP_Labels_Alternative.pdf", width = 8, height = 8)
#n <- 60
#qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',] 
#col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
ggplot(df, aes(x=UMAP_1, y=UMAP_2)) +
ggrastr::geom_point_rast(aes(colour = Cell),size=0.5) +
ggrepel::geom_text_repel(data = label_2, aes(label = label),
              color = "black",
              #fontface = 'bold',
              segment.colour = "grey60",
                box.padding = unit(0.25, "lines"),
                point.padding = unit(0.5, "lines"),
                nudge_x = .15,
                nudge_y = 1,
                size = 6) + 
    #scale_color_viridis(discrete=TRUE,option="inferno")
    #scale_colour_manual(values = colorRampPalette(brewer.pal(13, "Paired"))(31))
    #scale_colour_manual(values = met.brewer(name="Klimt",n=13))+    
    scale_colour_manual(values = col)+
    theme_void() +
    theme(legend.position="none")

dev.off()

# Color table
colors_used <- data.frame(Cells = unique(df$Cell), Colors = col)
write.table(colors_used, "output_relabel/Colors_Used_Initial_Data.txt",sep="\t",quote=F)

# Proportion by Genotype

order <- c("Exc_L2_Cux2","Exc_L3_Calb1","Exc_L3_Ntng1","Exc_L5_Il1rapl2","Exc_L6_Foxp2","Exc_L6_Fezf2",
  "Inh_Pvalb","Inh_Sst","Inh_Vip","Inh_Lamp5",
  "Astro_Gja1","Astro_Id3","Oligo_Mobp","ImmOligo_Bcas1","OPC_Vcan","Micro_C1qc","Micro_Mrc1",
  "Peri_Vtn","VLMC_Slc7a11","VLMC_Cped1","Endo_Flt1","Endo_Cldn5")

# From VanGogh2 MetBrewer
colore <- c("#bd3106","#454b87")


prop_cell <-  meta %>%
  as.data.frame() %>%
  mutate(Cell=fct_relevel(Cell,order)) %>%
  arrange(Cell) %>%
    group_by(Cell,Genotype) %>% 
    summarise(cnt = n()) %>% 
    group_by(Cell) %>%  
     mutate(percent = 100*round(cnt/sum(cnt), 2)) %>% 
     group_by(Cell) %>%
      ggbarplot("Cell", "percent",
    fill = "Genotype", color = "Genotype", palette = colore,
    label = TRUE,lab.pos = "in", lab.col = "white") +
  scale_y_continuous(labels = dollar_format(suffix = "%", prefix = "")) +
rotate_x_text(angle = 45) +
theme(legend.position="right") +
xlab("")
ggsave("output_relabel/CellProportion_By_Genotype.pdf", plot = prop_cell, width = 9, height = 5, units = "in", dpi = 150)

# UMAP by Genotype
pdf("output_relabel/Data_Integrated_ByGenotype.pdf", width = 6, height = 5)
#n <- 60
#qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',] 
#col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
ggplot(df, aes(x=UMAP_1, y=UMAP_2)) +
ggrastr::geom_point_rast(aes(colour = Genotype),size=0.5) + 
    #scale_color_viridis(discrete=TRUE,option="inferno")
    #scale_colour_manual(values = colorRampPalette(brewer.pal(13, "Paired"))(31))
    #scale_colour_manual(values = met.brewer(name="Klimt",n=13))+    
    scale_colour_manual(values = colore)+
    theme_void() +
    theme(legend.position="none") +
    facet_wrap(.~Genotype)
dev.off()


# Dot plot
Idents(seuObject_slim_nodoub) <- factor(Idents(seuObject_slim_nodoub), levels = order)

pdf("output_relabel/Dot_Plot_Markers.pdf", width = 7, height = 6)


cd_genes <- c("Cux2","Calb1","Ntng1","Il1rapl2","Foxp2","Fezf2",
              "Pvalb","Sst","Vip","Lamp5",
              "Gja1","Id3",
              "Mobp","Bcas1","Vcan",
              "C1qc","Vtn",
              "Slc7a11","Cped1","Flt1",
              "Snap25","Elavl2","Syt1")
DotPlot(object = seuObject_slim_nodoub, features = cd_genes) + 
      theme(axis.text.x = element_text(angle = 45, hjust=1)) + 
          theme(legend.position="top") 
dev.off()


# Do heatmap
markers<- presto::wilcoxauc(seuObject_slim_nodoub, 'Cell', assay = 'data')
markers<- presto::top_markers(markers, n = 20, auc_min = 0.5, pct_in_min = 20, pct_out_max = 20)

all_markers<- markers %>%
  select(-rank) %>% 
  unclass() %>% 
  stack() %>%
  pull(values) %>%
  unique() %>%
   .[!is.na(.)]

cols <- viridis::viridis(100)[c(1, 50, 100)]

downsamples <- subset(seuObject_slim_nodoub, downsample = 300)


mat<- downsamples[["RNA"]]@data[all_markers, ] %>% as.matrix()
mat<- t(scale(t(mat)))
cluster_anno<- downsamples@meta.data$Cell

cluster_anno <- factor(cluster_anno, levels = order)

# Colors
col_fun = circlize::colorRamp2(c(-1, 0, 3), c("#440154FF", "#2D708EFF", "#FDE725FF"))

pdf("output_relabel/heatmap_markers.pdf", width = 8, height = 3)
ComplexHeatmap::Heatmap(mat, name = "Expression",  
        column_split = cluster_anno,
        cluster_columns = FALSE,
        show_column_dend = FALSE,
        cluster_column_slices = TRUE,
        column_title_gp = grid::gpar(fontsize = 8),
        column_gap = unit(0.5, "mm"),
        cluster_rows = FALSE,
        show_row_dend = FALSE,
        col = col_fun,
        row_names_gp = grid::gpar(fontsize = 4),
        column_title_rot = 90,
        top_annotation = ComplexHeatmap::HeatmapAnnotation(foo = anno_block(gp = grid::gpar(fill = col))),
        show_column_names = FALSE,
        use_raster = TRUE,
        show_row_names = FALSE,
        raster_quality = 4)
dev.off()

