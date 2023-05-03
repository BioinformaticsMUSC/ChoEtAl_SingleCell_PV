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
library(schex)
library(ggradar)
library(scGate)
library(UCell)
library(UpSetR)
library(Nebulosa)
})

source("utils/Utils.R")

dir.create("output_Figure1")

#dge <- read.table("output_sct_reclust/output_mast_manual/SalineCPPvsSalineNPAS4/Final_DEGS_CocaineCPPvsCocaineNPAS4.txt",header=T,sep="\t")
#full <- read.table("output_sct_reclust/output_mast_manual/SalineCPPvsSalineNPAS4/Final_FullStats_CocaineCPPvsCocaineNPAS4.txt",header=T,sep="\t")

load("output_relabel/JenniferCho_SeuratObj_Final_NoDoublet_Relabel.RData")
load("output_relabel/DGE_MAST/JennCho_DGE.RData")

# Order the cells

order <- c("Exc_L2_Cux2","Exc_L3_Calb1","Exc_L3_Ntng1","Exc_L5_Il1rapl2","Exc_L6_Foxp2","Exc_L6_Fezf2",
  "Inh_Pvalb","Inh_Sst","Inh_Vip","Inh_Lamp5",
  "Astro_Gja1","Astro_Id3","Oligo_Mobp","ImmOligo_Bcas1","OPC_Vcan","Micro_C1qc","Micro_Mrc1",
  "Peri_Vtn","VLMC_Slc7a11","VLMC_Cped1","Endo_Flt1","Endo_Cldn5")

# From VanGogh2 MetBrewer
colore <- c("#bd3106","#454b87")


meta <- as.data.frame(seuObject_slim_nodoub@meta.data)


# Additional Quality Check

pdf("output_relabel/Quality_Control_Genes_Cells.pdf",width=8,height=4,useDingbats=FALSE)
meta %>%
mutate(Cell=fct_relevel(Cell,order)) %>%
ggviolin("Cell", "nFeature_RNA", fill = "Cell",
   add = "boxplot", add.params = list(fill = "white"))  +
        theme(legend.position="none") +
        rotate_x_text(angle = 45) +
        xlab("") + 
        ylab("Genes Detected")
dev.off()


pdf("output_relabel/Quality_Control_UMI_Cells.pdf",width=8,height=4,useDingbats=FALSE)
meta %>%
mutate(Cell=fct_relevel(Cell,order)) %>%
ggviolin("Cell", "nCount_RNA", fill = "Cell",
   add = "boxplot", add.params = list(fill = "white"))  +
        theme(legend.position="none") +
        rotate_x_text(angle = 45) +
        xlab("") + 
        ylab("UMI")
dev.off()


## Visualize number of Mef2c Expressing cells 
Mef2c = GetAssayData(object=seuObject_slim_nodoub,assay="RNA",slot="data")["Mef2c",]

seuObject_slim_nodoub@meta.data$Mef2c <- as.numeric(log2(Mef2c+1))

umap <- as.data.frame(Embeddings(seuObject_slim_nodoub, reduction = "umap"))
meta <- as.data.frame(seuObject_slim_nodoub@meta.data)

# Npas4 expression violin plot
pdf("output_Figure1/Mef2c_ViolinPlot.pdf", width = 6, height = 6)
meta %>%
filter(Mef2c > 0) %>%
mutate(Cell=fct_relevel(Cell,order)) %>%
ggboxplot(x = "Genotype", y = "Mef2c",
            color = "Genotype", 
            palette=c("#bd3106","#454b87"),
            outlier.shape = NA)+ 
        stat_compare_means(comparisons = list(c("HET_F", "WT_F")),label.y = c(2),method = "wilcox.test") +
        theme(legend.position="none") +
        rotate_x_text(angle = 45) +
        xlab("") + 
        ylab("Expression") +
        facet_wrap(.~Cell)
dev.off()


# Cell number total
colore <- c("#bd3106","#454b87")
grouped_data <- table(meta$Cell, meta$Genotype) %>%
                as.matrix() %>%
                as.data.frame() %>%
                arrange(Var1,Var2) %>%
                dplyr::rename(Cell = Var1, Genotype = Var2, Cells = Freq)  %>%
                mutate(Cell=fct_rev(fct_relevel(Cell,order))) %>%
                as.data.frame() %>%
                droplevels()

pdf("output_Figure1/Barplot_NumberOfCells.pdf", width = 4, height = 5)
ggbarplot(grouped_data, "Cell", "Cells",
  fill = "Genotype", color = "Genotype", palette = colore,
  label = FALSE, lab.col = "white", lab.pos = "in",
  rotate = TRUE) + 
#scale_y_continuous(labels = dollar_format(suffix = "%", prefix = "")) +
rotate_x_text(angle = 45) +
theme(legend.position="top") +
xlab("")
dev.off()

# Mef2c only
df <- cbind(umap,meta)%>% 
    filter(Mef2c > 0)

grouped_data <- table(df$Cell, df$Genotype) %>%
                as.matrix() %>%
                as.data.frame() %>%
                arrange(Var1,Var2) %>%
                dplyr::rename(Cell = Var1, Genotype = Var2, Cells = Freq)  %>%
                mutate(Cell=fct_rev(fct_relevel(Cell,order))) %>%
                as.data.frame() %>%
                droplevels()

pdf("output_Figure1/Barplot_NumberOfCells_Mef2c.pdf", width = 4, height = 5)
ggbarplot(grouped_data, "Cell", "Cells",
  fill = "Genotype", color = "Genotype", palette = colore,
  label = FALSE, lab.col = "white", lab.pos = "in",
  rotate = TRUE) + 
#scale_y_continuous(labels = dollar_format(suffix = "%", prefix = "")) +
rotate_x_text(angle = 45) +
theme(legend.position="top") +
xlab("")
dev.off()


# Radarplot for significant
tmp1 <- dge %>% 
        group_by(cell_type) %>%
        mutate(LOG = -log10(p_val_adj), ABS = abs(avg_logFC)) %>% 
        mutate(Threshold = if_else(p_val_adj < 0.05 & ABS > 0.2, "TRUE","FALSE")) %>%
        mutate(Direction = case_when(avg_logFC > 0.2 & p_val_adj < 0.05 ~ "UpReg", avg_logFC < -0.2 & p_val_adj < 0.05 ~ "DownReg")) %>%
        filter(Threshold == "TRUE") 

df1 <- table(tmp1$cell_type) %>%
    as.data.frame() %>% 
    pivot_wider(names_from = Var1, values_from = Freq) %>%
    as.data.frame()

df1$Class <- "Mef2c_Het"
labels <- c("Class",order)

df1 <- df1[,match(labels,colnames(df1))]

pdf("output_Figure1/RadarPlot_Degs.pdf", width = 8, height = 8)
ggradar::ggradar(
  df1, 
  values.radar = c("0", "400", "800"),
  grid.min = 0, 
  grid.mid = 400, 
  grid.max = 800,
  group.line.width = 1, 
  group.point.size = 3,
  group.colours = c("#bd3106","#454b87"),
  background.circle.colour = "white",
  gridline.mid.colour = "grey",
  legend.position = "bottom"
  )
dev.off()

# Barplot all genes
tmp1 <- dge %>% 
        group_by(cell_type) %>%
        mutate(LOG = -log10(p_val_adj), ABS = abs(avg_logFC)) %>% 
        mutate(Threshold = if_else(p_val_adj < 0.05 & ABS > 0.2, "TRUE","FALSE")) %>%
        mutate(Direction = case_when(avg_logFC > 0.2 & p_val_adj < 0.05 ~ "UpReg", avg_logFC < -0.2 & p_val_adj < 0.05 ~ "DownReg")) %>%
        filter(Threshold == "TRUE") 

df1 <- table(tmp1$cell_type) %>%
    as.data.frame() %>%
    mutate(Var1=fct_rev(fct_relevel(Var1,order)))


df1 <- df1[match(order,df1$Var1),]

pdf("output_Figure1/Barplot_Number_DEGS.pdf", width = 4, height = 5)
ggbarplot(df1, "Var1", "Freq",
  fill = "#404788FF", color = "#404788FF",
  label = FALSE, lab.col = "white", lab.pos = "in",
  rotate = TRUE) + 
#scale_y_continuous(labels = dollar_format(suffix = "%", prefix = "")) +
rotate_x_text(angle = 45) +
theme(legend.position="top") +
xlab("") +
ylim(0,1000)+
ylab("")
dev.off()

# Radarplot for neurons
df1 <- table(tmp1$cell_type) %>%
    as.data.frame() %>% 
    pivot_wider(names_from = Var1, values_from = Freq) %>%
    as.data.frame()

df1$Class <- "Mef2c_Het"
labels <- c("Class",order)

df1 <- df1[,match(labels,colnames(df1))]

df1 <- df1[,c(1:11)]

pdf("output_Figure1/RadarPlot_Degs_OnlyNeurons.pdf", width = 5, height = 5)
ggradar(
  df1, 
  values.radar = c("0", "200", "400"),
  grid.min = 0, grid.mid = 200, grid.max = 400,
  # Polygons
  group.line.width = 1, 
  group.point.size = 3,
  group.colours = c("#da7901","#f0be3d","#da7901","#247d3f"),
  # Background and grid lines
  background.circle.colour = "white",
  gridline.mid.colour = "grey",
  legend.position = "bottom"
  )
dev.off()


# Upset plot
tmp_upset <- tmp1 %>%
    as.data.frame() %>%
    filter(cell_type %in% c("Inh_Pvalb","Inh_Sst","Inh_Vip","Inh_Lamp5")) %>%
    select(gene, cell_type) %>%
    droplevels()

pdf("output_Figure1/Upset_Plot_Inhibitory_Neurons.pdf", width = 6, height = 4)
l <- split(as.character(tmp_upset$gene),tmp_upset$cell_type)
Class <- names(l)
ToTGene <- as.numeric(sapply(l, length))
metadata <- as.data.frame(cbind(Class, ToTGene))
names(metadata) <- c("Class", "ToTGene")
metadata$ToTGene <- as.numeric(as.character(metadata$ToTGene))
upset(fromList(l),,nsets = 6, set.metadata = list(data = metadata, plots = list(list(type = "hist", 
    column = "ToTGene", assign = 20), 
    list(type = "matrix_rows", column = "sets", colors = c(Inh_Pvalb = "#89E651", Inh_Sst = "#DB8BE4",Inh_Vip = "#DDA187", Inh_Lamp5="#DDD3DD"), 
    alpha = 0.5))))
dev.off()



###########################################
            # Vulcano plots  #
###########################################
df1 <- dge %>% 
        group_by(cell_type) %>%
        mutate(LOG = -log10(p_val_adj), ABS = abs(avg_logFC)) %>% 
        mutate(Threshold = if_else(p_val_adj < 0.05 & ABS > 0.2, "TRUE","FALSE")) %>%
        mutate(Direction = case_when(avg_logFC > 0.2 & p_val_adj < 0.05 ~ "UpReg", avg_logFC < -0.2 & p_val_adj < 0.05 ~ "DownReg")) 

# cell colors original
col <- data.table::fread("output_relabel/Colors_Used_Initial_Data.txt") %>% as.data.frame()
col <- col[match(order,col$Cells),]

colors <- as.character(col$Colors)


pdf("output_Figure1/Vulcano_AllCells.pdf", width = 4, height = 4)
df1 %>%
ggplot(aes(x=avg_logFC, y=LOG)) +
ggrastr::geom_point_rast(aes(colour = cell_type),size=0.5) +
      xlab("log2(Fold Change)")+ 
      ylab("-log10(FDR)")+
      geom_vline(xintercept = 0, colour = "grey",linetype="dotted",size=1,alpha=0.5) +
      geom_vline(xintercept = 0.2, colour = "grey",linetype="dashed",size=1,alpha=0.5) +
      geom_vline(xintercept = -0.2, colour = "grey",linetype="dashed",size=1,alpha=0.5) + 
      geom_hline(yintercept = 1.3, colour = "grey",linetype="dotted",size=1,alpha=0.5) +
      theme(legend.position="none")+
       scale_colour_manual(values = colors) +   # customized color palette
      ylim(0,100) + xlim(-5,5)
dev.off()

pdf("output_Figure1/Vulcano_Pvalb.pdf", width = 4, height = 4)
top_labelled <- df1 %>% 
                  filter(cell_type == "Inh_Pvalb", gene %in% c("Slit3", "Grik2", "Nrg1", "Ctnna2", "Tac1", "Chrm3","Grin2a")) 

df1 %>%
  filter(cell_type == "Inh_Pvalb") %>%
  ggscatter( 
            x = "avg_logFC", 
            y = "LOG",
            color = "Threshold",
            palette=c("grey","red"),
            size = 1,
            alpha=0.3,
            shape=19)+
      xlab("log2(Fold Change)")+ 
      ylab("-log10(FDR)")+
      geom_vline(xintercept = 0, colour = "grey",linetype="dotted",size=1,alpha=0.5) +
      geom_vline(xintercept = 0.2, colour = "grey",linetype="dashed",size=1,alpha=0.5) +
      geom_vline(xintercept = -0.2, colour = "grey",linetype="dashed",size=1,alpha=0.5) + 
      geom_hline(yintercept = 1.3, colour = "grey",linetype="dotted",size=1,alpha=0.5) +
      geom_text_repel(data = top_labelled, 
                        mapping = aes(label = gene), 
                        size = 5,
                        nudge_x = .15,
                        box.padding = 0.5,
                        nudge_y = 1,
                        segment.curvature = -0.1,
                        segment.ncp = 3,
                        segment.angle = 20)+
      theme(legend.position="none")+
      ylim(0,20) + xlim(-2,2)
dev.off()

# Density Plot MEF2C


pdf("output_Figure1/Mef2c_UMAP.pdf", width = 6, height = 6)

df <- cbind(umap,meta)

ggplot(df, aes(x=UMAP_1, y=UMAP_2,color = Mef2c)) +
    #geom_point() + 
    ggrastr::geom_point_rast(aes(colour = Mef2c),size=0.1) + 
    scale_colour_gradient2(low = "gray60", high = "red", mid="gray60",na.value = NA) +##FDE725FF #"#440154FF"
    theme_void() +
    theme(legend.position="none")

dev.off()

# ComplexHeatmap Neurons
load("output_relabel/DGE_MAST/JennCho_DGE.RData")

tmp1 <- dge %>% 
        group_by(cell_type) %>%
         mutate(LOG = -log10(p_val_adj), ABS = abs(avg_logFC)) %>% 
         mutate(Threshold = if_else(p_val_adj < 0.05 & ABS > 0.2, "TRUE","FALSE")) %>%
         mutate(Direction = case_when(avg_logFC > 0.2 & p_val_adj < 0.05 ~ "UpReg", avg_logFC < -0.2 & p_val_adj < 0.05 ~ "DownReg")) %>%
         filter(Threshold == "TRUE") %>%
         filter(grepl("Inh_|Exc_",cell_type)) %>% 
         as.data.frame()


tab<-table(tmp1$cell_type, tmp1$Direction)

mat<-matrix(nrow=length(unique(tmp1$gene)),ncol=10)
rownames(mat)<-unique(tmp1$gene)
colnames(mat)<-unique(tmp1$cell_type)

for (i in 1:nrow(mat)){
   for (j in 1:ncol(mat)){
     gene_tmp<-tmp1[which(tmp1$gene==rownames(mat)[i]),]
     gene_tmp<-gene_tmp[which(gene_tmp$cell_type==colnames(mat)[j]),]
     mat[i,j]<-ifelse(nrow(gene_tmp)>0, gene_tmp$avg_logFC,0)
   }
 }

noDup<-tmp1[!duplicated(tmp1$gene),]
mat<-mat[order(noDup$cell_type, noDup$avg_logFC),]
mat[ , colSums(is.na(mat))==0]

tab <- tab[rownames(tab) %in% colnames(mat),]

 # Heatmap
 ha2<-HeatmapAnnotation(celltype=colnames(mat), 
      col= list(celltype=c("Exc_L3_Calb1"="#DF5AD2",
                           "Exc_L5_Il1rapl2"="#847E74",
                           "Exc_L6_Foxp2"="#DCDE4A",
                            "Exc_L2_Cux2"="#E36284",
                           "Exc_L4_Ntng1"="#4DA4C3",
                           "Exc_L6_Fezf2"="#D9A6D1",
                           "Inh_Pvalb"="#89E651",
                           "Inh_Sst"="#DB8BE4",
                           "Inh_Lamp5"="#DDD3DD",
                           "Inh_Vip"="#DDA187"
                           )), 
      show_legend=F)
 
 tab <- tab[rownames(tab) %in% colnames(mat),]
 tab <- tab[match(colnames(mat), rownames(tab)),]

# Top annotaiton
 bar1<-HeatmapAnnotation(up=anno_barplot(tab[,2], gp=gpar(fill="red"), axis_param=list(labels=c("","","","")), ylim=c(0,200)),
                         down=anno_barplot((tab[,1] * -1), gp=gpar(fill="blue"), axis_param=list(labels=c("","","")), ylim=c(-200,0),
                         show_annotation_name=F,
                         show_legend=F))
 ha<-c( bar1, ha2)

# Annotate cool ASD genes
asd_genes <- read.table("utils/ASD_SC.txt",header=T,sep="\t")
lab<-which(noDup$gene %in% firstup(asd_genes$Gene)) #& abs(noDup$avg_logFC)>0.5) 



pdf("output_Figure1/DGE_Neurons_heatmap_foldCs.pdf", width=6, height=8)
ht <- Heatmap(mat, cluster_rows=F,show_row_dend=F,cluster_columns=F,col=circlize::colorRamp2(c(1,0,-1),c("blue","white","red")), 
        show_column_names=T, show_row_names=F, column_title=NULL,row_names_gp = gpar(fontsize = 3), top_annotation=ha)
rowAnnotation(gene=anno_mark(at=lab, labels=noDup[lab,]$gene, side="left", labels_gp=gpar(fontsize=7)))+ht
dev.off()


# ComplexHeatmap Inhibitory
load("output_relabel/DGE_MAST/JennCho_DGE.RData")

tmp1 <- dge %>% 
        group_by(cell_type) %>%
         mutate(LOG = -log10(p_val_adj), ABS = abs(avg_logFC)) %>% 
         mutate(Threshold = if_else(p_val_adj < 0.05 & ABS > 0.2, "TRUE","FALSE")) %>%
         mutate(Direction = case_when(avg_logFC > 0.2 & p_val_adj < 0.05 ~ "UpReg", avg_logFC < -0.2 & p_val_adj < 0.05 ~ "DownReg")) %>%
         filter(Threshold == "TRUE") %>%
         filter(grepl("Inh_",cell_type)) %>% 
         as.data.frame()


tab<-table(tmp1$cell_type, tmp1$Direction)

mat<-matrix(nrow=length(unique(tmp1$gene)),ncol=4)
rownames(mat)<-unique(tmp1$gene)
colnames(mat)<-unique(tmp1$cell_type)

for (i in 1:nrow(mat)){
   for (j in 1:ncol(mat)){
     gene_tmp<-tmp1[which(tmp1$gene==rownames(mat)[i]),]
     gene_tmp<-gene_tmp[which(gene_tmp$cell_type==colnames(mat)[j]),]
     mat[i,j]<-ifelse(nrow(gene_tmp)>0, gene_tmp$avg_logFC,0)
   }
 }

noDup<-tmp1[!duplicated(tmp1$gene),]
mat<-mat[order(noDup$cell_type, noDup$avg_logFC),]
mat[ , colSums(is.na(mat))==0]

tab <- tab[rownames(tab) %in% colnames(mat),]

 # Heatmap
 ha2<-HeatmapAnnotation(celltype=colnames(mat), 
      col= list(celltype=c(
                           "Inh_Pvalb"="#89E651",
                           "Inh_Sst"="#DB8BE4",
                           "Inh_Lamp5"="#DDD3DD",
                           "Inh_Vip"="#DDA187"
                           )), 
      show_legend=F)
 
 tab <- tab[rownames(tab) %in% colnames(mat),]
 tab <- tab[match(colnames(mat), rownames(tab)),]

# Top annotaiton
 bar1<-HeatmapAnnotation(up=anno_barplot(tab[,2], gp=gpar(fill="red"), axis_param=list(labels=c("","","","")), ylim=c(0,200)),
                         down=anno_barplot((tab[,1] * -1), gp=gpar(fill="blue"), axis_param=list(labels=c("","","")), ylim=c(-200,0),
                         show_annotation_name=F,
                         show_legend=F))
 ha<-c( bar1, ha2)

# Annotate cool ASD genes
asd_genes <- read.table("utils/ASD_SC.txt",header=T,sep="\t")
lab<-which(noDup$gene %in% DIZutils::firstup(asd_genes$Gene)) #& abs(noDup$avg_logFC)>0.5) 



pdf("output_Figure1/DGE_Inhibitory_heatmap_foldCs.pdf", width=4, height=6)
ht <- Heatmap(mat, cluster_rows=F,show_row_dend=F,cluster_columns=F,col=circlize::colorRamp2(c(1,0,-1),c("red","white","blue")), 
        show_column_names=T, show_row_names=F, column_title=NULL,row_names_gp = gpar(fontsize = 3), top_annotation=ha)
rowAnnotation(gene=anno_mark(at=lab, labels=noDup[lab,]$gene, side="left", labels_gp=gpar(fontsize=7)))+ht
dev.off()


# ComplexHeatmap for non neuronal 

tmp1 <- dge %>% 
        group_by(cell_type) %>%
         mutate(LOG = -log10(p_val_adj), ABS = abs(avg_logFC)) %>% 
         mutate(Threshold = if_else(p_val_adj < 0.05 & ABS > 0.2, "TRUE","FALSE")) %>%
         mutate(Direction = case_when(avg_logFC > 0.2 & p_val_adj < 0.05 ~ "UpReg", avg_logFC < -0.2 & p_val_adj < 0.05 ~ "DownReg")) %>%
         filter(Threshold == "TRUE") %>%
         filter(!grepl("Inh_|Exc_",cell_type)) %>% 
         as.data.frame()


tab<-table(tmp1$cell_type, tmp1$Direction)

mat<-matrix(nrow=length(unique(tmp1$gene)),ncol=10)
rownames(mat)<-unique(tmp1$gene)
colnames(mat)<-unique(tmp1$cell_type)

for (i in 1:nrow(mat)){
   for (j in 1:ncol(mat)){
     gene_tmp<-tmp1[which(tmp1$gene==rownames(mat)[i]),]
     gene_tmp<-gene_tmp[which(gene_tmp$cell_type==colnames(mat)[j]),]
     mat[i,j]<-ifelse(nrow(gene_tmp)>0, gene_tmp$avg_logFC,0)
   }
 }

noDup<-tmp1[!duplicated(tmp1$gene),]
mat<-mat[order(noDup$cell_type, noDup$avg_logFC),]
mat[ , colSums(is.na(mat))==0]

tab <- tab[rownames(tab) %in% colnames(mat),]

 # Heatmap
 ha2<-HeatmapAnnotation(celltype=colnames(mat), 
      col= list(celltype=c("Astro_Gja1"="#A5D996",
                           "Micro_C1qc"="#8259D9",
                           "Oligo_Mobp"="#D8CE7A",
                            "OPC_Vcan"="#95ADDC",
                           "VLMC_Slc7a11"="#BB37EA",
                           "Peri_Vtn"="#D9E4C7",
                           "VLMC_Cped1"="#E3853D",
                           "Endo_Flt1"="#90DBE0",
                           "Micro_Mrc1"="#65E8C8",
                           "ImmOligo_Bcas1"="#6DE688"
                           )), 
      show_legend=F)
 
 tab <- tab[rownames(tab) %in% colnames(mat),]
 tab <- tab[match(colnames(mat), rownames(tab)),]

# Top annotaiton
 bar1<-HeatmapAnnotation(up=anno_barplot(tab[,2], gp=gpar(fill="red"), axis_param=list(labels=c("","","","")), ylim=c(0,500)),
                         down=anno_barplot((tab[,1] * -1), gp=gpar(fill="blue"), axis_param=list(labels=c("","","")), ylim=c(-500,0),
                         show_annotation_name=F,
                         show_legend=F))
 ha<-c( bar1, ha2)

# Annotate cool ASD genes
asd_genes <- read.table("utils/ASD_SC.txt",header=T,sep="\t")
lab<-which(noDup$gene %in% firstup(asd_genes$Gene)) #& abs(noDup$avg_logFC)>0.5) 



pdf("output_Figure1/DGE_Glia_heatmap_foldCs.pdf", width=6, height=8)
ht <- Heatmap(mat, cluster_rows=F,show_row_dend=F,cluster_columns=F,col=circlize::colorRamp2(c(1,0,-1),c("blue","white","red")), 
        show_column_names=T, show_row_names=F, column_title=NULL,row_names_gp = gpar(fontsize = 3), top_annotation=ha)
rowAnnotation(gene=anno_mark(at=lab, labels=noDup[lab,]$gene, side="left", labels_gp=gpar(fontsize=5)))+ht
dev.off()

