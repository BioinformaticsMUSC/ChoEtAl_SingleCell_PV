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
library(WGCNA)
library(hdWGCNA)
library(igraph)
library(patchwork)
library(JASPAR2020)
library(motifmatchr)
library(TFBSTools)
library(EnsDb.Mmusculus.v79)
library(GeneOverlap)
library(enrichR)
})
source("utils/Utils.R")

dir.create("output_hdWGCNA")

load("output/JenniferCho_SeuratObj_Final_NoDoublet.RData")

DefaultAssay(seuObject_slim_nodoub) <- "RNA"

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



seuObject_slim_nodoub@meta.data <- seuObject_slim_nodoub@meta.data %>%
  mutate(Class = case_when(grepl("Exc", Cell) ~ "Excitatory", 
                       grepl("Inh", Cell) ~ "Inhibitory",
                       grepl("Astro|VLMC|Endo|Oligo|OPC|Micro|Peri", Cell) ~ "Glia"))


##################
## Set Up WGCNA ##
##################

wgcna_obj <- SetupForWGCNA(
  seuObject_slim_nodoub,
  gene_select = "fraction", # the gene selection approach
  fraction = 0.05, # fraction of cells that a gene needs to be expressed in order to be included
  wgcna_name = "Mef2c_Het" # the name of the hdWGCNA experiment
)

wgcna_obj <- MetacellsByGroups(
  seurat_obj = wgcna_obj,
  group.by = c("Class","Genotype"), # specify the columns in seurat_obj@meta.data to group by
  reduction = 'umap', # select the dimensionality reduction to perform KNN on
  k = 25, # nearest-neighbors parameter
  max_shared = 10, # maximum number of shared cells between two metacells
  ident.group = 'Class',
  min_cells = 30 # set the Idents of the metacell seurat object
)

Inh_obj <- SetDatExpr(
  wgcna_obj,
  group_name = "Inhibitory", # the name of the group of interest in the group.by column
  group.by='Class', # the metadata column containing the cell type info. This same column should have also been used in MetacellsByGroups
  assay = 'RNA', # using RNA assay
  slot = 'data' # using normalized data
)

# Power Test
Inh_obj <- TestSoftPowers(
  Inh_obj,
  networkType = 'signed',# you can also use "unsigned" or "signed hybrid"
  corFnc = "bicor"
)

# plot the results:
pdf("output_hdWGCNA/Power_Plots.pdf",width=6,height=6)
plot_list <- PlotSoftPowers(Inh_obj)
patchwork::wrap_plots(plot_list, ncol=2)
dev.off()

# Blockwise Module Analysis
Inh_obj <- ConstructNetwork(
  Inh_obj, 
  soft_power=12,
  corType = "bicor",
  deepSplit = 4,
  mergeCutHeight = 0.10,
  networkType = "signed",
  setDatExpr=FALSE,
  tom_name = 'Inh',
  overwrite_tom = TRUE# name of the topoligical overlap matrix written to disk
)

pdf("output_hdWGCNA/Dendrogram.pdf",width=6,height=4)
PlotDendrogram(Inh_obj, main='INH Dendrogram')
dev.off()


# Compute TOM 
Inh_Tom <- GetTOM(Inh_obj)
Inh_obj <- ScaleData(Inh_obj, features=VariableFeatures(Inh_obj),model.use = "linear")


# compute all MEs in the full single-cell dataset
Inh_obj <- ModuleEigengenes(
 Inh_obj,
 group.by.vars="Genotype",
 exclude_grey = FALSE
)

# harmonized module eigengenes:
hMEs <- GetMEs(Inh_obj)
MEs <- GetMEs(Inh_obj, harmonized=FALSE)

# compute eigengene-based connectivity (kME):
Inh_obj <- ModuleConnectivity(
  Inh_obj,
  group.by = 'Class', 
  group_name = 'Inhibitory',
  corFnc = "bicor",
  assay = 'RNA', # using RNA assay
  slot = 'data' # using normalized data
)

# rename the modules
Inh_obj <- ResetModuleNames(
  Inh_obj,
  new_name = "INH-M"
)

pdf("output_hdWGCNA/kMEs_Modules.pdf",width=6,height=4)
PlotKMEs(Inh_obj, ncol=4)
dev.off()


# compute gene scoring for the top 25 hub genes by kME for each module
# with Seurat method
Inh_obj <- ModuleExprScore(
  Inh_obj,
  n_genes = 25,
  method='Seurat'
)

plot_list <- ModuleFeaturePlot(
  Inh_obj,
  features='hMEs', # plot the hMEs
  order=TRUE # order so the points with highest hMEs are on top
)

pdf("output_hdWGCNA/Modules_UMPAS.pdf",width=8,height=8)
patchwork::wrap_plots(plot_list, ncol=3)
dev.off()

# Run UMAP for the modules and visualize it
Inh_obj <- RunModuleUMAP(
  Inh_obj,
  n_hubs = 25, # number of hub genes to include for the UMAP embedding
  n_neighbors=30, # neighbors parameter for UMAP
  min_dist=0.2 # min distance between points in UMAP space
)

pdf("output_hdWGCNA/Module_UMAP_Network.pdf",width=8,height=8)
ModuleUMAPPlot(
  Inh_obj,
  edge.alpha=0.25,
  sample_edges=TRUE,
  edge_prop=0.1, # proportion of edges to sample (20% here)
  label_hubs=2 ,# how many hub genes to plot per module?
  keep_grey_edges=FALSE
)
dev.off()

# Save the modules
ModuleNetworkPlot(Inh_obj,outdir ="output_hdWGCNA/ModuleNetworks")


#################################
# Differential Network Analysis #
#################################
p <- t(MEs)

pd <- Inh_obj@meta.data %>%
      dplyr::filter(Class == "Inhibitory")

p <- p[,colnames(p) %in% rownames(pd)]


design <- model.matrix(~Genotype,pd) # describe model to be fit
fitLM <- lmFit(p, design)

# Check only condition
cont_1 <- contrasts.fit(fitLM, coef = 2) # Directly test second coefficient
fitEb_1 <- eBayes(cont_1,robust = TRUE)


#
FullTab = topTable(fitEb_1, coef = "GenotypeWT_F",number=nrow(p)) %>%
             rownames_to_column("Module") %>%
             mutate(log = -log10(adj.P.Val))

FullTab <- FullTab %>% 
            filter(Module != "grey") %>% 
            mutate(Module = as.factor(Module)) %>%
            distinct()


pdf("output_hdWGCNA/Vulcano_Plot_DiffMEs.pdf",width=4,height=4,useDingbats=FALSE)
ggscatter(FullTab, 
            x = "logFC", 
            y = "log",
            color = "Module",
            palette=c("black","blue","brown","green","pink","red","turquoise","yellow"),
            size = "logFC",
            alpha=0.3,
            shape=19,label = "Module", repel = TRUE)+
      xlab("log2(Fold Change)")+ 
      ylab("-log10(FDR)")+
      geom_vline(xintercept = 0, colour = "grey",linetype="dotted",size=1,alpha=0.5) + 
      theme(legend.position="none")+
      ylim(0,14) + xlim(-2,+2)
dev.off()


###########################
# Correlation with TRAITS #
###########################
cur_traits <- c('Genotype', 'pMito')
Inh_obj$Genotype <- as.factor(Inh_obj$Genotype)

Inh_obj <- ModuleTraitCorrelation(
  Inh_obj,
  traits = cur_traits,
  group.by='Cell'
)


pdf("output_hdWGCNA/CorrTrait_Plot_DiffMEs.pdf",width=4,height=6,useDingbats=FALSE)
PlotModuleTraitCorrelation(
  Inh_obj,
  label = 'fdr',
  label_symbol = 'stars',
  text_size = 2,
  text_digits = 2,
  text_color = 'white',
  high_color = 'yellow',
  mid_color = 'black',
  low_color = 'purple',
  plot_max = 0.2,
  combine=TRUE
)
dev.off()



#################
# TF Enrichment #
#################
pfm_core <- TFBSTools::getMatrixSet(
  x = JASPAR2020,
  opts = list(collection = "CORE", tax_group = 'vertebrates', all_versions = FALSE)
)

# run the motif scan with these settings for the mouse dataset
Inh_obj <- MotifScan2(
  Inh_obj,
  species_genome = 'mm10',
  pfm = pfm_core,
  EnsDb = EnsDb.Mmusculus.v79
)


target_genes <- GetMotifTargets(Inh_obj)
Inh_obj<- OverlapModulesMotifs(Inh_obj)


MotifOverlapBarPlot(
  Inh_obj,
  motif_font = 'xkcd_regular',
  outdir = 'output_hdWGCNA/MotifOverlaps/',
  plot_size=c(5,6)
)

Inh_obj <- MotifTargetScore2(
  Inh_obj,
  method='Seurat'
)


df <- GetMotifOverlap(Inh_obj)
cur_df <- df %>% 
          subset(tf == 'MEF2C') %>%
          mutate(Significance = stars.pval(pval))


pdf("output_hdWGCNA/MEF2C_Binding_Barplot.pdf",width=3,height=4,useDingbats=FALSE)
plot_var <- 'odds_ratio'
cur_df %>%
  ggplot(aes(y=reorder(module, odds_ratio), x=odds_ratio)) +
  geom_bar(stat='identity', fill=cur_df$color) +
  geom_vline(xintercept = 1, linetype='dashed', color='gray') +
  geom_text(aes(label=Significance), color='black', size=3.5, hjust='center') +
  ylab('') +
  xlab("Odds Ratio") +
  ggtitle("MEF2C overlap") +
  theme(
    plot.title = element_text(hjust = 0.5)
  )+ 
  theme_classic() +
  xlim(0,2) 
dev.off()


#################
## Bubble Plot ##
#################

Inh_obj <- RenameIdents(Inh_obj, 'Inh_Adarb2' = 'Inh_Lamp5')
Inh_obj@meta.data <- Inh_obj@meta.data %>%
mutate(Cell = str_replace(Cell, "Inh_Adarb2", "Inh_Lamp5"))

order <- c("Exc_L2_Cux2","Exc_L3_Calb1","Exc_L4_Ntng1","Exc_L5_Il1rapl2","Exc_L6_Foxp2","Exc_L6_Fezf2",
  "Inh_Pvalb","Inh_Sst","Inh_Vip","Inh_Lamp5",
  "Astro_Gja1","Astro_Id3","Oligo_Mobp","ImmOligo_Bcas1","OPC_Vcan","Micro_C1qc","Micro_Mrc1",
  "Peri_Vtn","VLMC_Slc7a11","VLMC_Cped1","Endo_Flt1","Endo_Cldn5")

Idents(Inh_obj) <- factor(Idents(Inh_obj), levels = order)

Inh_obj@meta.data$Cell <- factor(Inh_obj@meta.data$Cell, levels = order)


MEs <- GetMEs(Inh_obj, harmonized=TRUE)
mods <- colnames(MEs); mods <- mods[mods != 'grey']
Inh_obj@meta.data <- cbind(Inh_obj@meta.data, MEs)

# Plot INH-M4 hME using Seurat VlnPlot function
p <- DotPlot(Inh_obj, features=mods, group.by = 'Cell')

# flip the x/y axes, rotate the axis labels, and change color scheme:
p <- p +
  coord_flip() +
  RotatedAxis() +
  scale_color_gradient2(high='red', mid='grey95', low='blue')

pdf("output_hdWGCNA/Module_Expression_Bubble.pdf",width=8,height=4,useDingbats=FALSE)
print(p)
dev.off()




# Save files! 
modules <- GetModules(Inh_obj)
mt_cor <- GetModuleTraitCorrelation(Inh_obj)


save(Inh_obj,mt_cor, modules, file='output_hdWGCNA/hdWGCNA_Inh.RData')



