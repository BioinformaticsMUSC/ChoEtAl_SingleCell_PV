library(WGCNA);
library(cluster);
library(ggplot2)
library(reshape2)
library(RColorBrewer)
library(dplyr)
library(magrittr)
library(tidyverse)
library(ggpubr)
enableWGCNAThreads()

wgcna = list.files(path = "output/",pattern = 'Markers_Significant_DGE.txt',full.names = TRUE)
tab=read.table(wgcna,sep="\t",header=T)
colnames(tab)=c("Gene","DEFINITION")
tab$DEFINITION <- as.factor(tab$DEFINITION)
Genes=as.data.frame(table(tab$DEFINITION))

# Loop to make the overlap
# The loop goes along the length of the GeneSets lists
load("utils/GeneSets_Allen_Mouse_PFC.RData")

ln=length(GeneSets)
cl=length(Genes$Var1)
TEMP=list()
INT=list()
for (i in 1:ln)
{
TEMP[[i]]=tab[tab$Gene %in% GeneSets[[i]]$Gene,]
INT[[i]]=as.data.frame(table(TEMP[[i]]$DEFINITION))
}
names(INT)=names(GeneSets)
names(TEMP)=names(GeneSets)

# Create the matrix for each GeneSet
NROWS <- sapply(GeneSets,nrow)

#
#               GeneSet
#              +       -
#          +-------+-------+
#        + |   a   |   b   |
#  Module  +-------+-------+
#        - |   c   |   d   |
#          +-------+-------+

for (i in 1:length(INT))
{
INT[[i]]$b <- NROWS[[i]]-INT[[i]]$Freq
INT[[i]]$c <- Genes$Freq-INT[[i]]$Freq 
INT[[i]]$d <- 15585-Genes$Freq-nrow(GeneSets[[i]]) #19776 #15585
}

# sum(Genes$Freq)
RunFisher <- function(row, alt = 'greater', cnf = 0.85) {
  f <- fisher.test(matrix(row, nrow = 2), alternative = alt, conf.level = cnf)
  return(c(row,
           P_val = f$p.value,
           LogP = -log10(f$p.value), 
           OR = f$estimate[[1]],
           OR_Low = f$conf.int[1],
           OR_Up = f$conf.int[2]))
}

# run
FisherMat=list()
for (i in 1:length(INT))
{
FisherMat[[i]] <- t(apply(INT[[i]][,2:5], 1, RunFisher))
rownames(FisherMat[[i]]) <- INT[[i]]$Var1
FisherMat[[i]] <- FisherMat[[i]][rownames(FisherMat[[i]]) != "grey",]
}
names(FisherMat)<-names(INT)
#save(FisherMat,TEMP,file="Allen_MTG_Stats.RData")

# Create matrices of Pval
tmp<-list()
FisherP<-matrix()
rowNames <- rownames(FisherMat[[1]])
colNames <- names(FisherMat)
for (i in 1:length(INT))
{
tmp[[i]] <- cbind(as.numeric(FisherMat[[i]][,5]))
FisherP <- do.call(cbind,tmp)
}
rownames(FisherP) <- rowNames
colnames(FisherP) <- colNames

# Create matrices of OR
tmp<-list()
FisherOR<-matrix()
rowNames <- rownames(FisherMat[[1]])
colNames <- names(FisherMat)
for (i in 1:length(INT))
{
tmp[[i]] <- cbind(as.numeric(FisherMat[[i]][,7]))
FisherOR <- do.call(cbind,tmp)
}
rownames(FisherOR) <- rowNames
colnames(FisherOR) <- colNames

# Labels by OR
tmp <- as.data.frame(FisherOR)
tmp$Rows <- rownames(tmp)
tmp <- melt(tmp)

tmp2 <- tmp %>% 
	group_by(Rows) %>%
    filter(value == max(value)) %>%
    arrange(Rows,variable,value) %>%
    as.data.frame()

write.table(tmp2,"output/Labels_Clusters.txt",sep="\t",quote=F)

#library(reshape2)
p <- FisherP %>% 
    as.data.frame() %>%
    rownames_to_column("Module") %>%
    melt() %>%
    group_by(Module) %>%
    mutate(FDR = p.adjust(value,method="BH")) %>%
    mutate(log = -log10(FDR)) %>%
    rename(Pval = value) %>%
    as.data.frame()
    
OR <- FisherOR %>%
    as.data.frame() %>%
    rownames_to_column("Module") %>%
    melt()%>%
    rename(OR = value)

p <- Reduce(dplyr::left_join, list(p, OR))
p$OR[!is.finite(p$OR)] <- max(p$OR[is.finite(p$OR)])
p$log[p$log < 1.3] <- NA
p$OR<- ifelse(is.na(p$log), p$log, p$OR)

#p$variable <- factor(p$variable,levels=c("Micro_L1-3","Endo_L2-6","Astro_L1-2", "Astro_L1-6","OPC_L1-6","Oligo_L1-6","Exc_L2","Exc_L2-3","Exc_L2-4","Exc_L3-4","Exc_L3-5","Exc_L4-5","Exc_L4-6","Exc_L5-6",
#                                          "Exc_L6","Inh_L1","Inh_L1-2","Inh_L1-3","Inh_L1-4","Inh_L2-3","Inh_L2-4", "Inh_L2-5","Inh_L2-6","Inh_L3-5","Inh_L3-6","Inh_L4-5","Inh_L4-6","Inh_L5-6"))

pdf("output/Allen_PFC_Enrichment.pdf",width=18,height=8)
ggscatter(p, 
        x = "variable",
        y = "Module",
        size="OR",
        color="log",
        alpha = 0.8,
        xlab = "",
            ylab = "",) +
        theme_minimal() +
        rotate_x_text(angle = 45)+
        #coord_flip()+
      scale_size(range = c(0.5, 10))+ 
      gradient_color(c("red","darkred"))
dev.off()
