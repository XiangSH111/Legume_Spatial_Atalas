library(dplyr)
library(tidyr)
library(Seurat)
library(stringr)
library(patchwork)
library(ggplot2)
library(cowplot)
library(gridExtra)

####input rds
my_rds <-readRDS("DD_nodule.RDS")
my_rds$organ <- paste0(my_rds$phase,'_',Idents(my_rds))

###input gene list
gene_list <- read.table('../data/Figure1_data/Figure1E/DD_Interacting_gene.txt',header = T, sep = '\t',quote="")
gene <- gsub(" ","",gene_list[,1])

###AverageExpression
expr1 <- AverageExpression(my_rds,features = unique(gsub("_","-",gene)),group.by = "organ",slot = 'data')[['Spatial']]
expr2 <- AverageExpression(my_rds,features = unique(gsub("_","-",gene)),slot = 'data')[['Spatial']]
expr3 <- AverageExpression(my_rds,features = unique(gsub("_","-",gene)),group.by = "phase",slot = 'data')[['Spatial']]

write.table(expr1,'DD_marker_gene_organ_phase_exp.txt',sep="\t",quote=FALSE,row.names=TRUE)
write.table(expr2,'DD_marker_gene_organ_exp.txt',sep="\t",quote=FALSE,row.names=TRUE)
write.table(expr3,'DD_marker_gene_phase_exp.txt',sep="\t",quote=FALSE,row.names=TRUE)

