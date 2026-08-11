library(dplyr)
library(tidyr)
library(Seurat)
library(stringr)
library(patchwork)
library(ggplot2)
library(cowplot)
library(gridExtra)

####input data
my_rds <-readRDS("Medicago_nodule.RDS")
my_rds$organ <- paste0(my_rds$phase,'_',Idents(my_rds))

###input gene list
MX_rootmarker <- read.table('../data/MX_Interacting_gene.txt',header = F, sep = '\t',quote="")
MX_rootID <- gsub(" ", "", MX_rootmarker[,1], fixed = TRUE)

###AverageExpression
expr1 <- AverageExpression(my_rds,features = unique(gsub("_","-",MX_rootID)),group.by = "organ",slot = 'data')[['SCT']]
expr2 <- AverageExpression(my_rds,features = unique(gsub("_","-",MX_rootID)),slot = 'data')[['SCT']]
expr3 <- AverageExpression(my_rds,features = unique(gsub("_","-",MX_rootID)),group.by = "phase",slot = 'data')[['SCT']]

write.table(expr1,'MX_Interacting_gene_organ_phase_exp(rm_Senescence).txt',sep="\t",quote=FALSE,row.names=TRUE)
write.table(expr2,'MX_Interacting_gene_organ_exp(rm_Senescence).txt',sep="\t",quote=FALSE,row.names=TRUE)
write.table(expr3,'MX_Interacting_gene_phase_exp(rm_Senescence).txt',sep="\t",quote=FALSE,row.names=TRUE)

