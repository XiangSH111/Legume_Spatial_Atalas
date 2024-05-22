library(dplyr)
library(tidyr)
library(Seurat)
library(stringr)
library(patchwork)
library(ggplot2)
library(cowplot)

###load rds data
adata <- readRDS('Medicago_nodule.RDS')

###plot
pdf("01_MX_umap.pdf", width = 20,height=20)
p1 <- DimPlot(object = adata, reduction = "umap", label = T,pt.size=1,group.by = "ident")
p2 <- DimPlot(object = adata, reduction = "umap", label = T,pt.size=1,group.by = "group")
p3 <- DimPlot(object = adata, reduction = "umap", pt.size=1,group.by = "group",split.by = "group")
p4 <- SpatialDimPlot(adata, stroke = 0, pt.size=8, label=T,label.size=2)
print(p1)
print(p2)
print(p3)
print(p4)
dev.off()



