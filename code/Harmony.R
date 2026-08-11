library(harmony)
library(tidyr)
library(Seurat)
library(stringr)
library(patchwork)
library(ggplot2)
library(cowplot)
library(readr)
library(gridExtra)
library(RColorBrewer)

rdsf <- read_tsv("../data/MX_RDS.info")
nodule<-list()
id<-c()
vf<-vector()

for (i in 1:length(rdsf$Rds)){
	nodule[[i]]<-readRDS(rdsf$Rds[i])
	nodule[[i]]@meta.data$group<-rdsf$Sample[i]
	vf <- c(vf,VariableFeatures(nodule[[i]]))
	id[i]<- rdsf$Sample[i]
}

nodule<-merge(nodule[[1]], y=nodule[-1], add.cell.ids = id, project="nodule")

VariableFeatures(nodule) <- vf
nodule <- ScaleData(object = nodule,verbose = FALSE,features=VariableFeatures(nodule))
nodule <- RunPCA(object = nodule, verbose = FALSE)
nodule <-RunHarmony(nodule,'group',plot_convergenc=TRUE,assay.use="Spatial")
nodule <- RunUMAP(nodule, reduction = "harmony", dims = 1:30)

nodule <- FindNeighbors(nodule, reduction = "harmony", dims = 1:30)
nodule <- FindClusters(nodule, resolution = 0.5)

pdf("umap.pdf", width = 25)
plot1 <- DimPlot(nodule, group.by = c("ident","group"))
plot1
dev.off()

mypalette<-c(brewer.pal(11,"BrBG"),brewer.pal(11,"PiYG"))

nodules<- SplitObject(nodule, split.by = "group")

p<-list()
for(i in 1:length(rdsf$Rds))
    {
    nodules[[i]]@images<- nodules[[i]]@images[i]
    p[[i]]<-SpatialDimPlot(nodules[[i]],stroke=0,pt.size=4,label=T,label.size=2)+scale_fill_manual(breaks=as.character(c(0:21)),values=mypalette)
}
sum<-grid.arrange(grobs=p,ncol=3)
ggsave2("Spatial_dim_plot.png",sum,dpi=600,width=20,height=20)
saveRDS(nodule,"data.RDS")

