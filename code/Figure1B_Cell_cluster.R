###seurat using spatial data
library(dplyr)
library(tidyr)
library(Seurat)
library(stringr)
library(patchwork)
library(ggplot2)
library(cowplot)
dir.create("./gem_2_rds.output")
setwd("./gem_2_rds.output")
nodule.data <- read.csv("BIN60FP200000343BR_C4_5_P4.csv", header = T, sep=",")
rownames(nodule.data) <- nodule.data[,1]
nodule.data <- nodule.data[,2:dim(nodule.data)[2]]
meta <- data.frame(Barcodes = colnames(nodule.data), Sample = "Medi", Tissue = "nodule", Annotation = NA, Celltype = NA)
adata <- CreateSeuratObject(counts = nodule.data,
                              project = "nodule",
                              min.cells = 0,
                              min.features = 0,
                              assay = 'Spatial')

###adding location information
locations<-data.frame(substring(colnames(nodule.data),2))
locations[,1:2]<-str_split_fixed(locations$substring.colnames.nodule.data...2.,"_",2)
names(locations)<-c("x","y")
rownames(x = locations) <-colnames(nodule.data)
locations$x<-as.numeric(locations$x)
locations$y<-as.numeric(locations$y)
adata[['image']] <- new(Class = 'SlideSeq',assay = "Spatial",
    coordinates = locations)

pdf("QCVlnplot.pdf", width = 7, height = 3)
plot1 <- VlnPlot(adata, features = 'nCount_Spatial', pt.size = 0, log = TRUE) + NoLegend()
adata$log_nCount_Spatial <- log(adata$nCount_Spatial)
plot2 <- SpatialFeaturePlot(adata, features = 'log_nCount_Spatial', pt.size.factor = 2, stroke = 0.1) + theme(legend.position = "right")
wrap_plots(plot1, plot2)
dev.off()
adata <- subset(x = adata, subset = nFeature_Spatial >= 10)

###SCTransform to normalize the data
adata <- SCTransform(adata, assay = "Spatial", verbose = FALSE)
adata <- RunPCA(adata,assay="SCT")

pdf("ElbowPlot.pdf", width = 10)
ElbowPlot(object = adata)
dev.off()

adata <- FindNeighbors(adata, dims = 1:20)
adata <- FindClusters(adata, resolution = 0.5, verbose = FALSE)
adata <- RunUMAP(adata, dims = 1:20)

pdf("umap.pdf", width = 16)
plot1 <- DimPlot(object = adata, reduction = "umap", label = T,pt.size=1)
plot2 <- SpatialDimPlot(adata, stroke = 0, pt.size=2, label=T)
plot1+plot2
dev.off()

DefaultAssay(adata) <- "SCT"
adata <- FindSpatiallyVariableFeatures(adata, assay = "SCT", slot = "scale.data", features = VariableFeatures(adata)[1:2000],
    selection.method = "moransi", x.cuts = 100, y.cuts = 100)

###plot variable features with and without labels
pdf("VariableFeatures.pdf", height = 15, width = 20)
plot1 <- SpatialFeaturePlot(adata, features = head(SpatiallyVariableFeatures(adata, selection.method = "moransi"),10), ncol = 5, alpha = c(0.1, 1), max.cutoff = "q95",pt.size.factor = 2)
plot1
dev.off()

pdf("cell_number.pdf", width = 12, height = 8)
library(ggplot2)
dat<-as.data.frame(table(adata@active.ident))
p<-ggplot(dat,aes(Var1,Freq,fill=Var1))+geom_bar(stat = "identity",position
= "dodge")+geom_text(aes(label=Freq,vjust=-0.5)
)
p
dev.off()

sample.cluster<-AverageExpression(adata)$SCT
hc = hclust(dist(t(sample.cluster)))
hcd = as.dendrogram(hc)
pdf("hclust_integrated.pdf")
plot(hcd)
dev.off()

DefaultAssay(object = adata) <- "SCT"
adata.markers <- FindAllMarkers(object = adata,
                                 only.pos = T,
                                 min.pct = 0.5,
                                 logfc.threshold = 0.5)
write.table(adata.markers,file="markers_0.5.txt",sep="\t",quote=F)

adata.markers.spatial <- FindAllMarkers(object = adata,
                                 assay="Spatial",
                                 only.pos = T,
                                 min.pct = 0.5,
                                 logfc.threshold = 0.5)
write.table(adata.markers.spatial,file="markers_0.5.spatial.txt",sep="\t",quote=F)

pdf ("markergeneheatmap.pdf")
top10 <- adata.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
DoHeatmap(adata, features = top10$gene) + NoLegend()
dev.off()

saveRDS(adata,file="data.RDS")

