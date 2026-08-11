
library(spacexr)
library(Matrix)
library(STdeconvolve)
library(ggplot2)
library(ggsci)
packageVersion("STdeconvolve")

setwd("./")
dir.create('Result')
savedir = './Result'

sc_matrix <- "../data/scRNA_expression_matrix.csv"
sc_barcode <- "../data/scRNA_coordinates.csv"
spatial_matrix <- "../data/spatial_expression_matrix.csv"
spatial_barcode <- "../data/spatial_coordinates.csv"

counts <- read.table(sc_matrix,sep="\t", header=T)
#rownames(counts) <- counts[,1]; counts[,1] <- NULL
meta_data <- read.table(sc_barcode,header=T,sep="\t")
cell_types <- meta_data$cluster; names(cell_types) <- meta_data$barcodes
cell_types <- as.factor(cell_types)
nUMI <- colSums(counts)
reference <- Reference(counts, cell_types, nUMI, min_UMI = 0)

counts1 <- read.table(spatial_matrix,sep="\t", header=T)
coords <- read.table(spatial_barcode,header=T,sep="\t")
#rownames(counts1) <- counts1[,1]; counts1[,1] <- NULL
rownames(coords) <- coords$barcodes; coords$barcodes <- NULL
nUMI <- colSums(counts1)
puck <- SpatialRNA(coords, counts1, nUMI)

myRCTD <- create.RCTD(puck, reference, max_cores = 1)
myRCTD <- run.RCTD(myRCTD, doublet_mode = 'doublet')
results <- myRCTD@results
chars <- capture.output(print(results))
writeLines(chars, con = file("DD_sc-spatial_bin60_imputation.result"))


## result
barcodes <- colnames(myRCTD@spatialRNA@counts)
weights <- myRCTD@results$weights
norm_weights <- normalize_weights(weights)

# plot Dentate weights
for(i in 0:(length(unique(meta_data$cluster))-1)){
	p <- plot_puck_continuous(myRCTD@spatialRNA, barcodes, norm_weights[,(i+1)], ylimit = c(0,0.5), 
                     title =paste0('plot of cluster_',i,' weights'), size=2, alpha=0.8) 
	ggsave(paste0(savedir, "/Spaital_weights_cluster_",i,".png"), width=8, height=6, plot=p,bg="white")

}
#STdeconvolve
m <- as.matrix(norm_weights)
p <- coords
colnames(p) <- c('x','y')

write.table(m,file = "output_norm_weights.txt",quote = F)
write.table(p,file = "output_coords.txt",quote = F)

plt <- vizAllTopics(theta = m,
             pos = p,
             topicOrder=seq(ncol(m)),
             topicCols=rainbow(ncol(m)),
             groups = NA,
             group_cols = NA,
             r = 0.5, # size of scatterpies; adjust depending on the coordinates of the pixels
             lwd = 0.3,
             showLegend = TRUE,
             plotTitle = "scatterpies")

## function returns a `ggplot2` object, so other aesthetics can be added on:
plt <- plt + ggplot2::guides(fill=ggplot2::guide_legend(ncol=2))
ggsave(paste0(savedir, "/Spaital_scatterpies.png"), width=12, height=6, plot=plt, bg="white")


