
library(getopt)

arg <- matrix(c("input", "i","1","character","input file1",
                "minC","c","1","integer","genes expressed by at least this number of cells",
                "minG","g","1","integer","cells expressing at least this number of genes",
                "dims","d","1","integer","dims option for FindNeighbors,default=15",
                "cluster","l","2","character","annotation csv",
                "highgenes","h","2","integer","threshold to determine high expressed genes",
                "knum","k","2","integer","number of expected spatial clusters",
                "help","e","0","logical", "Usage: Rscript runiDrop.R -i <input> [-s SAMPLE -t TISSUE]"
                ),byrow=T,ncol=5)
opt = getopt(arg)
if(!is.null(opt$help) || is.null(opt$input)){
        cat(paste(getopt(arg, usage = T), "\n"))
 q()
}
if (is.null(opt$minC)){
        opt$minC <- 3
}
if (is.null(opt$minG)){
        opt$minG <- 100
}
if (is.null(opt$knum)){
    opt$knum<-8
}

library(Seurat)
library(remotes)
library(ggplot2)
library(cowplot)
library(dplyr)
library(Giotto)
library(stringr)

data<-read.csv(opt$input,sep=",")
rownames(data)<-data[,1]
data<-data[,-1]

locations<-data.frame(substring(names(data),2))
locations[,1:2]<-str_split_fixed(locations$substring.names.data...2.,"_",2)
names(locations)<-c("V1","V2")

locations$V1<-as.numeric(locations$V1)
locations$V2<-as.numeric(locations$V2)
h <- as.numeric((max(locations$V2) - min(locations$V2) + 1)/10)
w <- as.numeric((max(locations$V1) - min(locations$V1) + 1)/10)
my_giotto_object = createGiottoObject( raw_exprs = data,
                                      spatial_locs = locations,
                                      instructions = my_instructions
                                     )

#processing data
my_giotto_object <- filterGiotto(gobject = my_giotto_object,
                            expression_threshold = 1,
                            gene_det_in_min_cells = opt$minC,
                            #gene_det_in_min_cells = 3,     
                            min_det_genes_per_cell = opt$minG,
                            #min_det_genes_per_cell =100,
                            expression_values = c('raw'),
                            verbose = T)
my_giotto_object <- normalizeGiotto(gobject = my_giotto_object,verbose = T)
my_giotto_object <- addStatistics(gobject = my_giotto_object)
my_giotto_object <- adjustGiottoMatrix(gobject = my_giotto_object, 
                                       expression_values = c('normalized'),
                                       covariate_columns = c('nr_genes', 'total_expr'))


# run PCA
my_giotto_object <- runPCA(gobject = my_giotto_object)

my_giotto_object <- createNearestNetwork(gobject = my_giotto_object, dimensions_to_use = 1:opt$dims, k = 30)

#空间网络图
# get information about the Delaunay network
my_giotto_object = createSpatialNetwork(gobject = my_giotto_object, minimum_k=2,k = 8)
km_spatialgenes = binSpect(my_giotto_object,bin_method="rank")

km_genes=km_spatialgenes[km_spatialgenes$adj.p.value<=quantile(km_spatialgenes$adj.p.value,0.1),]$genes
spat_cor_netw_DT = detectSpatialCorGenes(my_giotto_object,
                                         method = 'network', spatial_network_name = 'Delaunay_network',
                                         subset_genes = km_genes)
#write.csv(x=spat_cor_netw_DT$cor_DT,"cor_qu0.1.csv")
top_cor_genes=spat_cor_netw_DT$cor_DT[spat_cor_netw_DT$cor_DT$spat_cor>0.8&spat_cor_netw_DT$cor_DT$spat_cor<1,]
#pdf("cor_qu0.1.pdf",width=20,height=20)
spat_cor_netw_DT = clusterSpatialCorGenes(spat_cor_netw_DT, name = 'spat_netw_clus', k = opt$knum)
#heatmSpatialCorGenes(my_giotto_object, spatCorObject = spat_cor_netw_DT, use_clus_name = 'spat_netw_clus')
#dev.off()

top_netw_spat_cluster = showSpatialCorGenes(spat_cor_netw_DT, use_clus_name = 'spat_netw_clus',
                                            selected_clusters = 1:opt$knum, show_top_genes = 1)

temp2<-merge(top_netw_spat_cluster[,c("gene_ID","clus")],top_cor_genes,by="gene_ID")
#write.csv(x=top_netw_spat_cluster,"top_netw_spat_cluster.csv")
#write.csv(x=temp2,"top_netw_spat_genes.csv")
netw_ranks = rankSpatialCorGroups(my_giotto_object, 
                                  spatCorObject = spat_cor_netw_DT, 
                                  use_clus_name = 'spat_netw_clus')
cluster_genes_DT = showSpatialCorGenes(spat_cor_netw_DT, use_clus_name = 'spat_netw_clus', show_top_genes = 1)
cluster_genes = cluster_genes_DT$clus; names(cluster_genes) = cluster_genes_DT$gene_ID
my_giotto_object = createMetagenes(my_giotto_object, gene_clusters = cluster_genes, name = 'cluster_metagene')
pdf("spatial_genes.pdf",width=w*3+2,height=h*6+11)
spatCellPlot(my_giotto_object,
             spat_enr_names = 'cluster_metagene',
             cell_annotation_values = netw_ranks$clusters,
             point_size = 3, cow_n_col = 2,
             point_shape = 'no_border',cell_color_gradient=c('blue','#FFFF99','red'))
dev.off()
