library(dplyr)
library(stringr)
library(ggplot2)
library(Seurat)
library(monocle3)

#
my_rds <- readRDS(Medicago_part_nodule.RDS)
my_rds$organ <- Idents(my_rds)

#load data
data <- as(as.matrix(my_rds@assays$Spatial@counts), 'sparseMatrix')
pd <- data.frame(my_rds@meta.data,Idents(my_rds))
fData <- data.frame(gene_short_name = row.names(data), row.names = row.names(data))

cds <- new_cell_data_set(data,
                         cell_metadata  = pd,
                         gene_metadata  = fData)
#remove batch effect
cds <- preprocess_cds(cds, num_dim = 100)
#cds <- align_cds(cds, alignment_group = "group")
cds <- align_cds(cds)

#reduce dimensions
cds <- reduce_dimension(cds)

#determine how many dimensions to use
pdf("00_dim_s1.pdf")
plot_pc_variance_explained(cds)
dev.off()

pdf("01_UMAP_s1.pdf", width = 7, height = 7)
plot_cells(cds, label_groups_by_cluster=TRUE,  color_cells_by = "Idents.my_rds.",group_label_size = 5,cell_size = 1)
dev.off()

#clustering
cds <- cluster_cells(cds)
pdf("02_partition_s1.pdf", width = 7, height = 7)
plot_cells(cds, label_groups_by_cluster=TRUE,  color_cells_by = "partition")
dev.off()

#learn graph
cds <- learn_graph(cds)
pdf("03_trajectorise_s1.pdf", width = 7, height = 7)
plot_cells(cds,
           color_cells_by = "Idents.my_rds.",
           label_groups_by_cluster=FALSE,
           group_label_size = 4,
           label_leaves=TRUE,
           label_branch_points=TRUE)
plot_cells(cds,
           color_cells_by = "phase",
           label_groups_by_cluster=FALSE,
           label_cell_groups=TRUE,
           group_label_size = 4,
           label_leaves=TRUE,
           label_branch_points=TRUE)
plot_cells(cds,
           color_cells_by = "group",
           label_groups_by_cluster=FALSE,
           label_cell_groups=TRUE,
           group_label_size = 4,
           label_leaves=TRUE,
           label_branch_points=TRUE)
plot_cells(cds,
           color_cells_by = "organ",
           label_groups_by_cluster=FALSE,
           label_cell_groups=TRUE,
           group_label_size = 4,
           label_leaves=TRUE,
           label_branch_points=TRUE)
dev.off()

get_earliest_principal_node <- function(cds, time_bin= "Meristematic zone"){
	cell_ids <- which(colData(cds)[,"Idents.my_rds."] == time_bin)
	closest_vertex <-
	cds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
	closest_vertex <- as.matrix(closest_vertex[colnames(cds), ])
	root_pr_nodes <-
	igraph::V(principal_graph(cds)[["UMAP"]])$name[as.numeric(names
	(which.max(table(closest_vertex[cell_ids,]))))]
	root_pr_nodes
}

cds = order_cells(cds, root_pr_nodes=get_earliest_principal_node(cds))

pdf("trajectorise_set_nodule_s1.pdf", width = 7, height = 7)
plot_cells(cds,
           color_cells_by = "pseudotime",
           label_cell_groups=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           graph_label_size=1.5)
dev.off()

pseudotime <- pseudotime(cds, reduction_method = 'UMAP')
coordinate <- str_split_fixed(names(pseudotime),"X",2)[,2] %>% str_split_fixed("_",2)
a <- coordinate[,1]
b <- coordinate[,2]
a<-as.numeric(a)
b<-as.numeric(b)

c<-data.frame(a,b,pseudotime)
names(c)<-c("x","y","pseudotime")

h <- as.numeric(max(c$y) - min(c$y) + 1)
w <- as.numeric(max(c$x) - min(c$x) + 1)
pdf("spatial_plot.pdf", width = w/15, height = h/15)
ggplot(c) + geom_point(mapping=aes(x=x, y=y, colour = pseudotime),size = 1.5)+
  scale_color_gradient(low = "blue",high = "#eaed18")
dev.off()

#Finding genes that change as a function of pseudotime
ciliated_cds_pr_test_res <- graph_test(cds, neighbor_graph="principal_graph", cores=4)
write.csv(subset(ciliated_cds_pr_test_res, q_value < 0.05),"cds_pr_test_res.csv")
pr_deg_ids <- row.names(subset(ciliated_cds_pr_test_res, q_value < 0.05))

pdf("cds_pr_test_res.pdf",width=10,height=5)
plot_cells(cds, genes=head(pr_deg_ids),
           show_trajectory_graph=FALSE,
           label_cell_groups=FALSE,
           label_leaves=FALSE)
dev.off()

gene_module_df <- find_gene_modules(cds[pr_deg_ids,], resolution=c(10^seq(-6,-1)))
cell_group_df <- tibble::tibble(cell=row.names(colData(cds)),  cell_group=colData(cds)$Idents.my_rds.)
agg_mat <- aggregate_gene_expression(cds, gene_module_df, cell_group_df)
row.names(agg_mat) <- stringr::str_c("Module ", row.names(agg_mat))
pdf("gene_module.pdf")
pheatmap::pheatmap(agg_mat,
                   scale="column", clustering_method="ward.D2")
dev.off()

write.table(gene_module_df, "module_gene.txt",col.names = T, append = F,sep = '\t',quote = F)

