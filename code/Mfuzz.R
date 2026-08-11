rm(list = ls())
library(Seurat)
library(Mfuzz)

protein <- read.delim("/data/MX_anno_exp.txt", row.names = 1, check.names = FALSE)

protein <- as.matrix(protein)
mfuzz_class <- new('ExpressionSet', exprs = protein)
mfuzz_class <- filter.NA(mfuzz_class, thres = 0.25)
mfuzz_class <- fill.NA(mfuzz_class, mode = 'mean')
mfuzz_class <- filter.std(mfuzz_class, min.std = 0)
mfuzz_class <- standardise(mfuzz_class)

set.seed(123)
cluster_num <- 18
sum(is.na(exprs(mfuzz_class)))
mfuzz_class <- mfuzz_class[-na_rows, ]
mfuzz_cluster <- mfuzz(mfuzz_class, c = cluster_num, m = mestimate(mfuzz_class))
mfuzz.plot2(mfuzz_class, cl = mfuzz_cluster, mfrow = c(3,6), time.labels = colnames(protein))

pdf("/output/MX_anno_exp_Mfuzz18.pdf", width = 20, height = 5, family = "GB1")
par(mar = c(0.5, 0.5, 0.2, 0.2), oma = c(0,0,0,0))

p1 <- mfuzz.plot2(
  mfuzz_class,
  cl = mfuzz_cluster,
  mfrow = c(3,6),
  time.labels = abbreviate(colnames(protein), 3),
  col.lab = "gray30",
  cex.lab = 0.8,
  x11 = FALSE
)
print(p1)
dev.off()

cluster_size <- mfuzz_cluster$size
names(cluster_size) <- 1:cluster_num
cluster_size

protein_cluster <- mfuzz_cluster$cluster
protein_cluster <- cbind(protein[names(protein_cluster), ], protein_cluster)
head(protein_cluster)
write.table(protein_cluster, "/output/MX_anno_exp_gene_list_cluster18.txt", sep = '\t', col.names = NA, quote = FALSE)

protein_cluster <- mfuzz_cluster$cluster
protein_standard <- mfuzz_class@assayData$exprs
protein_standard_cluster <- cbind(protein_standard[names(protein_cluster), ], protein_cluster)
head(protein_standard_cluster)
write.table(protein_standard_cluster, "/output/MX_anno_exp_gene_list_standard_cluster18.txt", sep = '\t', col.names = NA, quote = FALSE)

