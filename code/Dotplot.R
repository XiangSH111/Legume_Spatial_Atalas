library(dplyr)
library(tidyr)
library(Seurat)
library(stringr)
library(patchwork)
library(ggplot2)
library(cowplot)
library(gridExtra)

my_rds <-readRDS("/data/work/05.MX/MX_rds_normalized_SCT.rds")
DefaultAssay(my_rds) <- "SCT"
my_gene_list <- read.table('/data/work/05.MX/11.hormone_nodule_dotplot/01.MX_hormone/1_9.plot_gene.txt',header = T,sep = '\t')

new_list_ID <- my_gene_list[,4]
new_list_name <- my_gene_list[,2]

unique(Idents(my_rds))
my_rds$condition2 <- my_rds$phase_organ
my_rds <- subset(my_rds,subset = condition2 == "P5_Meristematic zone" | condition2 == "P4_Meristematic zone" | condition2 == "P3_Meristematic zone" | condition2 == "P2_Meristematic zone" | condition2 == "P1_Meristematic zone"  | condition2 == "P5_Proximal infection zone" | condition2 == "P4_Proximal infection zone" | condition2 == "P3_Proximal infection zone" | condition2 == "P5_Distal infection zone" | condition2 == "P4_Distal infection zone" | condition2 == "P3_Distal infection zone")
my_rds$condition2 <- factor(my_rds$condition2, levels = c("P5_Proximal infection zone","P4_Proximal infection zone","P3_Proximal infection zone","P5_Distal infection zone","P4_Distal infection zone","P3_Distal infection zone","P5_Meristematic zone","P4_Meristematic zone","P3_Meristematic zone","P2_Meristematic zone","P1_Meristematic zone"))


pdf("MX_hormone_dotplot_ID_SCT.pdf",width=15,height=6)
p1<-DotPlot(my_rds, group.by = "condition2", cols=c("blue","red"), features = unique(gsub("_","-",new_list_ID)))+theme(axis.text.x = element_text(angle = 90, hjust = 1))
#ggsave2("1_9_select_plot_gene.png",p1,dpi=300,width=15,height=15)
p1
dev.off()

pdf("MX_hormone_dotplot_name_SCT.pdf",width=15,height=6)
p2<-DotPlot(my_rds, group.by = "condition2", cols=c("blue","red"), features = unique(gsub("_","-",new_list_ID)))+theme(axis.text.x = element_text(angle = 90, hjust = 1))+scale_x_discrete(labels=new_list_name)
#ggsave2("1_9_select_plot_gene.png",p1,dpi=300,width=15,height=15)
p2
dev.off()

