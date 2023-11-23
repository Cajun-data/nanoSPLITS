library(tidyverse)
library(PCAtools)
library(SingleCellExperiment)
library(DEsingle)
library(edgeR)
library(BiocParallel)
library(ComplexHeatmap)
######################## PCA
##########################
RNA_norm3 <- read.csv("C10Treated_NormCounts.csv",
                      check.names = TRUE) %>%
  column_to_rownames(var = "X")

RNA_norm3[RNA_norm3 == 0] <- NA
RNA_norm3 <- na.omit(RNA_norm3)

pcatest2 <- pca(RNA_norm3, scale = TRUE)

namesc2 <- as.character(colnames(RNA_norm3)) %>%
  str_sub(1,4)


#pdf('PCA_RNAseq.pdf', width = 7, height = 7)
biplot(pcatest2, x = "PC1", y = "PC2",
       showLoadings = FALSE,
       ntopLoadings = 5,
       labSize = 0.01,
       drawConnectors = FALSE,
       lab = namesc2,
       encircle = TRUE,
       legendPosition = "right"
)+
  theme_minimal(base_size = 22)+
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(),
        text=element_text(family="Helvetica"))+
  scale_fill_discrete(name = "")
#dev.off()


#################################
###############################
RNA_norm3 <- read.csv("C10Treated_NormCounts.csv",
                      check.names = TRUE) %>%
  column_to_rownames(var = "X")

RNA_norm3$Percent_0s<-(rowSums(RNA_norm3==0)/ncol(RNA_norm3))*100

RNA_norm3 <- RNA_norm3 %>%
  filter(Percent_0s <= 50) %>%
  dplyr::select(-Percent_0s)


type <- RNA_norm3 %>%
  pivot_longer(everything(), names_to = "SampleID") %>%
  dplyr::select(SampleID) %>%
  distinct(SampleID) %>%
  mutate(cell_type1 = case_when(grepl("Ctrl", SampleID) ~ "Ctrl",
                                grepl("Trtd", SampleID) ~ "Trtd")) %>%
  remove_rownames() %>%
  column_to_rownames(var="SampleID")

type$cell_type1 <- as.factor(type$cell_type1)

sce2 <- SingleCellExperiment(
  assays = list(
    counts = as.matrix(RNA_norm3)
  )
)

results <- DEsingle(counts = sce2, group = type$cell_type1,
                    parallel = TRUE, BPPARAM = bpparam())

results2 <- results %>%
  mutate(log2FC = -log2(norm_foldChange)) %>%
  filter(total_mean_1 != 0) %>%
  filter(pvalue.adj.FDR != 0) %>%
  mutate(logpvalue = -log10(pvalue.adj.FDR)) %>%
  rownames_to_column(var = "Gene") %>%
  select(Gene, log2FC, logpvalue, pvalue.adj.FDR) 

results3 <- results2 %>%
  mutate(Significant = case_when(pvalue.adj.FDR < 0.01 & !between(log2FC, -0.5,0.5) ~ TRUE,
                                 TRUE ~ FALSE)) %>%
  filter(Significant == TRUE)


RNA_mat <- RNA_norm3 




RNA_mat<- cpm(RNA_mat, log = T)
RNA_mat <- as.data.frame(t(scale(t(RNA_mat)))) %>%
  filter(row.names(.) %in% results3$Gene) %>%
  as.matrix()


col_limit = 2.5

set.seed(10)
 RNA_heat <- Heatmap(RNA_mat, cluster_rows = T,
        clustering_distance_rows = "pearson",
        clustering_method_rows = "ward.D",
        clustering_distance_columns = "euclidean",
         clustering_method_columns = "single",
        cluster_row_slices = T,
        column_km = 2,
        top_annotation = HeatmapAnnotation(foo = anno_block(gp = gpar(fill = c(4,2)))),
        row_km = 4,
        show_row_names = FALSE,
        show_column_dend = F,
        show_row_dend = F,
        border = TRUE,
       use_raster = T,
        row_names_side = "left", column_names_side = "top", 
        col = circlize::colorRamp2(
          seq(-col_limit, col_limit,(col_limit/5)),
                                  rev(RColorBrewer::brewer.pal(11, "RdBu"
                                                                ))))


 
#pdf("RNA_heatmap.pdf", height = 10, width = 10)
RNA_heat
#dev.off()

#save(results2, file = "scRNAseq_DE.RData")


