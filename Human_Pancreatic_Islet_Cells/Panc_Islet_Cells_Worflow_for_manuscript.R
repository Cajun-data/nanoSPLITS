library(tidyverse)
#library(MSnSet.utils)
library(circlize)
library(DreamAI)
library(impute)
library(mclust)
library(ComplexHeatmap)

rotate <- function(x) t(apply(x, 2, rev))


###https://www.rdocumentation.org/packages/graph4lg/versions/0.3.0/topics/reorder_mat
reorder_mat <- function(mat, order){
  n <- length(order)
  if(!inherits(mat, "matrix")){
    stop("'mat' must be a matrix")
  } else if (!inherits(order, "character")){
    stop("'order' must be a character vector")
  } else if(!(isSymmetric(mat))){
    stop("The matrix 'mat' must be symmetric")
  } else if (n != length(colnames(mat))){
    stop("'order' must have as many elements as there are rows and
         columns in 'mat'")
  } else if(length(which(colnames(mat) %in% order)) != n){
    stop("The column names of the matrix you want to reorder must
         be present in the vector 'order'")
  } else if (length(which(row.names(mat) %in% order)) != n){
    print("The row names of the matrix you want to reorder must
          be present in the vector 'order'")
  } else {
    
    mat2 <- mat[order, order]
    
    return(mat2)
  }
}


####Impute and summarize to protein-level function

Summarize_pept <- function(data, Sample_groups = TRUE, group_patterns,
                             meta_cols = 1) {
  
  meta <- data[,1:meta_cols] 
  meta <- meta %>%
    mutate(Protein = case_when(Gene == "" ~ Protein.ID,
                               TRUE ~ Gene))
  dfList <- list()
  
  
  {
    for (group in group_patterns) {
      
      dat1 <- data %>%
        dplyr::select(-contains(colnames(meta))) %>%
        {if(Sample_groups) dplyr::select(.,contains(group)) else .} %>%
        mutate(Var1 = rowSums(is.na(.))) %>%
        filter(Var1 <= max(Var1)*0.6) %>%  ### Filter 40% missing values
        dplyr::select(-Var1) %>%
        as.matrix() %>%
        DreamAI::impute.KNN(.,k = 10) %>%  ### Impute with k of 10
        {as.data.frame(.) ->> dat2 } %>%
        as.data.frame() %>%
        rownames_to_column(var = "Modified.Sequence") %>%
        pivot_longer(-Modified.Sequence,names_to = "SampleID",
                     values_to = "Intensity") %>%
        inner_join(.,meta) %>%
        distinct(Protein, Modified.Sequence, SampleID, Intensity) %>%
        dplyr::rename(protein_list = Protein,
                      sample_list = SampleID,
                      id = Modified.Sequence,
                      quant = Intensity) %>%
        as.list() %>%
        iq::fast_MaxLFQ(.) %>%
        .[[1]]
      
      dfList[[group]] <- dat1
      
    }
    return(dfList)
  }
  
}


###Modified fast_correlation function

fast_cor <- function(data, use = "pairwise.complete.obs", method = "spearman") 
{
  {
    n = NULL
    n <- psych::pairwiseCount(data) - 2 ## Create matrix of degrees freedom
  }
  if (method == "spearman") {
    data <- base::apply(X = data, MARGIN = 2, data.table::frankv)
  }
  r = NULL
  p1 = NULL
  p2 = NULL
  r <- coop::pcor(x = data, use = use)
  {
    
    t <- sqrt(n)*r/sqrt(1 - r^2) ## T-test
    p1 <- stats::pt(t, n)  ## One side p-value
    p2 <- stats::pt(t, n, lower.tail = FALSE) ## Unfortunately run twice...
  }
  flt.Corr.Matrix <- function(cormat, pmat1 = NULL, pmat2 = NULL,
                              df = NULL) {
    ut <- base::upper.tri(cormat)
    flt_data <- data.frame(Var1 = base::rownames(cormat)[base::row(cormat)[ut]], 
                           Var2 = base::rownames(cormat)[base::col(cormat)[ut]], 
                           cor = cormat[ut])
    if (!is.null(pmat1)) 
      flt_data$p1 <- pmat1[ut]
    if (!is.null(pmat2)) 
      flt_data$p2 <- pmat2[ut]
    if (!is.null(df)) 
      flt_data$df <- df[ut]
    return(flt_data)
  }
  {
    result <- flt.Corr.Matrix(cormat = r,
                              df = n,
                              pmat1 = p1,
                              pmat2 = p2)
    result <-  base::transform(result, p = base::pmin(p1, p2)*2)
    result$FDR <- stats::p.adjust(result$p, method = "BH")
    result <- base::subset(result, select = -c(p1,p2))
    result <- list(result, r)
  }
  return(result)
}



###################################################




panc_pept <- read.delim("combined_modified_peptide_MaxLFQ.tsv") %>%
  filter(!grepl("contam", Protein)) %>%
  filter(grepl("HUMAN", Protein)) %>%
  filter(!grepl("Keratin", Protein.Description)) 
SampID <- data.frame(SampleID = colnames(panc_pept)) %>%
  mutate(SampleID = gsub("X", "" , SampleID))
annot <- read.csv("LDA_Celltypes.csv",
                  check.names = T)

annot <- annot %>%
  mutate(SampleID2 = paste(SampleID, Celltype, sep = "_")) %>%
  inner_join(., SampID)


meta <- panc_pept[c(2,9,10,11,13)]

panc_pept <- panc_pept %>%
  column_to_rownames(var = "Modified.Sequence") %>%
  dplyr::select(contains(annot$SampleID)) %>%
  log2() %>%
  rownames_to_column(var  = "Modified.Sequence") %>%
  pivot_longer(-`Modified.Sequence`,
               names_to = "SampleID",
               values_to = "Intensity") %>%
  filter(Intensity != -Inf) %>%
  mutate(SampleID = gsub("X","", SampleID)) %>%
  full_join(.,annot) %>%
  group_by(SampleID) %>%
  add_count(name = "TotalPeptideIDs") %>%
  filter(TotalPeptideIDs > 3000) %>% ### remove cells with fewer than 3000 peptides
  ungroup() %>%
  inner_join(.,meta) %>%
  dplyr::select(SampleID2,Assigned.Modifications, Modified.Sequence, Protein.ID, Gene, Intensity) %>%
  pivot_wider(names_from = "SampleID2",
              values_from = "Intensity") 

batch_meta <- data.frame(SampleID = colnames(panc_pept),
                         Batch = NA)
batch_meta <- batch_meta %>%
  filter(grepl("2", SampleID)) %>%
  mutate(Batch = case_when( grepl("22", SampleID) ~ "Donor 1",
                            grepl("23", SampleID)~ "Donor 2"
  )) 


panc_pept_batch <- panc_pept %>%
  column_to_rownames(var = "Modified.Sequence") %>%
  dplyr::select(-contains(colnames(meta))) %>%
  MSnSet.utils::ComBat.NA(., batch_meta$Batch) 

panc_pept_batch2 <- panc_pept_batch$`corrected data` %>%
  as.data.frame() %>%
  dplyr::select(contains(c("alpha","beta","delta"))) %>%
  rownames_to_column(var  = "Modified.Sequence") %>%
  inner_join(.,meta) %>%
  relocate(colnames(meta)) %>%
  dplyr::select(-Protein) %>%
  mutate(Modified.Sequence2 = Modified.Sequence) %>%
  column_to_rownames(var = "Modified.Sequence2")

  



bg <- unique(panc_pept$Gene)


group_patterns <- c("alpha","beta", "delta")


final_table <- Summarize_pept(panc_pept_batch2, Sample_groups = T,
                                group_patterns = group_patterns,
                                meta_cols = 4)

###################################
##########################

#######Alpha cells
alpha <- final_table[[1]]
 
alpha_cor <- fast_cor(t(alpha), method = "pearson")



mclust.options(subset = 4000)
alpha_clust <- Mclust(t(scale(t(alpha))),  G = 15,
                      modelNames = "VII")
clust_order <- data.frame(
  Gene = names(alpha_clust$classification),
  Cluster = alpha_clust$classification,
  uncertainty = alpha_clust$uncertainty) %>%
  arrange(desc(factor(Cluster)),
          uncertainty) %>%
  filter(Gene %in% rownames(alpha))


alpha_cor_mat <- alpha_cor[[2]] 

include_list <- clust_order$Gene

alpha_cor_mat <- alpha_cor_mat[include_list,include_list]


alpha_cor_mat <- rotate(reorder_mat(alpha_cor_mat,
                                      order = clust_order$Gene))


#png("alpha_heatmap_forpublication.png", width = 2100, height = 2100)
Heatmap(alpha_cor_mat, cluster_rows = F,
                cluster_columns = F,
                show_row_names = F,
                show_column_names = F,
                use_raster = F,
                show_heatmap_legend = T,
                heatmap_legend_param = list(title = "Pearson R")
)
#dev.off()



##########################
##########################
beta <- final_table[[2]]

beta_cor <- fast_cor(t(beta), method = "pearson")



mclust.options(subset = 4000)
beta_clust <- Mclust(t(scale(t(beta))),  G = 15,
                      modelNames = "VII")
clust_order2 <- data.frame(
  Gene = names(beta_clust$classification),
  Cluster = beta_clust$classification,
  uncertainty = beta_clust$uncertainty) %>%
  arrange(desc(factor(Cluster)),
          uncertainty) %>%
  filter(Gene %in% rownames(beta))


beta_cor_mat <- beta_cor[[2]] 

include_list <- clust_order2$Gene

beta_cor_mat <- beta_cor_mat[include_list,include_list]


beta_cor_mat <- rotate(reorder_mat(beta_cor_mat,
                                    order = clust_order2$Gene))


#png("beta_heatmap_forpublication.png", width = 2100, height = 2100)
Heatmap(beta_cor_mat, cluster_rows = F,
        cluster_columns = F,
        show_row_names = F,
        show_column_names = F,
        use_raster = F,
        show_heatmap_legend = T,
        heatmap_legend_param = list(title = "Pearson R")
)
#dev.off()


########################################################
########################################################
########################################################



alpha <- final_table[[1]]
beta <- final_table[[2]]


alpha_long <- alpha %>%
  as.data.frame() %>%
  rownames_to_column(var = "Gene") %>%
  pivot_longer(-Gene, values_to = "Intensity",
               names_to = "SampleID") 
beta_long <- beta %>%
  as.data.frame() %>%
  rownames_to_column(var = "Gene") %>%
  pivot_longer(-Gene, values_to = "Intensity",
               names_to = "SampleID") 


##### Cluster Comparison

beta_cell <- clust_order2 %>%
  filter(Cluster == 15) %>%
  mutate(Type = "Beta")

alpha_cell <- clust_order %>%
  filter(Cluster == 8) %>%
  mutate(Type = "Alpha")

beta_long_check <- beta_long %>%
  distinct(Gene)

alpha_long_check <- alpha_long %>%
  distinct(Gene)

both_cell <- beta_cell %>%
  full_join(., alpha_cell) %>%
  group_by(Gene) %>%
  filter(Gene != "GCG") %>%
  filter(Gene != "PTPRN2") %>%
  add_count(name = "n") %>%
  mutate(Type = case_when(n == 2 ~ "Both",
                          Type == "Alpha" & n == 1 ~ "Alpha",
                          Type == "Beta" & n == 1 ~ "Beta")) %>%
  distinct(Gene,Type) %>%
  ungroup() %>%
  arrange(Type = fct_relevel(Type, "Both", "Beta", "Alpha"))



correlations_GCG_INS <- beta_cor[[2]]

include_list <- both_cell$Gene

correlations_GCG_INS <- correlations_GCG_INS[include_list,include_list]


correlations_GCG_INS <- reorder_mat(correlations_GCG_INS,
                                      order = both_cell$Gene)




col_fun = colorRamp2(c(-0.5,0,1), c( "blue","white", "red"))

#png(filename = "beta_cells_secretory_proteins.png", width = 6,
#    height = 4,
#    units = "in",
 #   res = 300)
Heatmap(correlations_GCG_INS, cluster_rows = T,
        clustering_method_rows = "ward.D",
        clustering_method_columns = "ward.D",
        cluster_columns = T,
        col = col_fun,
        show_row_names = T,
        show_column_names = T,
        use_raster = F,
        show_heatmap_legend = F,
        row_names_centered = T,
        column_names_centered = T,
        row_names_gp = gpar(col = c(rep("#000000", 16),
                                    rep("#00BA38", 11),
                                    rep("#93AA00", 5)),
                            fontsize = 9),
        column_names_gp = gpar(col = c(rep("#000000", 16),
                                       rep("#00BA38", 11),
                                       rep("#93AA00", 5)),
                               fontsize = 9),
        heatmap_legend_param = list(title = "Pearson R")
)
#dev.off()





alpha_long_check <- alpha_long %>%
  distinct(Gene)

both_cell <- beta_cell %>%
  full_join(., alpha_cell) %>%
  group_by(Gene) %>%
  filter(Gene != "INS") %>%
  add_count(name = "n") %>%
  mutate(Type = case_when(n == 2 ~ "Both",
                          Type == "Alpha" & n == 1 ~ "Alpha",
                          Type == "Beta" & n == 1 ~ "Beta")) %>%
  distinct(Gene,Type) %>%
  inner_join(., alpha_long_check) %>%
  arrange(Type = fct_relevel(Type, "Both", "Beta", "Alpha"))




correlations_GCG_INS <- alpha_cor[[2]]

include_list <- both_cell$Gene

correlations_GCG_INS <- correlations_GCG_INS[include_list,include_list]


correlations_GCG_INS <- reorder_mat(correlations_GCG_INS,
                                           order = both_cell$Gene)





col_fun = colorRamp2(c(-0.5,0,1), c( "blue","white", "red"))

#png(filename = "alpha_cells_secretory_proteins.png", width = 6,
#    height = 4,
#    units = "in",
#    res = 300)
Heatmap(correlations_GCG_INS, cluster_rows = T,
        clustering_method_rows = "ward.D2",
        clustering_method_columns = "ward.D2",
        cluster_columns = T,
        col = col_fun,
        show_row_names = T,
        show_column_names = T,
        use_raster = F,
        show_heatmap_legend = F,
        row_names_centered = T,
        column_names_centered = T,
        row_names_gp = gpar(col = c(rep("#000000", 16),
                                    rep("#00BA38", 9),
                                    rep("#93AA00", 7)),
                            fontsize = 9),
        column_names_gp = gpar(col = c(rep("#000000", 16),
                                       rep("#00BA38", 9),
                                       rep("#93AA00", 7)),
                               fontsize = 9),
        heatmap_legend_param = list(title = "Pearson R")
)
#dev.off()

 