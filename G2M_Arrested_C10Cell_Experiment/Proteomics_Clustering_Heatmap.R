library(tidyverse)
library(proDA)
library(DreamAI)
library(PCAtools)
library(DEP)
library(ComplexHeatmap)


C10protein <- read.delim("C10Treated_Protein_intensities.tsv") 


meta <- C10protein[,1:4] %>%
  filter(grepl("_MOUSE", Protein)) 


C10protein <- C10protein %>%
  dplyr::select(-Gene, -Protein, -Protein.ID ) %>%
  column_to_rownames(var = "Entry.Name")

C10protein[C10protein == 0] <- NA
C10protein[,1:96] <- log2(C10protein[,1:96]) 
C10protein <- as.data.frame(median_normalization(as.matrix(C10protein))) %>%
  rownames_to_column(var = "Protein")


C10protein_long <- pivot_longer(C10protein, -Protein, names_to = "SampleID",
                           values_to = "Intensity",
                           names_repair = "check_unique") %>%
  filter(!is.na(Intensity)) %>%
  filter(grepl("_MOUSE", Protein))

filterS<-C10protein_long %>%
  distinct(Protein, SampleID) %>%
  mutate(Group1 = case_when(grepl("Trtd", SampleID) ~ "Treated",
                           grepl("Ctrl", SampleID)~ "Control"
                           )) %>%
  group_by(SampleID) %>%
  add_count(name = "n") %>%
  distinct(SampleID, n, Group1) %>%
  filter(if_else(Group1 == "Control", n >= 2000, n >= 2500)) %>%
  ungroup()

# C10protein_long %>%
#   distinct(Protein, SampleID) %>%
#   mutate(Group1 = case_when(grepl("Trtd", SampleID) ~ "Treated",
#                             grepl("Ctrl", SampleID)~ "Control"
#   )) %>%
#   group_by(SampleID) %>%
#   add_count(name = "n") %>%
#   distinct(SampleID, n, Group1) %>%
#   ungroup()  %>% 
#   ggplot()+
#   aes(y = n, x = Group1, color = Group1)+
#   geom_jitter(position = position_jitter(w = 0.1, h = 0), size = 3)+
#   theme_bw(base_size = 22) +
#   theme(legend.position = "none",
#         panel.background = element_rect(fill= 'white'),
#         axis.text.y=element_text(color = 'black', size = 22),
#         axis.text.x=element_text(angle = 0, vjust = 0.5, hjust= 0.5, color='black',size = 22),
#         panel.grid.minor = element_blank(),
#         panel.grid.major = element_blank(),
#         panel.border = element_blank(),
#         axis.line = element_line(),
#         text=element_text(family="Helvetica")) +
#   ylab("Proteins Identified (n)")+
#   scale_color_manual(values = c("#88419d", "#8c96c6"))+
#   xlab("Cell Type")+
#   geom_segment(aes(x=0.5,xend=1.5),y=2000,yend=2000,
#                linetype='dotted', col = 'red', size = 1)+
#   geom_segment(aes(x=2.5,xend=1.5),y=2500,yend=2500,
#                linetype='dotted', col = 'red', size = 1)+
#   annotate("text", x = 1.5, y = 2000, label = "QC Cutoff", vjust = -1)
#  # ggsave("QC_proteins.png", width = 8, height = 5)

C10protein_wide <- C10protein_long %>%
  filter(SampleID %in% filterS$SampleID) %>%
  spread(SampleID, Intensity) %>% 
  column_to_rownames(var= "Protein") %>%
  mutate(Var1 = rowSums(is.na(.))) %>%
  filter(Var1 <= 34) %>%
  dplyr::select(-Var1)

C10protein_impute <-  DreamAI(C10protein_wide, k = 10, maxiter_MF = 10, ntree = 100, 
                         maxnodes = NULL, maxiter_ADMIN = 30, tol = 10^(-2), 
                         gamma_ADMIN = NA, gamma = 50, CV = FALSE, 
                         fillmethod = "row_mean", maxiter_RegImpute = 10, 
                         conv_nrmse = 1e-06, iter_SpectroFM = 40, method = "KNN", 
                         out = c("KNN")) 

C10protein_impute <- as.data.frame(C10protein_impute$Ensemble)

pca_j<- PCAtools::pca(C10protein_impute, scale = TRUE, center = T)

Group <- as.character(filterS$Group1)


PCAtools::biplot(pca_j, x = "PC1", y =  "PC2", lab = Group,
                 showLoadings = F,
                 ntopLoadings = 2,
                 labSize = 0, drawConnectors = FALSE,
                 legendPosition = "right",
                 encircle = TRUE
                 )
  

  
meta2 <- filterS %>%
  distinct(SampleID, Group1)

C10protein2 <- C10protein_impute %>%
  rownames_to_column(var = "Entry.Name")

C10protein2[,2:69] <- 2^(C10protein2[,2:69])


C10protein2 <- inner_join(C10protein2, meta)


C10protein3 <- C10protein2 %>%
  group_by( Protein) %>%
  add_count(name  = "n") %>%
  filter(n == 1) %>%
  dplyr::select(-n) %>%
  mutate(ID =  Protein) %>%
  mutate(name =  Protein)



LFQ_columns <- grep("[1-9]",colnames(C10protein3)) # get LFQ column numbers
experimental_design <- meta2 %>%
  mutate(label = SampleID) %>%
  mutate(condition = Group1) %>%
  group_by(condition) %>%
  mutate(replicate = dense_rank(dplyr::desc(label))) %>%
  ungroup() %>%
  dplyr::select(label, condition, replicate) %>%
  arrange(desc(label))

data_se <- make_se(C10protein3, LFQ_columns, experimental_design)

data_diff <- DEP::test_diff(data_se, type = "all")

data_diff@elementMetadata@listData <- data_diff@elementMetadata@listData %>%
  as.data.frame() %>%
  mutate(Treated_vs_Control_p.adj = p.adjust(Treated_vs_Control_p.val, method = "BH", n = length(Treated_vs_Control_p.adj))
  ) %>%
  as.list()

dep <- add_rejections(data_diff, alpha = 0.01, lfc = 0.5)



#pdf("cluster_heatmap.pdf", height = 10, width = 5)
p_heat <- DEP::plot_heatmap(dep, type = "centered", kmeans = T,
                  clustering_distance	= "euclidean",
             k = 6, col_limit = 2, show_row_names = F,
             #indicate = "condition",
             column_km = 2,
             show_column_names = FALSE,
             show_column_dend = T,
             show_row_dend = T,
             border = T,
             use_raster = F)
#dev.off()


data_results <- as.data.frame(dep@elementMetadata@listData)

### Clusters from the heatmap are extracted below, along with protein (gene) names 
r.dend <- row_dend(p_heat)


clu_df <- lapply(names(r.dend), function(i){
  out <- data.frame(Protein = names(dendextend:::cutree.dendrogram(r.dend[[i]],k =1)),
                    Cluster = paste0("cluster", i),
                    stringsAsFactors = FALSE)
  return(out)
}) %>%  
  do.call(rbind, .)

clu_df <- inner_join(clu_df, meta)

C10protein_wide2 <-C10protein_long %>%
  filter(SampleID %in% filterS$SampleID) %>%
  dplyr::select(SampleID, Intensity, Protein) %>%
  spread(SampleID, Intensity) %>% 
  column_to_rownames(var= "Protein")

C10protein2 <- C10protein_wide2 %>%
  rownames_to_column(var = "Entry.Name")

C10protein2[,2:69] <- 2^(C10protein2[,2:69])


C10protein2 <- inner_join(C10protein2, meta)

C10protein3 <- C10protein2 %>%
  group_by( Entry.Name) %>%
  add_count(name  = "n") %>%
  filter(n == 1) %>%
  dplyr::select(-n) %>%
  mutate(ID =  Entry.Name) %>%
  mutate(name =  Entry.Name) %>%
  ungroup()


LFQ_columns <- grep("[1-9]",colnames(C10protein3)) # get LFQ column numbers
experimental_design <- meta2 %>%
  mutate(label = SampleID) %>%
  mutate(condition = Group1) %>%
  group_by(condition) %>%
  mutate(replicate = dense_rank(dplyr::desc(label))) %>%
  ungroup() %>%
  dplyr::select(label, condition, replicate) %>%
  arrange(desc(label))

data_se <- make_se(C10protein3, LFQ_columns, experimental_design)

data_diff <- DEP::test_diff(data_se, type = "all")


data_diff@elementMetadata@listData <- data_diff@elementMetadata@listData %>%
  as.data.frame() %>%
  mutate_at(vars(matches("p.val")), ~ p.adjust(.x, method = "BH")) %>%
  dplyr::select(-contains("p.adj"),
                -contains("CI.L"),
                -contains("CI.R")) %>%
  rename_with(.fn = ~ str_replace(., "p.val", "p.adj"),
              .cols = contains("p.val")) %>%
  as.list()

dep <- add_rejections(data_diff, alpha = 0.01, lfc = 0.5)



data_results2 <- as.data.frame(dep@elementMetadata@listData)

save(data_results2, clu_df, file = "Proteomics_results.RData")
save(filterS, file = "SampleID_Filter.RData")


