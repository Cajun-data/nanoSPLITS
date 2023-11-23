library(tidyverse)
library(proDA)
library(DreamAI)
library(PCAtools)
library(DEP)
library(ComplexHeatmap)

load("SampleID_Filter.RData")

phospho <- read.delim("C10Treated_Peptide_intensities.tsv") 


meta <- phospho[,1:8]

meta2 <- meta %>%
  dplyr::select(Modified.Sequence, Entry.Name)

phospho <- phospho %>%
  dplyr::select(-Peptide.Sequence, -Assigned.Modifications, -Protein, -ID, -Entry.Name,
         -Gene, -Protein.Description) %>%
  column_to_rownames(var = "Modified.Sequence")

phospho[phospho == 0] <- NA
phospho[,1:96] <- log2(phospho[,1:96]) 
phospho <- as.data.frame(median_normalization(as.matrix(phospho))) %>%
  rownames_to_column(var = "Modified.Sequence")


phospho_long <- pivot_longer(phospho, -Modified.Sequence, names_to = "SampleID",
                             values_to = "Intensity",
                             names_repair = "check_unique") %>%
  filter(!is.na(Intensity)) %>%
  inner_join(., meta2) %>%
  filter(grepl("_MOUSE", Entry.Name)) %>%
  dplyr::select(-Entry.Name)


phospho_wide <- phospho_long %>%
  filter(SampleID %in% filterS$SampleID) %>%
  spread(SampleID, Intensity) %>% 
  column_to_rownames(var= "Modified.Sequence") %>%
  mutate(Var1 = rowSums(is.na(.))) %>%
  filter(Var1 <= 34) %>%
  dplyr::select(-Var1)

phospho_impute <-  DreamAI(phospho_wide, k = 10, maxiter_MF = 10, ntree = 100, 
                           maxnodes = NULL, maxiter_ADMIN = 30, tol = 10^(-2), 
                           gamma_ADMIN = NA, gamma = 50, CV = FALSE, 
                           fillmethod = "row_mean", maxiter_RegImpute = 10, 
                           conv_nrmse = 1e-06, iter_SpectroFM = 40, method = "KNN", 
                           out = c("KNN")) 

phospho_impute <- as.data.frame(phospho_impute$Ensemble)


meta2 <- filterS %>%
  distinct(SampleID, Group1)

phospho2 <- phospho_impute %>%
  rownames_to_column(var = "Modified.Sequence")

phospho2[,2:69] <- 2^(phospho2[,2:69])


phospho2 <- inner_join(phospho2, meta)

phospho3 <- phospho2 %>%
  group_by(Modified.Sequence) %>%
  add_count(name  = "n") %>%
  filter(n == 1) %>%
  dplyr::select(-n) %>%
  mutate(ID =  Modified.Sequence) %>%
  mutate(name =  Modified.Sequence)



LFQ_columns <- grep("[1-9]",colnames(phospho3)) # get LFQ column numbers
experimental_design <- meta2 %>%
  mutate(label = SampleID) %>%
  mutate(condition = Group1) %>%
  group_by(condition) %>%
  mutate(replicate = dense_rank(dplyr::desc(label))) %>%
  ungroup() %>%
  dplyr::select(label, condition, replicate) %>%
  arrange(desc(label))

data_se <- make_se(phospho3, LFQ_columns, experimental_design)

data_diff <- DEP::test_diff(data_se, type = "all")

dep <- add_rejections(data_diff, alpha = 0.01, lfc = log2(1.0))

dep@elementMetadata@listData <- dep@elementMetadata@listData %>%
  as.data.frame() %>%
  mutate(Treated_vs_Control_p.adj = p.adjust(Treated_vs_Control_p.val, method = "BH", n = length(Treated_vs_Control_p.adj))
  ) %>%
  mutate(Treated_vs_Control_significant = case_when(Treated_vs_Control_p.adj <=0.01 & !between(Treated_vs_Control_diff, -0.5, 0.5) ~ TRUE,
                                                    TRUE ~ FALSE)) %>%
  mutate(significant = case_when(Treated_vs_Control_significant == TRUE ~ TRUE,
                                 TRUE ~ FALSE)) %>%
  as.list()


peptide_results <- as.data.frame(dep@elementMetadata@listData)

save(peptide_results, file = "Phosphopeptide_DE.RData")
