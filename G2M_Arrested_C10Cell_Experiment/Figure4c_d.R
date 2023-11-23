library(tidyverse)
library(edgeR)
library(proDA)
library(ggpubr)
library(rstatix)



#############################
load("Proteomics_results.RData")
C10protein <- read.delim("C10Treated_Protein_intensities.tsv")
load("SampleID_Filter.RData")
load("Phosphopeptide_DE.RData")
phospho <- read.delim("C10Treated_Peptide_intensities.tsv") 
RNA_norm3 <- read.csv("C10Treated_NormCounts.csv",
                      check.names = TRUE) %>%
  column_to_rownames(var = "X")
load("scRNAseq_DE.RData")


#################################
meta <- C10protein[,1:4] %>%
  filter(grepl("_MOUSE", Protein)) 


C10protein <- C10protein %>%
  distinct(Gene, .keep_all = TRUE) %>%
  column_to_rownames(var = "Gene") %>%
  dplyr::select(-Protein, -Protein.ID, -Entry.Name) %>%
  dplyr::select(one_of(filterS$SampleID))

C10protein[C10protein == 0] <- NA
C10protein[,1:68] <- log2(C10protein[,1:68]) 
C10protein <- median_normalization(as.matrix(C10protein)) 

C10protein <- as.data.frame(t(scale(t(C10protein)))) %>%
  rownames_to_column(var = "Gene")




C10protein_long <- pivot_longer(C10protein, -Gene, names_to = "SampleID",
                                values_to = "Protein",
                                names_repair = "check_unique") %>%
  filter(Gene %in% data_results2$Gene) %>%
  dplyr::select(Gene, SampleID, Protein)


#############################

meta2 <- phospho[,1:8]

phospho <- phospho %>%
  dplyr::select(-Peptide.Sequence, -Assigned.Modifications, -Protein, -ID, -Entry.Name,
                -Gene, -Protein.Description) %>%
  column_to_rownames(var = "Modified.Sequence") %>%
  dplyr::select(one_of(filterS$SampleID))

phospho[phospho == 0] <- NA
phospho[,1:68] <- log2(phospho[,1:68]) 
phospho <- median_normalization(as.matrix(phospho)) 
phospho <- as.data.frame(t(scale(t(phospho)))) %>%
  rownames_to_column(var = "Modified.Sequence") %>%
  filter(grepl("79.9", Modified.Sequence)) 

phospho_long <- pivot_longer(phospho, -Modified.Sequence, names_to = "SampleID",
                                values_to = "Phosphopeptide",
                                names_repair = "check_unique") %>%
  inner_join(.,meta2) %>%
  filter(Gene %in% peptide_results$Gene) %>%
  dplyr::select(Modified.Sequence, Gene, SampleID,Phosphopeptide)


check <- full_join(C10protein_long, phospho_long)


############################
RNA_mat <- RNA_norm3 

RNA_mat<- cpm(RNA_mat, log = T)
RNA_mat[RNA_mat <= 0.004] <- NA
RNA_mat <- as.data.frame(t(scale(t(RNA_mat)))) %>%
  rownames_to_column(var = "Gene")


mRNA_long <- pivot_longer(RNA_mat, -Gene, names_to = "SampleID",
                             values_to = "mRNA",
                             names_repair = "check_unique") %>%
  filter(!is.na(mRNA)) %>%
  filter(Gene %in% results2$Gene) %>%
  dplyr::select(Gene, SampleID,mRNA)




############################ Combine
gene1 <- data_results2$Gene
gene2 <- peptide_results$Gene
gene3 <- results2$Gene
gene4 <- Reduce(intersect, list(gene1,gene2,gene3))

check2 <- full_join(check, mRNA_long) %>%
  filter(Gene %in% gene4) %>%
  pivot_longer(cols = c("Protein", "mRNA", "Phosphopeptide"),
               names_to = "Type",
               values_to = "Zscore",
               values_drop_na = T) %>%
  group_by(Gene, SampleID, Type) %>%
  mutate(Rcheck = order(Modified.Sequence)) %>%
  mutate(Gene = paste(Gene, Rcheck, sep ="_")) %>%
  mutate(Group = case_when(grepl("Ctrl", SampleID) ~ "Untreated",
                           TRUE ~ "G2/M Arrested")) %>%
  mutate(Group2 = paste(Group, Type, sep = "_"))




remove1 <- check2 %>%
  mutate(Group = case_when(grepl("Ctrl", SampleID) ~ "Untreated",
                           TRUE ~ "G2/M Arrested")) %>%
  group_by(Gene,Type,Group) %>%
  add_count(name = "n") %>%
  filter(n >= 5) %>%
  dplyr::select(-n) %>%
  distinct(Gene,Type, Group) %>%
  group_by(Gene, Type) %>%
  tally() %>%
  filter(n == 1)


stat.test <- check2 %>%
  mutate(Group = case_when(grepl("Ctrl", SampleID) ~ "Untreated",
                           TRUE ~ "G2/M Arrested")) %>%
  group_by(Gene,Type,Group) %>%
  add_count(name = "n") %>%
  filter(n >= 5) %>%
  dplyr::select(-n) %>%
  filter(!Gene %in% remove1$Gene) %>%
  group_by(Gene, Type) %>%
  wilcox_test(Zscore ~ Group) %>%
  adjust_pvalue(method = "BH") %>%
  add_xy_position(x = "Group", fun ="median_iqr") %>%
  add_significance()

stat.test <- stat.test %>%
  mutate(y.position = y.position + 0.1)

stat.test$type_f <- factor(stat.test$Type,
                        levels = c("Protein", "mRNA",
                                   "Phosphopeptide"))


check2$type_f <- factor(check2$Type,
                        levels = c("Protein", "mRNA",
                                   "Phosphopeptide"))

#pdf("Vimentin_compareplot.pdf", width = 7, height = 3.7)
check2 %>%
  filter(Gene == "Vim_1") %>%
  mutate(Group = case_when(grepl("Ctrl", SampleID) ~ "Untreated",
                           TRUE ~ "G2/MArrested")) %>%
  mutate(Group2 = paste(Group, Type, sep = "_")) %>%
  group_by(Group) %>%
  ggplot()+
  aes(x = Group, y = Zscore)+
  geom_boxplot(alpha = 0.5, aes(fill = Group),
               outlier.shape = NA)+
  theme_bw(base_size = 20)+
  xlab(NULL)+
  ylab("z-score")+
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        axis.line = element_line(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        text=element_text(family="Helvetica"))+
  scale_y_continuous(limits = c(-3, 3))+
  ggplot2::facet_wrap(.~type_f)+
  stat_pvalue_manual(size = 6, subset(stat.test, Gene == "Vim_1"))
#dev.off()


#pdf("Hnrpu1_compareplot.pdf", width = 7, height = 3.7)
check2 %>%
  filter(Gene == "Hnrnpu_1") %>%
  mutate(Group = case_when(grepl("Ctrl", SampleID) ~ "Untreated",
                           TRUE ~ "G2/M Arrested")) %>%
  mutate(Group2 = paste(Group, Type, sep = "_")) %>%
  group_by(Group) %>%
  ggplot()+
  aes(x = Group, y = Zscore)+
  geom_boxplot(alpha = 0.5, aes(fill = Group),
               outlier.shape = NA)+
  theme_bw(base_size = 20)+
  xlab(NULL)+
  ylab("z-score")+
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        axis.line = element_line(),
        legend.position = "n",
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        text=element_text(family="Helvetica"))+
  scale_y_continuous(limits = c(-3, 3))+
  ggplot2::facet_wrap(.~type_f)+
  stat_pvalue_manual(size = 6, subset(stat.test, Gene == "Hnrnpu_1"))
#dev.off()

