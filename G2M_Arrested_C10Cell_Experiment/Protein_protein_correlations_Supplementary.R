library(tidyverse)
library(proDA)
library(corrplot)
library(ggpubr)

##############################
load("SampleID_Filter.RData")
load("Proteomics_results.RData")

C10treated <- read.delim("C10Treated_Protein_intensities.tsv")

meta <- C10treated[,1:4] %>%
  filter(grepl("_MOUSE", Protein)) 


C10treated <- C10treated %>%
  dplyr::select(-Entry.Name, -Protein, -Protein.ID ) %>%
  distinct(Gene, .keep_all = TRUE) %>%
  column_to_rownames(var = "Gene") %>%
  dplyr::select(one_of(filterS$SampleID))

  

C10treated[C10treated == 0] <- NA
C10treated[,1:68] <- log2(C10treated[,1:68]) 
C10treated <- as.data.frame(median_normalization(as.matrix(C10treated))) %>%
  rownames_to_column(var = "Gene")


C10treated_long <- pivot_longer(C10treated, -Gene, names_to = "SampleID",
                                values_to = "Intensity",
                                names_repair = "check_unique") %>%
  filter(!is.na(Intensity))  


#pdf(file = "Sparc_vs_Col1a1_protein_correlations.pdf", height = 4, width = 4)
C10treated_long %>%
  mutate(Type = case_when(grepl("Ctrl", SampleID) ~ "Untreated",
                                grepl("Trtd", SampleID) ~ "G2/M Arrested")) %>%
  filter(Gene ==  "Psap" | Gene == "Grn") %>%
  spread(Gene, Intensity) %>%
  ggscatter( x = "Psap", y = "Grn",
             add = "reg.line",                         # Add regression line
             conf.int = TRUE,                          # Add confidence interval
             color = "Type", 
             palette = "jco",           # Color by groups "cyl"
              shape = "Type"                             # Change point shape by groups "cyl"
  )+
  stat_cor(aes(color = Type),
           label.x.npc = "middle",
           label.y.npc = "bottom")+
  stat_cor(aes())+
  xlab("log2(Intensity, Psap)")+
  ylab("log2(Intensity, Grn)")
#dev.off()




C10treated_long %>%
  mutate(Type = case_when(grepl("Ctrl", SampleID) ~ "Untreated",
                          grepl("Trtd", SampleID) ~ "G2/M Arrested")) %>%
  filter(Gene ==  "Ctsa" | Gene == "Ctsb") %>%
  spread(Gene, Intensity) %>%
  ggscatter( x = "Ctsa", y = "Ctsb",
             add = "reg.line",                         # Add regression line
             conf.int = TRUE,                          # Add confidence interval
             color = "Type", 
             palette = "jco",           # Color by groups "cyl"
             shape = "Type"                             # Change point shape by groups "cyl"
  )+
  stat_cor(aes(color = Type),
           label.x.npc = "middle",
           label.y.npc = "bottom")+
  stat_cor(aes())+
  xlab("log2(Intensity, Ctsa)")+
  ylab("log2(Intensity, Ctsb)")





C10treated_long %>%
  mutate(Type = case_when(grepl("Ctrl", SampleID) ~ "Untreated",
                          grepl("Trtd", SampleID) ~ "G2/M Arrested")) %>%
  filter(Gene ==  "Kpna2" | Gene == "Top2a") %>%
  spread(Gene, Intensity) %>%
  ggscatter( x = "Kpna2", y = "Top2a",
             add = "reg.line",                         # Add regression line
             conf.int = TRUE,                          # Add confidence interval
             color = "Type", 
             palette = "jco",           # Color by groups "cyl"
             shape = "Type"                             # Change point shape by groups "cyl"
  )+
  stat_cor(aes(color = Type),
           label.x.npc = "middle",
           label.y.npc = "bottom")+
  stat_cor(aes())+
  xlab("log2(Intensity, Kpna2)")+
  ylab("log2(Intensity, Top2a)")







C10treated_long %>%
  mutate(Type = case_when(grepl("Ctrl", SampleID) ~ "Untreated",
                          grepl("Trtd", SampleID) ~ "G2/M Arrested")) %>%
  filter(Gene ==  "Tpx2" | Gene == "Top2a") %>%
  spread(Gene, Intensity) %>%
  ggscatter( x = "Tpx2", y = "Top2a",
             add = "reg.line",                         # Add regression line
             conf.int = TRUE,                          # Add confidence interval
             color = "Type", 
             palette = "jco",           # Color by groups "cyl"
             shape = "Type"                             # Change point shape by groups "cyl"
  )+
  stat_cor(aes(color = Type),
           label.x.npc = "middle",
           label.y.npc = "bottom")+
  stat_cor(aes())+
  xlab("log2(Intensity, Tpx2)")+
  ylab("log2(Intensity, Top2a)")
