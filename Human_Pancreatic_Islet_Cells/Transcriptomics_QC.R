library(tidyverse)


###

### Donor 2
secondrun <- read.delim("scRNAseq_Donor2_TPM.tsv", check.names =F) %>%
  mutate(Gene = gsub("gene-","", Gene)) %>%
  pivot_longer(!Gene, names_to = "SampleID", values_to = "TPM") %>%
  mutate(TPM = log2(TPM)) %>%
  filter(TPM != -Inf) %>%
  filter(TPM >= 0) 

### Donor 1
firstrun <- read.delim("scRNAseq_Donor1_TPM.tsv", check.names =F) %>%
  mutate(Gene = gsub("gene-","", Gene)) %>%
  pivot_longer(!Gene, names_to = "SampleID", values_to = "TPM") %>%
  mutate(TPM = log2(TPM)) %>%
  filter(TPM != -Inf) %>%
  filter(TPM >= 0) 


runs_combined <- full_join(firstrun, secondrun)

combined <- runs_combined %>%
  pivot_wider(names_from = "SampleID", values_from = "TPM") %>%
  column_to_rownames(var = "Gene")

check2 <- data.frame(GeneCount = nrow(combined)-colSums(is.na(combined))) %>%
  rownames_to_column(var = "SampleID")


mtgenes <- read.delim("MTgenes.txt",
                      header = F) %>%
  rename(Gene = V1) %>%
  mutate(Gene = gsub("GENE-","", toupper(Gene))) %>%
  mutate(Gene = gsub("MT-","", Gene))


firstrunc <- read.delim("scRNAseq_Donor1_Counts.tsv",
                        check.names =F) %>%
  mutate(Gene = gsub("gene-","", Gene)) %>%
  pivot_longer(!Gene, names_to = "SampleID", values_to = "Counts") %>%
  mutate(Type = case_when(Gene %in% mtgenes$Gene ~ "Mito",
                          TRUE ~ "Not_Mito")) %>%
  filter(Counts >  0)

secondrunc <- read.delim("scRNAseq_Donor2_Counts.tsv",
                         check.names =F)%>%
  mutate(Gene = gsub("gene-","", Gene)) %>%
  pivot_longer(!Gene, names_to = "SampleID", values_to = "Counts") %>%
  mutate(Type = case_when(Gene %in% mtgenes$Gene ~ "Mito",
                          TRUE ~ "Not_Mito")) %>%
  filter(Counts >  0)

runs_combinedc <- full_join(firstrunc, secondrunc)

mitocheck <- runs_combinedc  %>%
  group_by(SampleID) %>%
  mutate(Total = sum(Counts)) %>%
  group_by(SampleID, Type) %>%
  mutate(MitoT = sum(Counts)) %>%
  filter(Type == "Mito") %>%
  ungroup() %>%
  distinct(SampleID, Total, MitoT) %>%
  mutate(percentMito = round(MitoT/Total*100, digits = 1)) %>%
  select(percentMito, SampleID)

check3 <- inner_join(check2, mitocheck)

filterS <- check3 %>%
  mutate(Type = case_when(percentMito <= 0.99 ~"PoorMito%",
                          percentMito >= 30 ~ "PoorMito%",
                          GeneCount < 400 ~ "LowDepth",
                          TRUE ~ "QCPassed")) %>%
  filter(Type == "QCPassed")


subset2 <- runs_combined %>%
  filter(SampleID %in% filterS$SampleID) %>%
  pivot_wider(names_from = "SampleID", values_from = "TPM") %>%
  column_to_rownames(var = "Gene")


check3 %>%
  mutate(`QC Metric` = case_when(percentMito <= 0.99 ~"Low Mitochondrial Reads",
                          percentMito >= 30 ~ "High Mitochondrial Reads",
                          GeneCount < 400 ~ "Low Identification Depth",
                          TRUE ~ "QC Passed")) %>%
  mutate(Donor = case_when(grepl("227", SampleID)~ "Donor 1",
                           TRUE ~ "Donor 2")) %>%
  ggplot()+
  aes(x = percentMito, y = GeneCount, color = `QC Metric`,
      shape = Donor)+
  geom_point(size = 3)+
  theme_classic(base_size = 18)+
  ylab("Gene Identifications (n)")+
  xlab("Mitochondrial Reads (%)")

geneF <- runs_combined %>%
  filter(SampleID %in% filterS$SampleID) %>%
  group_by(Gene) %>%
  tally() %>%
  ungroup() %>%
  distinct(Gene) %>%
  filter(Gene != "GTF2H2C-2")

 seuratislet <- runs_combinedc %>%
  filter(SampleID %in% filterS$SampleID) %>%
  filter(Gene %in% geneF$Gene) %>%
  select(SampleID, Counts, Gene) %>%
  pivot_wider(names_from = "SampleID", values_from = "Counts") %>%
  column_to_rownames(var = "Gene")
 
# save(seuratislet, file = "SeuratIslet.RData")
 
