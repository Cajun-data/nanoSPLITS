library(tidyverse)
library(SingleCellExperiment)
library(scp)
library(ggplot2)


mod_pep <- read_tsv("C10SVEC_singlecells_Peptide_intensities.tsv") 


peptide1 <- mod_pep %>%
  filter(grepl("MOUSE", Protein)) %>%
  dplyr::select(`Modified.Sequence` | contains("_")) %>%
  column_to_rownames(var="Modified.Sequence")
###################
mod_pep <- read_tsv("C10Treated_Peptide_intensities.tsv") 
meta1 <- read.csv("C10Treated_NormCounts.csv", check.names = T) %>%
  select(contains("Ctrl")| contains("Trtd")) %>%
  colnames(.)
peptide2 <- mod_pep %>%
  filter(grepl("MOUSE", Protein)) %>%
  dplyr::select(`Modified.Sequence` | contains("Ctrl")| contains("Trtd")) %>%
  column_to_rownames(var="Modified.Sequence") %>%
  dplyr::select((meta1))
######################
mod_pep <- read_tsv("combined_modified_peptide_MaxLFQ.tsv") 
meta2 <- read.csv("annotations.csv")
peptide3 <- mod_pep %>%
  filter(grepl("HUMAN", Protein)) %>%
  dplyr::select(`Modified.Sequence` | contains("_")) %>%
  column_to_rownames(var="Modified.Sequence") %>%
  select(contains(meta2$SampleID))

peptide1 <- peptide1 != 0
peptide2 <- peptide2 != 0
peptide3 <- peptide3 != 0

sampleGroup1 <- data.frame(SampleID = colnames(peptide1)) %>%
  mutate(Group = case_when(grepl("_C10_", SampleID)~ "C10",
                           grepl("_SVEC_", SampleID)~ "SVEC"))
sampleGroup2 <- data.frame(SampleID = colnames(peptide2)) %>%
  mutate(Group = "CDK1_Experiment")
sampleGroup3 <- data.frame(SampleID = colnames(peptide3)) %>%
  mutate(Group = "Pancreatic_Cells")

sce <- SingleCellExperiment(assays = peptide1)
sce2 <- SingleCellExperiment(assays = peptide2)
sce3 <- SingleCellExperiment(assays = peptide3)
sim <- QFeatures(experiments = List(Assay1 = sce,
                                    Assay2 = sce2,
                                    Assay3 = sce3))


pep_Missing1 <- reportMissingValues(sim, "Assay1", by = sampleGroup1$Group) 
pep_Missing2 <- reportMissingValues(sim, "Assay2", by = sampleGroup2$Group)
pep_Missing3 <- reportMissingValues(sim, "Assay3", by = sampleGroup3$Group)

csc1 <- cumulativeSensitivityCurve(
  sim, "Assay1", by = sampleGroup1$Group,
) 
csc2 <- cumulativeSensitivityCurve(
  sim, "Assay2", by = sampleGroup2$Group,
) 
csc3 <- cumulativeSensitivityCurve(
  sim, "Assay3", by = sampleGroup3$Group,
) 

ji1 <- jaccardIndex(sim, "Assay1", by = sampleGroup1$Group)
ji2 <- jaccardIndex(sim, "Assay2", by = sampleGroup2$Group) 
ji3 <- jaccardIndex(sim, "Assay3", by = sampleGroup3$Group)

csc <- do.call("rbind", list(csc1,csc2,csc3))


csc %>%
  ggplot()+
  aes(x = SampleSize, y = Sensitivity, colour = by) +
  geom_point()+
  geom_hline(yintercept = pep_Missing1$LocalSensitivityMean,
             linetype = "dashed",
             color = c("#F8766D","#C77CFF"))+
  #geom_hline(yintercept = pep_Missing1$TotalSensitivity,
   #          color = c("#F8766D","#C77CFF"))+
  geom_hline(yintercept = pep_Missing3$LocalSensitivityMean,
             linetype = "dashed",
             color = "#00BFC4")+
 # geom_hline(yintercept = pep_Missing3$TotalSensitivity,
    #         color = "#00BFC4")+
  geom_hline(yintercept = pep_Missing2$LocalSensitivityMean,
             linetype = "dashed",
             color = "#7CAE00")+
 # geom_hline(yintercept = pep_Missing2$TotalSensitivity,
    #         color = "#7CAE00")+
  scale_y_continuous(limits = c(0,26000))+
  theme_bw(base_size = 20) +
  theme(#legend.position = "none",
        panel.background = element_rect(fill= 'white'),
        axis.text.y=element_text(color = 'black', size = 20),
        axis.text.x=element_text(angle = 0, vjust = 0.5, hjust= 0.5, color='black',size = 20),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line())+
  ylab("Sensitivity (Peptides ,n)")+
  xlab("Number of Cells")+
  ggsave("sensitivity_index.png", width = 8, height = 4.7)

  
  ji <- do.call("rbind", list(ji1,ji2,ji3))
  

ji %>%
  ggplot()+
  aes(y = jaccard, x = by, fill = by) +
  geom_violin()+
  scale_y_continuous(limits = c(0,1))+
    theme_bw(base_size = 15) +
    theme(legend.position = "none",
          panel.background = element_rect(fill= 'white'),
          axis.text.y=element_text(color = 'black', size = 12),
          axis.text.x=element_text(angle = 0, vjust = 0.5, hjust= 0.5, color='black',size = 15),
          panel.grid.minor = element_blank(),
          panel.grid.major = element_blank(),
          panel.border = element_blank(),
          axis.line = element_line())+
  ylab("Jaccard Index")+
  ggsave("jaccard_index.png", width = 8, height = 4.7)



data_completedness <- do.call("rbind", list(pep_Missing1,pep_Missing2,pep_Missing3)) %>%
  mutate(Completeness = Completeness * 100) %>%
  mutate(Type = rownames(.))


 data_completedness %>%
  ggplot()+
  aes(x = Completeness, y = LocalSensitivityMean, color = Type,
      group = Type,
      size = NumberCells)+
  geom_point(size =10,
             shape = 5)+
  geom_point()+
  theme_bw()+
  scale_x_continuous(limits = c(0,90))+
  scale_y_continuous(limits = c(0,20000))+
  geom_pointrange(aes(ymin=LocalSensitivityMean-LocalSensitivitySd,
                    ymax=LocalSensitivityMean+LocalSensitivitySd),
                position=position_dodge(0.05))+
   theme_bw(base_size = 20) +
   theme(panel.background = element_rect(fill= 'white'),
         axis.text.y=element_text(color = 'black', size = 20),
         axis.text.x=element_text(angle = 0, vjust = 0.5, hjust= 0.5, color='black',size = 20),
         panel.grid.minor = element_blank(),
         panel.grid.major = element_blank(),
         panel.border = element_blank(),
         axis.line = element_line(),
         text=element_text(family="Helvetica"))+
   ylab("Local Sensitivty")+
   xlab("Data completeness (%)")+
   scale_size(range = c(0.27,1.2),
              breaks = c(31,39,62,128))+
   ggsave("data_completeness.png", width = 8, height = 4.7)

