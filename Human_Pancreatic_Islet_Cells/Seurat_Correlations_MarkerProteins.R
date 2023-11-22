library(tidyverse)
library(Seurat)
library(Azimuth)
library(SeuratData)



impute_normal <- function(object, width=0.3, downshift=1.8, seed=100) {
  
  if (!is.matrix(object)) object <- as.matrix(object)
  mx <- max(object, na.rm=TRUE)
  mn <- min(object, na.rm=TRUE)
  if (mx - mn > 20) warning("Please make sure the values are log-transformed")
  
  set.seed(seed)
  object <- apply(object, 2, function(temp) {
    temp[!is.finite(temp)] <- NA
    temp_sd <- stats::sd(temp, na.rm=TRUE)
    temp_mean <- mean(temp, na.rm=TRUE)
    shrinked_sd <- width * temp_sd   # shrink sd width
    downshifted_mean <- temp_mean - downshift * temp_sd   # shift mean of imputed values
    n_missing <- sum(is.na(temp))
    temp[is.na(temp)] <- stats::rnorm(n_missing, mean=downshifted_mean, sd=shrinked_sd)
    temp
  })
  return(object)
}



#####################################
pancREF <- readRDS(file = "ref.Rds")

load(file  = "SeuratIslet.RData")
seuratislet  <- seuratislet %>%
  as.sparse()
seuratislet[is.na(seuratislet)] <- 0

# Create Seurat object
isletData <- CreateSeuratObject(seuratislet)
isletData <- NormalizeData(isletData)
isletData <- FindVariableFeatures(isletData)
isletData <- ScaleData(isletData)
# The RunAzimuth function can take a Seurat object as input
isletData <- RunAzimuth(isletData, reference = "pancreasref",
                         k.weight = 5,
                        n.trees = 20,
                         mapping.score.k = 100)

isletData2 <-  subset(isletData, subset = predicted.annotation.l1.score >= 0.8)


p1 <- DimPlot(isletData2, group.by = "predicted.annotation.l1", label = TRUE,
              pt.size = 3,
              label.size = 4,
              repel = T) + NoLegend() + ggtitle("nanoSPLITS Cells")

p3 <- DimPlot(pancREF, reduction = "refUMAP",
              #group.by = 'predicted.annotation.l1',
              cols = c('alpha' = '#93AA00',
                       'acinar' = "#F8766D",
                       'beta' = '#00BA38',
                       "delta" = "#00C19F",
                       "ductal" = "#00B9E3" ,
                       "endothelial" = "#619CFF",
                       "immune" = "#FF61C3",
                       "gamma" = "#DB72FB",
                       "activated_stellate" = "#D39200"),
              label = TRUE, label.size = 4,
              repel = T) + NoLegend()+ ggtitle("Azimuth Pancreas Reference")



annotations_csv <- as.data.frame(isletData@assays[["prediction.score.annotation.l1"]]@data)

annotations_csv <-  annotations_csv %>%
  rownames_to_column(var = "Celltype") %>%
  pivot_longer(-Celltype, names_to = "SampleID",
               values_to = "AnnotationScore") %>%
  filter(AnnotationScore != 0) %>%
  group_by(SampleID) %>%
  slice_max(order_by = AnnotationScore, n = 1) %>%
  ungroup()

#write.csv(annotations_csv, file = "annotations.csv")





########## Combine with proteomics
### start with panc_pept_batch object from "Panc_Islet_Cells_Workflow.R"
panc_pept <- read.delim("combined_modified_peptide_MaxLFQ.tsv") %>%
  filter(!grepl("contam", Protein)) %>%
  filter(grepl("HUMAN", Protein)) %>%
  filter(!grepl("Keratin", Protein.Description)) 

annot <- read.csv("LDA_Celltypes.csv",
                  check.names = T)

meta <- panc_pept[c(2,9,10,11,13)]
###139 samples
panc_pept <- panc_pept %>%
  column_to_rownames(var = "Modified.Sequence") %>%
  dplyr::select(contains("2")) %>%
  log2() %>%
  rownames_to_column(var  = "Modified.Sequence") %>%
  pivot_longer(-`Modified.Sequence`,
               names_to = "SampleID",
               values_to = "Intensity") %>%
  filter(Intensity != -Inf) %>%
  mutate(SampleID = gsub("X","", SampleID)) %>%
  group_by(SampleID) %>%
  add_count(name = "TotalPeptideIDs") %>%
  filter(TotalPeptideIDs > 3000) %>% ### remove cells with fewer than 3000 peptides
  ungroup() %>%
  dplyr::select(SampleID,Modified.Sequence,Intensity) %>%
  pivot_wider(names_from = "SampleID",
              values_from = "Intensity") %>%
  column_to_rownames(var = "Modified.Sequence") 

###126 samples
batch_meta <- data.frame(SampleID = colnames(panc_pept),
                         Batch = NA)
batch_meta <- batch_meta %>%
  mutate(Batch = case_when( grepl("22", SampleID) ~ "Donor 1",
                            grepl("23", SampleID)~ "Donor 2"
  )) 


panc_pept_batch <- panc_pept %>%
  MSnSet.utils::ComBat.NA(., batch_meta$Batch) 

combined <- panc_pept_batch$`corrected data` %>%
  as.data.frame()  %>%
  rownames_to_column(var  = "Modified.Sequence") %>%
  pivot_longer(-`Modified.Sequence`,
               names_to = "SampleID",
               values_to = "Intensity") %>%
  full_join(.,meta) %>%
  dplyr::select(SampleID,Assigned.Modifications, Modified.Sequence, Protein.ID, Gene, Intensity) %>%
  mutate(Protein = case_when(Gene == "" ~ Protein.ID,
                             TRUE ~ Gene)) %>%
  distinct(Protein, Modified.Sequence, SampleID, Intensity) %>%
  dplyr::rename(protein_list = Protein,
                sample_list = SampleID,
                id = Modified.Sequence,
                quant = Intensity) %>%
  filter(!is.na(quant)) %>%
  as.list() %>%
  iq::fast_MaxLFQ(.) %>%
  .[[1]] 

isletData3 <- subset(isletData2, cells =  colnames(combined))

combined2 <- combined %>%
  as.data.frame() %>%
  dplyr::select(ends_with(colnames(isletData3))) %>%
  as.sparse()

prot_assay <- CreateAssayObject(counts = combined2)
isletData3[["Prot"]] <- prot_assay

DefaultAssay(isletData3) <- 'Prot'
isletData3 <-  ScaleData(isletData3)

FeaturePlot(isletData3, label = F,
            features = c("prot_SST",
                         "prot_INS",
                         "prot_GCG",
                         "prot_PPY",
                         "prot_ACP5",
                         "prot_COL1A1"),
            keep.scale = NULL,
            slot = "data",
           cols = c("gray100","blue"),
            order = T,
            pt.size = 3)


FeaturePlot(isletData3, label = F,
            features = c("prot_TIMP1",
                         "prot_TIMP2",
                         "prot_PTPRN",
                         "prot_PTPRN2"),
            keep.scale = NULL,
            slot = "data",
            cols = c("gray100","blue"),
            order = T,
            pt.size = 3)


#########################
#########################

TPM1 <- read.delim("scRNAseq_Donor1_TPM.tsv", row.names = 1,
                   check.names = F) %>%
  dplyr::select(ends_with(colnames(isletData3))) %>%
  log2() %>%
  rownames_to_column(var = "Gene") %>%
  pivot_longer(-Gene,
               names_to = "SampleID",
               values_to = "TPM") %>%
  filter(TPM != -Inf) %>%
  inner_join(., annot) %>%
  dplyr::select(Gene, SampleID, Celltype, TPM) 
TPM2 <- read.delim("scRNAseq_Donor2_TPM.tsv", row.names = 1,
                   check.names = F) %>%
  dplyr::select(ends_with(colnames(isletData3))) %>%
  log2() %>%
  rownames_to_column(var = "Gene") %>%
pivot_longer(-Gene,
             names_to = "SampleID",
             values_to = "TPM") %>%
  filter(TPM != -Inf) %>%
  inner_join(., annot) %>%
  dplyr::select(Gene, SampleID, Celltype, TPM)


combined2 <- combined %>%
  as.data.frame() %>%
  dplyr::select(ends_with(colnames(isletData3))) %>%
  rownames_to_column(var = "Gene") %>%
  pivot_longer(-Gene,
               names_to = "SampleID",
               values_to = "Intensity") %>%
  inner_join(., annot) %>%
  dplyr::select(Gene, SampleID, Celltype, Intensity) %>%
  filter(!is.na(Intensity))


TPM <- rbind(TPM1, TPM2) %>%
  mutate(Gene = gsub("gene-", "", Gene))

datatest <- inner_join(combined2, TPM) %>%
  group_by(Gene) %>%
  add_count(name = "Obs") %>%
  filter(Obs >= 4) %>%
  ungroup()


datatest <- as_tibble(datatest)


data_nest <- group_by(datatest, Gene) %>% nest()

cor_fun <- function(df) cor.test(df$Intensity, df$TPM, method = "pearson") %>%
  broom::tidy()

data_nest <- mutate(data_nest, model = map(data, cor_fun))
corr_pr <- dplyr::select(data_nest, -data) %>% unnest(cols = c(model))
corr_pr <- corr_pr %>%
  ungroup() %>%
  mutate(FDR = p.adjust(`p.value`, method = "BH")) %>%
  mutate(Type2 = "mRNA-Protein")

set.seed(420)
datatest <- inner_join(combined2, TPM) %>%
  group_by(Gene) %>%
  add_count(name = "Obs") %>%
  filter(Obs >= 4) %>%
  ungroup() %>%
  transform(TPM = sample(TPM)) %>%
  transform(Intensity = sample(Intensity)) %>%
  ungroup()

datatest <- as_tibble(datatest)

data_nest <- group_by(datatest, Gene) %>% nest()

cor_fun <- function(df) cor.test(df$Intensity, df$TPM, method = "pearson") %>%
  broom::tidy()

data_nest <- mutate(data_nest, model = map(data, cor_fun))
corr_pr2 <- dplyr::select(data_nest, -data) %>% unnest(cols = c(model))
corr_pr2 <- corr_pr2 %>%
  ungroup() %>%
  mutate(FDR = p.adjust(`p.value`, method = "BH")) %>%
  mutate(Type2 = "Randomized")

combined_pearson <- full_join(corr_pr, corr_pr2)


combined_pearson %>%
  ggplot()+
  aes(x = estimate, fill = Type2)+
  geom_histogram(alpha = 0.6,
                 position = "identity",
                 bins = 100)+
  xlab("Pearson Correlation (r)")+
  ylab("Gene-Protein Pair (n)")+
  theme_minimal(base_size = 16)+
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line())+
  scale_fill_brewer(name = "", palette = 2,type = "qual",
                    direction = -1)+
  theme(legend.position = "bottom",
        legend.direction = "vertical",
        legend.box = "horizontal")+
  geom_vline(xintercept = 0.0552) ### median of experimental mRNA-Protein correlations


wilcox.test(corr_pr$estimate, corr_pr2$estimate, alternative = "two.sided")


###########################
############################



########## Combine with proteomics
### start with panc_pept_batch object from "Panc_Islet_Cells_Workflow.R"


annot <- read.csv("LDA_Celltypes.csv",
                  check.names = T) %>%
  arrange(SampleID)



combined <- panc_pept_batch$`corrected data` %>%
  as.data.frame()  %>%
  rownames_to_column(var  = "Modified.Sequence") %>%
  pivot_longer(-`Modified.Sequence`,
               names_to = "SampleID",
               values_to = "Intensity") %>%
  inner_join(.,meta) %>%
  mutate(Protein = case_when(Gene == "" ~ Protein.ID,
                             TRUE ~ Gene)) %>%
  distinct(Protein, Modified.Sequence, SampleID, Intensity) %>%
  dplyr::rename(protein_list = Protein,
                sample_list = SampleID,
                id = Modified.Sequence,
                quant = Intensity) %>%
  filter(!is.na(quant)) %>%
  as.list() %>%
  iq::fast_MaxLFQ(.) %>%
  .[[1]]

combined <- combined[,order(colnames(x = combined))]

combined3 <- combined %>%
  impute_normal(width=0.3, downshift=3) %>% ###Need imputation...
  as.sparse()


isletProt <- CreateSeuratObject(combined3)
isletProt <- RenameAssays(
  isletProt,
  assay.name = "RNA",
  new.assay.name = "Prot",
  verbose = TRUE
)

isletProt <-  ScaleData(isletProt)


isletProt <- AddMetaData(
  object = isletProt,
  metadata = annot$Celltype,
  col.name = 'CellType'
)

Idents(object = isletProt, cells = 1:126) <- annot$Celltype

Prot.markers <- FindAllMarkers(isletProt, assay = "Prot", only.pos = TRUE,
                                  test.use = "wilcox",
                               slot = "scale.data")

top7_immune <- Prot.markers %>%  filter(cluster == "immune") %>%
  group_by(cluster) %>% top_n(n = 7, wt = avg_diff) %>%
  mutate(cluster = as.character(cluster)) %>%
  arrange(order(cluster),
          -avg_diff)
top6 <- Prot.markers %>% group_by(cluster) %>% 
  filter(cluster != "immune") %>%
  top_n(n = 6, wt = avg_diff) %>%
  ungroup() %>%
  filter(p_val_adj <= 0.001 | p_val_adj == 1) %>%
  mutate(cluster = as.character(cluster)) %>%
  filter(cluster != "unclassified") %>%
  arrange(order(cluster),
          -avg_diff)

extra <- Prot.markers  %>%
  filter(gene == "ERO1B" | gene == "TMEM176B")

top_all <- rbind(top6, top7_immune, extra)
top_all <- top_all[order(top_all$cluster,-top_all$avg_diff),]


my_levels <- c("alpha", "beta", "delta", "gamma", "immune", "unclassified",
               "ductal", "endothelial", "acinar", "activated-stellate"
               )

# Re-level object@ident
isletProt@active.ident <- factor(x = isletProt@active.ident, levels = my_levels)


library(viridis)
heat1 <- DoHeatmap(isletProt, assay  ="Prot", features = top_all$gene, slot = "scale.data",
                   group.bar.height = 0.05,
                   draw.lines = TRUE,
                   size = 4,
                   angle = 90,
                 #  hjust = 0.5,
                 #  vjust = 2,
                   label = F,
                   lines.width = 1
) + scale_fill_viridis(na.value = "white") + NoLegend()


heat1





#png(filename = "Azimuth_transfer_pancislets.png", width = 12.7,
#    height = 4, units = "in",
#    res = 300)
p3+p1
#dev.off()


#png(filename = "Protein_marker_Pancislets.png", width = 12.7,
#    height = 4, units = "in",
#    res = 300)
heat1
#dev.off()



