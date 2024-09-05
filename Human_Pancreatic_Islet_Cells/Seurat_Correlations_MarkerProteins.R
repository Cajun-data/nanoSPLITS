library(tidyverse)
library(Seurat)
library(Azimuth)
library(SeuratData)
library(limma)
################## Functions

MyFindAllMarkers <- function(
    object,
    assay = "Prot",
    slot = "data",
    logfc.threshold = 0,
    verbose = TRUE,
    ...
) {
  # Ensure the assay is available in the Seurat object
  DefaultAssay(object) <- assay
  
  # Extract data from the specified slot
  data <- GetAssayData(object, slot = slot)
  # check1 <- rowSums(is.na(data))
  # check1 <- check1[check1 < ncol(data)*min_obs]
  # data <- data[names(check1),]
  # Cluster information
  clusters <- Idents(object)
  
  # Initialize a data frame to store marker information
  all_markers <- data.frame()
  
  # Loop through each cluster and find markers
  for (cluster in levels(clusters)) {
    if (verbose) {
      message("Processing cluster: ", cluster)
    }
    
    # Cells in the current cluster
    cluster_cells <- WhichCells(object, idents = cluster)
    
    # Cells in other clusters
    other_cells <- setdiff(Cells(object), cluster_cells)
    
    # Create design matrix for limma
    design <- model.matrix(~ factor(c(rep(1, length(cluster_cells)), rep(0, length(other_cells)))))
    colnames(design) <- c("Other", "Cluster")
    
    # Combine expressions into a matrix
    expr_data <- cbind(data[, cluster_cells], data[, other_cells])
    
    cntrsts <- "Cluster-Other"
    
    # Fit linear model using limma
    fit <- lmFit(expr_data, design)
    made_contrasts <- makeContrasts(contrasts = cntrsts, levels = design)
    contrast_fit <- contrasts.fit(fit, made_contrasts)
    fit <- eBayes(fit)
    
    # Extract results
    results <- topTable(fit, coef = "Cluster", number = Inf, adjust.method = "BH")
    results <- results[!is.na(results$t), ]
    # Calculate log2 fold-change
    log2FC <- results$logFC
    
    # Create a data frame of results
    markers_df <- data.frame(
      gene = rownames(results),
      cluster = cluster,
      log2FC = log2FC,
      p_value = results$P.Value,
      p_value_adj = results$adj.P.Val
    )
    
    # Filter based on logfc.threshold and min.pct
    markers_df <- markers_df %>%
      filter(log2FC > logfc.threshold)
    
    # Append to the all_markers data frame
    all_markers <- rbind(all_markers, markers_df)
  }
  
  return(all_markers)
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


#png(filename = "Azimuth_transfer_pancislets.png", width = 12.7,
#    height = 4, units = "in",
#    res = 300)
p3+p1
#dev.off()

###########################
############################


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
  as.sparse()


isletProt <- CreateSeuratObject(combined3)
isletProt <- RenameAssays(
  isletProt,
  assay.name = "RNA",
  new.assay.name = "Prot",
  verbose = TRUE
)

isletProt <-  SetAssayData(isletProt,
                         layer = "data",
                         isletProt@assays$Prot$counts)


Seurat_scale <- combined %>%
  t() %>%
  scale() %>%
  t() %>%
  as.data.frame() %>%
  rownames_to_column(var = "Protein") %>%
  mutate(Protein = gsub("_","-", Protein)) %>%
  column_to_rownames(var ="Protein") %>%
  as.sparse()

isletProt <-  SetAssayData(isletProt,
                         layer = "scale.data",
                         Seurat_scale)



isletProt <- AddMetaData(
  object = isletProt,
  metadata = annot$Celltype,
  col.name = 'CellType'
)

Idents(object = isletProt, cells = 1:126) <- annot$Celltype

Prot.markers <- MyFindAllMarkers(isletProt, assay = "Prot", 
                                 slot = "data")

top_all <- c("GCG",
             "ALDH1A1",
             "GLS",
             "ELAPOR1",
             "ERP29",
             "SLC7A2",
             "INS",
             "HADH",
             "STX1A",
             "ERO1B",
             "GNAS",
             "SST",
             "RBP4",
             "PCSK1",
             "PCP4",
             "GRB2",
             "F13A1",
             "LSP1",
             "ACP5",
             "HLA-DRB3",
             "ALOX5AP"
             )

my_levels <- c("alpha", "beta", "delta", "gamma", "immune", "unclassified",
               "ductal", "endothelial", "acinar", "activated-stellate"
)

# Re-level object@ident
isletProt@active.ident <- factor(x = isletProt@active.ident, levels = my_levels)


library(viridis)
heat1 <- DoHeatmap(isletProt, assay  ="Prot", features = top_all, slot = "scale.data",
                   group.bar.height = 0.05,
                   draw.lines = T,
                   size = 4,
                   angle = 90,
                   #  hjust = 0.5,
                   #  vjust = 2,
                   label = T,
                   lines.width = 1
) + scale_fill_viridis(na.value = "black") + NoLegend()


#png(filename = "Protein_marker_Pancislets_updated.png", width = 12.7,
 #  height = 4, units = "in",
#   res = 300)
heat1
#dev.off()



x_long<- combined %>% as.data.frame() %>%
  rownames_to_column(var = "Gene") %>%
  pivot_longer(-Gene, names_to = "SampleID", values_to = "Intensity") %>%
  full_join(.,annot)


##############Identify markers based on missingvalues using hypergeometric test

markers_Missingvalue <- x_long %>%
  filter(!is.na(Intensity)) %>%
  group_by(Celltype, Gene) %>%
  add_tally(name = "GroupSize") %>%
  group_by(Celltype) %>%
  mutate(GroupSize = max(GroupSize)) %>%
  group_by(Gene) %>%
  add_count(name = "TotalN") %>%
  group_by(Celltype, Gene) %>%
  add_count(name = "GroupN") %>%
  filter(GroupSize == GroupN) %>%
  mutate(Check1 = TotalN - GroupSize) %>%
  ungroup() %>%
  mutate(Prob0 = (factorial(nrow(annot)-GroupSize))/(factorial(Check1)*factorial((nrow(annot)-GroupSize)-Check1))) %>%
  mutate(Prob1 = (factorial(nrow(annot)))/(factorial(TotalN)*factorial(nrow(annot)-TotalN))) %>%
  mutate(p_value = Prob0/Prob1) %>%
  mutate(Enrichment = GroupN/((GroupSize*TotalN)/nrow(annot))) %>%
  group_by(Celltype, Gene, p_value, Enrichment) %>%
  summarize(Avg = mean(Intensity))  %>%
  ungroup() %>%
  mutate(adj_p_value = p.adjust(p_value, method = "fdr")) %>%
  filter(p_value < 0.01) 

