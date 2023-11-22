library(mixOmics)
library(tidyverse)


####https://rdrr.io/github/jdreyf/jdcbioinfo/man/impute_normal.html
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

################################################

panc_pept <- read.delim("combined_modified_peptide_MaxLFQ.tsv") %>%
  filter(!grepl("contam", Protein)) %>%
  filter(grepl("HUMAN", Protein)) %>%
  filter(!grepl("Keratin", Protein.Description)) 

annot <- read.csv("annotations.csv",
                  check.names = T)


annot2 <- annot %>%
  mutate(Celltype = case_when(Celltype == "endothelial" ~ "Other",
                              Celltype == "ductal" ~ "Other",
                              Celltype == "acinar" ~ "Other",
                              Celltype == "activated-stellate" ~ "Other",
                              Celltype == "immune" ~ "Other",
                              TRUE ~ Celltype)) %>%
  group_by(Celltype) %>%
  filter(AnnotationScore >= 0.8) %>%
  add_count(name = "nObs") %>%
  ungroup() %>%
  filter(!is.na(Celltype)) %>%
  arrange(SampleID) %>%
  filter(nObs > 2)

annot_inv <- annot %>%
  filter(!SampleID %in% annot2$SampleID) %>%
  filter(AnnotationScore == 1)

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
  .[[1]] %>%
  impute_normal(width=0.3, downshift=3) %>%
  t()
unknowns <- combined[!row.names(combined)%in%annot2$SampleID,]
unknowns <- unknowns[!row.names(unknowns)%in%annot_inv$SampleID,]
unknowns <- unknowns[order(row.names(x = unknowns)),]

knowns <- combined[row.names(combined)%in%annot2$SampleID,]
knowns <- knowns[order(row.names(x = knowns)),]


color <- annot2 %>%
  filter(SampleID %in% rownames(knowns)) %>%
  arrange(SampleID)

pca_j <- PCAtools::pca(t(knowns))

PCAtools::biplot(pca_j,
                 x = "PC1",
                 y = "PC2",
                 showLoadings = F,
                 ntopLoadings = 5,
                 lab = color$Celltype
                 # colkey = color$Celltype
)

annot3 <- color %>%
  mutate(Celltype = as.factor(Celltype)) %>%
  dplyr::select(Celltype)


pca.srbct = pca(knowns, ncomp = 8, center = T, scale = T) 
plot(pca.srbct)
lda1 <- splsda(knowns, annot3$Celltype, ncomp = 8,
               scale = T)


plotIndiv(lda1 , comp = 1:2, 
          group = annot3$Celltype, ind.names = FALSE,  # colour points by class
          ellipse = TRUE, # include 95% confidence ellipse for each class
          legend = TRUE, title = '(a) PLSDA with confidence ellipses')



perf.splsda.srbct <- perf(lda1, validation = "Mfold", ,
                          folds = 5, nrepeat = 10, # use repeated cross-validation
                          progressBar = T, auc = TRUE) # include AUC values


plot(perf.splsda.srbct, col = color.mixo(5:7), sd = TRUE,
     legend.position = "horizontal")
perf.splsda.srbct$choice.ncomp


list.keepX <- c(1:20)


# undergo the tuning process to determine the optimal number of variables
tune.splsda.srbct <- tune.splsda(knowns, annot3$Celltype, ncomp = 7,
                                 validation = 'Mfold',
                                 scale =T,
                                 folds = 5, nrepeat = 10, # use repeated cross-validation
                                 dist = 'max.dist', 
                                 measure = "BER", # use balanced error rate of dist measure
                                 test.keepX = list.keepX,
                                 cpus = 3, progressBar = T) 

plot(tune.splsda.srbct, col = color.jet(7)) # plot output of variable number tuning


optimal.ncomp <- tune.splsda.srbct$choice.ncomp$ncomp
optimal.keepX <- tune.splsda.srbct$choice.keepX[1:optimal.ncomp]




final.splsda <- splsda(knowns,
                       annot3$Celltype, 
                       ncomp = optimal.ncomp, 
                       keepX = optimal.keepX,
                       scale =T)



predict.splsda.srbct <- predict(final.splsda, unknowns, 
                                dist = "max.dist")
predict.comp2 <- predict.splsda.srbct$class$max.dist[,3]
table(factor(predict.comp2))

celltypes <- as.data.frame(predict.comp2) 
celltypes2 <- celltypes %>%
  rename(Celltype = predict.comp2) %>%
  rownames_to_column(var = "SampleID") %>%
  mutate(AnnotationScore = 1.1) %>%
  rbind(., annot) %>%
  filter(SampleID %in% batch_meta$SampleID) %>%
  arrange(SampleID) %>%
  group_by(SampleID) %>%
  slice_max(., order_by = AnnotationScore, n = 1) %>%
  mutate(Celltype = case_when(Celltype == "Other" ~ "unclassified",
                              TRUE ~ Celltype)) %>%
  ungroup()

#### write.csv(celltypes2, file = "LDA_Celltypes.csv")


combined <- combined[order(row.names(x = combined)),]

color_post <- celltypes2 %>%
  filter(SampleID %in% rownames(combined)) %>%
  arrange(SampleID)

color_pre <- celltypes2 %>%
  filter(SampleID %in% rownames(combined)) %>%
  arrange(SampleID) %>%
  mutate(Annotation_Source = case_when(Celltype == "unclassified" ~ "unclassified",
                              AnnotationScore == 1.1 ~ "spls-DA",
                              TRUE ~ "Azimuth")) %>%
  column_to_rownames(var = "SampleID")

pca_j <- PCAtools::pca(t(combined),
                       metadata = color_pre)

PCAtools::biplot(pca_j,
                 x = "PC2",
                 y = "PC3",
                 showLoadings = F,
                 ntopLoadings = 5,
                 legendPosition = "right",
                 colby = 'Celltype',
                 drawConnectors = F,
                 lab = NULL
)


PCAtools::biplot(pca_j,
                 x = "PC2",
                 y = "PC3",
                 showLoadings = F,
                 ntopLoadings = 5,
                 legendPosition = "right",
                 colby = 'Annotation_Source',
                 drawConnectors = F,
                 lab = NULL,
)


