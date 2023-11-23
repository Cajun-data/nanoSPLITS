library(tidyverse)
library(EnhancedVolcano)

load("Proteomics_results.RData")

pp2 <- clu_df %>%
  filter(Cluster == "cluster2") %>%
  select(Gene)


pp3 <- clu_df %>%
  filter(Cluster == "cluster3") %>%
  select(Gene)

pp4 <- clu_df %>%
  filter(Cluster == "cluster4") %>%
  select(Gene)


label <- c("Vnn1", "Aldh3a1", "Acot13", "Serpinb6")

man1 <- c("Gsta1", "Gsta3", "Gsta4",
          "Gstk1", "Gsta2",
          "Ctsa", "Ctsb", "Ctsd")


pp6 <- do.call("rbind", list(pp2,pp3, pp4))


cs_x <- data_results2 %>%
  mutate(Protein = Gene) %>%
  mutate(adj.pvalue = Treated_vs_Control_p.adj) %>%
  mutate(log2FC = Treated_vs_Control_diff) %>%
  dplyr::select(Protein, adj.pvalue, log2FC) %>%
  group_by(Protein) %>%
  slice_min(adj.pvalue, n = 1) %>%
  mutate(Color = case_when(Protein %in% pp3$Gene ~ "Clust3",
                           Protein %in% pp2$Gene ~ "Clust2",
                           Protein %in% pp4$Gene ~ "Clust4",
                           TRUE ~ "NaN")) %>%
  mutate(Color = case_when(Color == "NaN" & adj.pvalue > 0.01 ~ "grey",
                           Color == "NaN" & between(log2FC, -0.5, 0.5)~ "grey",
                           Color == "NaN" & adj.pvalue <= 0.01 ~ "royalblue1",
                           TRUE ~ Color)) %>%
  arrange(desc(Color)) %>%
  ungroup()


rownames(cs_x) <- cs_x$Protein


keyvals <- case_when(grepl("Clust3", cs_x$Color) ~ "#984ea3",
                     grepl("Clust2", cs_x$Color) ~ "#e41a1c",
                     grepl("Clust4", cs_x$Color) ~ "#4daf4a",
                     grepl("grey", cs_x$Color) ~ "grey",
                     grepl("royalblue1", cs_x$Color) ~ "royalblue1" )

keyvals[is.na(keyvals)] <- 'NA'
names(keyvals)[keyvals == '#4daf4a'] <- '4'
names(keyvals)[keyvals == '#e41a1c'] <- '2'
names(keyvals)[keyvals == '#984ea3'] <- '3'
names(keyvals)[keyvals == 'royalblue1'] <- 'p-value'
names(keyvals)[keyvals == 'grey'] <- 'NS'
names(keyvals)[keyvals == 'NA'] <- ' '



pdf('volcano_plot_scProteomics.pdf', width = 3, height = 9)
EnhancedVolcano(cs_x,
                lab = cs_x$Protein,
                labSize = 5,
                drawConnectors = TRUE,
                boxedLabels = TRUE,
                arrowheads = F,
                selectLab = ifelse( cs_x$Protein %in% pp3$Gene & cs_x$log2FC > 1.2 | cs_x$Protein %in% man1 | cs_x$Protein %in% pp6$Gene & cs_x$log2FC < -1.5 , cs_x$Protein, ""),
                x = 'log2FC',
                y = 'adj.pvalue',
                pCutoff = 0.01,
                FCcutoff = 0.5,
                axisLabSize = 14,
                xlab = NULL,
                ylab = NULL,
                xlim = c(-3,3),
                ylim = c(0,25),
                colCustom = keyvals,
                colAlpha = c(ifelse(cs_x$Color == "royalblue1" | cs_x$Color == "grey", 0.3, 1)),
                pointSize = c(ifelse(cs_x$Color != "grey", 3, 3)),
                title = 'scProteomics',
                legendPosition = 'top',
                legendLabSize = 12,
                #shapeCustom = keyvals.shape,
                legendIconSize = 4.0)+
  theme(panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        axis.line = element_line())
dev.off()


load("scRNAseq_DE.RData")


cs_x <- results2 %>%
  filter(logpvalue != 0) %>%
  mutate(adjpvalue = 10^-logpvalue) %>%
  distinct(Gene, log2FC,  .keep_all = TRUE) %>%
  mutate(Color = case_when(Gene %in% pp2$Gene ~ "Clust2",
                           Gene %in% pp3$Gene ~ "Clust3",
                           Gene %in% pp4$Gene ~ "Clust4",
                           TRUE ~ "NaN")) %>%
  mutate(Color = case_when(Color == "NaN" & adjpvalue > 0.01 ~ "grey",
                           Color == "NaN" & between(log2FC, -0.5, 0.5)~ "grey",
                           Color == "NaN" & adjpvalue <= 0.01 ~ "royalblue1",
                           TRUE ~ Color)) %>%
 arrange(desc(Color))


rownames(cs_x) <- cs_x[,1]



keyvals <- case_when(grepl("Clust3", cs_x$Color) ~ "#984ea3",
                     grepl("Clust2", cs_x$Color) ~ "#e41a1c",
                     grepl("Clust4", cs_x$Color) ~ "#4daf4a",
                     grepl("grey", cs_x$Color) ~ "grey",
                     grepl("royalblue1", cs_x$Color) ~ "royalblue1" )

keyvals[is.na(keyvals)] <- 'NA'
names(keyvals)[keyvals == '#4daf4a'] <- '4'
names(keyvals)[keyvals == '#e41a1c'] <- '2'
names(keyvals)[keyvals == '#984ea3'] <- '3'
names(keyvals)[keyvals == 'royalblue1'] <- 'p-value'
names(keyvals)[keyvals == 'grey'] <- 'NS'
names(keyvals)[keyvals == 'NA'] <- ' '

pdf('volcano_plot_scRNAseq.pdf', width = 3, height = 9)
EnhancedVolcano(cs_x,
                lab = cs_x$Gene,
                x = 'log2FC',
                y = 'adjpvalue',
                labSize = 5,
                drawConnectors = TRUE,
                arrowheads = F,
                boxedLabels = TRUE,
                selectLab = ifelse(cs_x$Gene %in% pp2$Gene & cs_x$log2FC < -1 | cs_x$Gene %in% man1 & cs_x$log2FC > 2, cs_x$Gene, ""),
                pCutoff = 0.01,
                FCcutoff = 0.5,
                axisLabSize = 14,
                xlab = NULL,
                ylab = NULL,
                xlim = c(-4,9),
                ylim = c(0,15),
                colCustom = keyvals,
                colAlpha = c(ifelse(cs_x$Color == "royalblue1" | cs_x$Color == "grey", 0.3, 1)),
                pointSize = c(ifelse(cs_x$Color != "grey", 3, 3)),
                title = 'scRNAseq',
                legendPosition = 'top',
                legendLabSize = 12,
                #shapeCustom = keyvals.shape,
                legendIconSize = 4.0)+
  theme(panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        axis.line = element_line())
#  scale_y_reverse(limits = c(15,0))+
 # scale_x_continuous(position = "top", limits = c(-3,3)) 
dev.off()



load("Phosphopeptide_DE.RData")


cs_x <- peptide_results %>%
  mutate(Protein = Gene) %>%
  mutate(adj.pvalue = Treated_vs_Control_p.adj) %>%
  mutate(log2FC = Treated_vs_Control_diff) %>%
  dplyr::select(Protein, Modified.Sequence, adj.pvalue, log2FC) %>%
  group_by(Modified.Sequence) %>%
  slice_min(adj.pvalue, n = 1) %>%
  mutate(Color = case_when(grepl("79.9", Modified.Sequence) ~ "ab",
                           adj.pvalue > 0.01 ~ "grey",
                           between(log2FC, -0.5, 0.5)~ "grey",
                           adj.pvalue <= 0.01 ~ "royalblue1",
                          adj.pvalue <= 0.01 ~ "royalblue1"
                          )) %>%
  mutate(Protein = case_when(Color == "ab" ~ paste(" ",Protein, sep = ""),
                             TRUE ~ Protein)) %>%
  arrange(desc(Color)) %>%
  ungroup()

labelp <- c(" Top2a", " Vim", " Hnrnpu", " Npm1", " Hmgn1")
  
rownames(cs_x) <- cs_x$Modified.Sequence


keyvals <- case_when(grepl("ab", cs_x$Color) ~ "orange",
                     grepl("grey", cs_x$Color) ~ "grey",
                     grepl("royalblue1", cs_x$Color) ~ "royalblue1" )

keyvals[is.na(keyvals)] <- 'NA'
names(keyvals)[keyvals == 'orange'] <- 'Phosphorylated'
names(keyvals)[keyvals == 'royalblue1'] <- 'p-value'
names(keyvals)[keyvals == 'grey'] <- 'NS'
names(keyvals)[keyvals == 'NA'] <- ' '

pdf('volcano_plot_phosphopeptide.pdf', width = 3, height = 9)
EnhancedVolcano(cs_x,
                lab = cs_x$Protein,
                labSize = 5,
                drawConnectors = TRUE,
                boxedLabels = TRUE,
                arrowheads = F,
                selectLab = ifelse(cs_x$Protein %in% labelp, cs_x$Protein, ""),
                x = 'log2FC',
                y = 'adj.pvalue',
                pCutoff = 0.01,
                axisLabSize = 14,
                xlab = NULL,
                ylab = NULL,
                FCcutoff = 0.5,
                xlim = c(-3,3),
                ylim = c(0,30),
                colCustom = keyvals,
                colAlpha = c(ifelse(cs_x$Color == "royalblue1" | cs_x$Color == "grey", 0.25, 1)),
                pointSize = c(ifelse(cs_x$Color != "grey", 3, 3)),
                title = 'scProteomics',
                legendPosition = 'top',
                legendLabSize = 12,
                legendIconSize = 4.0)+
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line())
dev.off()

