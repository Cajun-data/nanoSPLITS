library(tidyverse)
library(proDA)
library(ggExtra)
library(viridis)
library(ggpointdensity)


#### Proteomics of both donor and acceptor chip, nanoSPLITS
#########Script for determining protein recovery from pooled C10 cells


######https://slowkow.com/notes/ggplot2-color-by-density/
get_density <- function(x, y, ...) {
  dens <- MASS::kde2d(x, y, ...)
  ix <- findInterval(x, dens$x)
  iy <- findInterval(y, dens$y)
  ii <- cbind(ix, iy)
  return(dens$z[ii])
}


singcell <- read_tsv("Protein_Recovery.tsv", show_col_types = FALSE)
prot_mass <- read.delim("protein_MW_list.tsv")

singcell2 <- singcell %>%
  filter(!grepl("contam_sp", PROTID)) %>%
  filter(PROTID %in% prot_mass$PROTID)

singcell2[singcell2 == 0] <- NA


singcell2[,2:13] <- log2(singcell2[,2:13])
singcell_long <- singcell2 %>%
  pivot_longer(!PROTID, names_to = "SampleID", values_to = "Intensity") %>%
  filter(!is.na(Intensity)) %>%
  mutate(Type = case_when(grepl("10cell", SampleID) ~ "10cell",
                          grepl("NoCell", SampleID) ~ "0cell"))

wide1 <- singcell_long %>%
  filter(Type == "10cell") %>%
  pivot_wider(names_from = "SampleID", values_from = Intensity) %>%
  select(-Type) %>%
  column_to_rownames(var = "PROTID")

long1 <- as.data.frame(median_normalization(as.matrix(wide1))) %>%
  rownames_to_column(var = "PROTID") %>%
  pivot_longer(!PROTID, names_to = "SampleID", values_to = "Intensity")

wide2 <- singcell_long %>%
  filter(Type == "0cell") %>%
  pivot_wider(names_from = "SampleID", values_from = Intensity) %>%
  select(-Type) %>%
  column_to_rownames(var = "PROTID")

long2 <- as.data.frame(median_normalization(as.matrix(wide2))) %>%
  rownames_to_column(var = "PROTID") %>%
  pivot_longer(!PROTID, names_to = "SampleID", values_to = "Intensity")

singcell_long1 <- full_join(long1, long2)


singcell_long1 %>%
  mutate(Type = case_when(grepl("10cell", SampleID) ~ "10cell",
                          grepl("NoCell", SampleID) ~ "0cell")) %>%
  ggplot()+
  aes(x = Type, y = Intensity)+
  geom_violin()



 figure_1 <- singcell_long1 %>%
  mutate(Type = case_when(grepl("10cell", SampleID) ~ "10cell",
                          grepl("NoCell", SampleID) ~ "0cell")) %>%
  group_by(PROTID, Type) %>%
  mutate(Intensity = 2^Intensity) %>%
  mutate(Avg = mean(Intensity)) %>%
  mutate(CV = sd(Intensity)/(Avg)*100) %>%
  distinct(PROTID, Avg, CV, Type) %>%
  filter(!is.na(CV)) %>%
  ungroup() %>%
  pivot_wider( names_from = Type, values_from = c(Avg, CV)) %>%
  mutate(Prop = ((Avg_10cell)/(Avg_10cell+Avg_0cell)*100)) %>%
   filter(!is.na(Prop)) %>%
  mutate(Avg_10cell = log2(Avg_10cell))

 fig1_1 <- figure_1 %>%
   ggplot()+
   aes( y = Prop, fill = "")+
   geom_boxplot(alpha = 0.5)+
   ylab("Relative Proportion Retained on Donor Chip (%)")+
   xlab("")+
   theme_bw(base_size = 16) +
   theme(legend.position = "none",
         panel.background = element_rect(fill= 'white'),
         axis.text.y=element_text(color = 'black', size = 14),
         axis.text.x=element_blank(),
         axis.ticks.x=element_blank(),
         axis.title.x=element_blank(),
         panel.grid.minor = element_blank(),
         panel.grid.major = element_blank(),
         panel.border = element_blank(),
         axis.line = element_line(),
         text=element_text(family="Helvetica"))+
   scale_fill_manual(values = c("#28819d" ))

  
  
figure_2 <- figure_1 %>%
  inner_join(., prot_mass)

  
  figure_2$density <- get_density(figure_2$`Prop`,
                                  figure_2$`Mass`, n = 2000)
  p1 <- figure_2 %>%
    ggplot()+
    aes(x = Mass, y = Prop, 
        )+
    geom_pointdensity(size = 4,
                      adjust = 10)+
    scale_color_viridis()+
    theme_minimal(base_size = 16)+
    geom_hline(yintercept = 50,
               linetype='dotted', col = 'red')+
    xlab("Protein Mass (kDa)")+
    ylab("Relative Proportion Retained (%)")+
    theme(legend.position = "none")+
    scale_x_continuous(limits = c(0,150000))
  
p1 <- p1 %>%
    ggMarginal(type = "histogram", fill = "white")

pdf("prop_vs_mass.pdf", height = 6, width = 6)
p1
dev.off()



xmed <- singcell_long1 %>%
  mutate(Type = case_when(grepl("10cell", SampleID) ~ "10cell",
                          grepl("NoCell", SampleID) ~ "0cell")) %>%
  filter(!is.na(Intensity)) %>%
  mutate(Intensity = 2^Intensity) %>%
  group_by(Type, PROTID) %>%
  add_count(name = "Obs") %>%
  filter(Obs >= 2) %>%
  mutate(CV = sd(Intensity)/mean(Intensity)) %>%
  distinct(Type, PROTID,CV) %>%
  group_by(Type) %>%
  mutate(med = median(CV)) %>%
  distinct(Type, med) %>%
  ungroup()

 
  
main1 <- singcell_long1 %>%
    mutate(Type = case_when(grepl("10cell", SampleID) ~ "10cell",
                            grepl("NoCell", SampleID) ~ "0cell")) %>%
    group_by(PROTID, Type) %>%
  add_count(name = "Obs") %>%
  filter(Obs >= 2) %>%
    mutate(Intensity = 2^Intensity) %>%
    mutate(Avg = mean(Intensity)) %>%
    mutate(CV = sd(Intensity)/(Avg)) %>%
    distinct(PROTID, Avg, CV, Type) %>%
  ggplot()+
  aes(x = Type, y = CV, fill = Type)+
  geom_violin(alpha = 0.5, show.legend = FALSE, size = 0)+
  scale_y_continuous(limits = c(0,1.6))+
  theme_bw(base_size = 16) +
  theme(legend.position = "none",
        panel.background = element_rect(fill= 'white'),
        axis.text.y=element_text(color = 'black', size = 14),
        axis.text.x=element_text(angle = 0, vjust = 0.5, hjust= 0.5, color='black',size = 14),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(),
        text=element_text(family="Helvetica")) +
  ylab("Coefficent of Variation (CV)") +
  xlab("Number of Cells")+
  scale_x_discrete(labels=function(x){sub("\\s", "\n", x)}) +
  geom_text(data = xmed, aes(x = Type, y = med, label = paste(round(med, digits = 2))), 
            size = 5, vjust = -2.2,hjust = -0.25, show.legend = FALSE) +
  scale_fill_manual(values = c("#88419d", "#8c96c6"))

data_summary <- function(x) {
  m <- median(x)
  ymin <- m-sd(x)
  ymax <- m+sd(x)
  return(c(y=m,ymin=ymin,ymax=ymax))
}

pdf("prop_recov_cv.pdf", height = 6, width = 3)
main1  +
  stat_summary(fun.data=data_summary)
dev.off()


pdf("proportion.pdf", height = 6, width = 3)
fig1_1
dev.off()





  