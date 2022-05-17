library(tidyverse)


Prot <- read.delim("C10BufferOpt_Proteomics.txt",
                   check.names = FALSE)

Prot[Prot == 0] <- NA

x_long <- Prot %>%
  pivot_longer(!PROTID, names_to = "SampleID", values_to = "Intensity") %>%
  mutate(Intensity = log2(Intensity)) %>%
  mutate(Condition = case_when(grepl("\\+DDM_\\-R", SampleID) ~ "+DDM, -RNase",
                               grepl(paste("\\-DDM_\\-RNase"), SampleID) ~ "-DDM, -RNase",
                               grepl("\\+DDM_\\+R", SampleID) ~ "+DDM, +RNase",
                               grepl("\\-DDM_\\+RNase", SampleID) ~ "-DDM, +RNase"))


x_long %>%
  filter(!is.na(Intensity)) %>%
  distinct(SampleID,PROTID, Condition) %>%
  group_by(SampleID) %>%
  add_count(name = "n") %>%
  ungroup() %>%
  group_by(Condition) %>%
  summarise_at("n", funs(mean, sd)) %>%
  ggplot()+
  aes(x = Condition, y = mean, fill = "white")+
  geom_bar(stat= "identity",show.legend = FALSE, alpha =  0.5)+
  theme_bw(base_size = 20) +
  theme(legend.position = "none",
        panel.background = element_rect(fill= 'white'),
        axis.text.y=element_text(color = 'black', size = 20),
        axis.text.x=element_text(angle = 60, vjust = 0.5, hjust= 0.5, color='black',size = 20),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(),
        text=element_text(family="Helvetica")) +
  ylab("Total Proteins (n)")+
  scale_x_discrete(labels=function(x){sub("\\s", "\n", x)}) +
  geom_errorbar( aes(x=Condition,
                     ymin=mean-sd, ymax=mean+sd), width=0.4, colour="black", alpha=1, size=0.8) +
  xlab(" ")+
  scale_fill_manual(values = c("#80b1d3"))



RNA <- read.delim(file = "C10BufferOpt_Transcriptomics.txt",
                  stringsAsFactors = FALSE, check.names = FALSE )



RNA <- RNA %>%
  pivot_longer(!Gene, names_to = "SampleID", values_to = "RPKM") %>%
  mutate(Condition = case_when(grepl("\\+DDM_\\-R", SampleID) ~ "+DDM, -RNase",
                               grepl(paste("\\-DDM_\\-RNase"), SampleID) ~ "-DDM, -RNase",
                               grepl("\\+DDM_\\+R", SampleID) ~ "+DDM, +RNase",
                               grepl("\\-DDM_\\+RNase", SampleID) ~ "-DDM, +RNase"))


RNA %>%
  filter(RPKM > 0) %>%
  distinct(SampleID,Gene, Condition) %>%
  group_by(SampleID) %>%
  add_count(name = "n") %>%
  ungroup() %>%
  group_by(Condition) %>%
  summarise_at("n", funs(mean, sd)) %>%
  ggplot()+
  aes(x = Condition, y = mean, fill = "white")+
  geom_bar(stat= "identity",show.legend = FALSE, alpha =  0.5)+
  theme_bw(base_size = 20) +
  theme(legend.position = "none",
        panel.background = element_rect(fill= 'white'),
        axis.text.y=element_text(color = 'black', size = 20),
        axis.text.x=element_text(angle = 60, vjust = 0.5, hjust= 0.5, color='black',size = 20),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(),
        text=element_text(family="Helvetica")) +
  ylab("Total Genes (n)")+
  scale_x_discrete(labels=function(x){sub("\\s", "\n", x)}) +
  geom_errorbar( aes(x=Condition,
                     ymin=mean-sd, ymax=mean+sd), width=0.4, colour="black", alpha=1, size=0.8) +
  xlab(" ")+
  scale_fill_manual(values = c("#d95f0e"))
