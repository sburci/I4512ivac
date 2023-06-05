# I 4,[5],12:i:- MICROBIOME ANALYSIS
# Selma Burciaga
# March 20, 2023
# Rscript: I4512ivac_nmds.R
# R v4.2.2/R v2022.12.0

# PART 1: Using the vegan R package to calculate ecological distances (CC188) - https://riffomonas.org/code_club/2022-02-17-distances
## anything ending with f is for Fecal samples
## anything ending with f2 is for Fecal samples without MVMC treatment
## anything ending with cc is for Cecal Content samples


# BETA DIVERSITY: NMDS___________________________________________________
library(phyloseq) #v1.42.0
library(tidyverse) #v1.3.2
library(vegan) #v2.6.4
library(ggplot2) #3.4.0
library(dplyr)
library(readr)
library(reshape2)

## Open I4512ivac_phylo_all.rds
PHYLO <- read_rds("./outputs/I4512ivac_phylo_all.rds")

## Set reference levels
sample_data(PHYLO)$Vaccine <- factor(sample_data(PHYLO)$Vaccine, levels=c('Mock', 'Enterisol'))
sample_data(PHYLO)$Challenge <- factor(sample_data(PHYLO)$Challenge, levels=c('Mock', 'I4512i'))

## Create subsets (Fecal and Cecal_contents)
PHYLOf <- prune_samples(samples = PHYLO@sam_data$Sample_Type=="Fecal", x = PHYLO)
PHYLOf2 <- prune_samples(samples = PHYLOf@sam_data$Challenge=="I4512i", x = PHYLOf)
PHYLOcc <- prune_samples(samples = PHYLO@sam_data$Sample_Type=="Cecal_contents", x = PHYLO)

## Create OTU R matrix and dataframe with metadata from phyloseq object
OTUmatf <- as(otu_table(PHYLOf), "matrix")
METf <- sample_data(PHYLOf) %>% 
  as(Class="data.frame")
OTUmatf2 <- as(otu_table(PHYLOf2), "matrix")
METf2 <- sample_data(PHYLOf2) %>% 
  as(Class="data.frame")
OTUmatcc <- as(otu_table(PHYLOcc), "matrix")
METcc <- sample_data(PHYLOcc) %>% 
  as(Class="data.frame")


## Rarefying distance calculations using avgdist (recommended) 
### If file doesn't exist, calculate dist and save, if it does read file and save as DIST (sample= min otu counts)
set.seed(19760620) #usually before func call, best for reproducibility

if (!file.exists("outputs/I4512ivac_fecal_dist.rds")) {
  DISTf <- avgdist(OTUmatf, dmethod="bray", sample=min(rowSums(OTUmatf))) 
  write_rds(DISTf, "outputs/I4512ivac_fecal_dist.rds")
} else {
  DISTf <- readRDS("outputs/I4512ivac_fecal_dist.rds")
}

if (!file.exists("outputs/I4512ivac_fecal_dist2.rds")) {
  DISTf2 <- avgdist(OTUmatf2, dmethod="bray", sample=min(rowSums(OTUmatf2))) 
  write_rds(DISTf2, "outputs/I4512ivac_fecal_dist2.rds")
} else {
  DISTf2 <- readRDS("outputs/I4512ivac_fecal_dist2.rds")
}

if (!file.exists("outputs/I4512ivac_cecalcontents_dist.rds")) {
  DISTcc <- avgdist(OTUmatcc, dmethod="bray", sample=min(rowSums(OTUmatcc)))
  write_rds(DISTcc, "outputs/I4512ivac_cecalcontents_dist.rds")
} else {
  DISTcc <- readRDS("outputs/I4512ivac_cecalcontents_dist.rds")
}

## Non-metric Multidimensional Scaling (NMDS) - create NMDS points and adhere to MET objects
set.seed(1)
NMDSf <- metaMDS(DISTf)
NMDSscoresf <- scores(NMDSf) %>% 
  as_tibble(rownames="SampleID")
METf <- METf %>% 
  left_join(NMDSscoresf)

set.seed(1)
NMDSf2 <- metaMDS(DISTf2)
NMDSscoresf2 <- scores(NMDSf2) %>% 
  as_tibble(rownames="SampleID")
METf2 <- METf2 %>% 
  left_join(NMDSscoresf2)

set.seed(1)
NMDScc <- metaMDS(DISTcc)
NMDSscorescc <- scores(NMDScc) %>% 
  as_tibble(rownames="SampleID")
METcc <- METcc %>% 
  left_join(NMDSscorescc)

    ###Calculation of ordination stress
NMDSf$stress  # = 0.1434989
NMDScc$stress # = 0.1324859

   ### View order of variables in METf columns
METf$Days_post_challenge %>% 
  unique()
METf$Treatment %>% 
  unique()

    ### Describing centroids: Group_by dpc and Treatment
METf <- METf %>% 
  group_by(dpc, Treatment) %>% 
  dplyr::mutate(centnmds1=mean(NMDS1),
         centnmds2=mean(NMDS2)) %>% 
  ungroup()
METf2 <- METf2 %>% 
  group_by(dpc, Treatment) %>% 
  dplyr::mutate(centnmds1=mean(NMDS1),
         centnmds2=mean(NMDS2)) %>% 
  ungroup()
METcc <- METcc %>% 
  group_by(dpc, Treatment) %>% 
  dplyr::mutate(centnmds1=mean(NMDS1),
         centnmds2=mean(NMDS2)) %>% 
  ungroup()

## Bringing metadata into ordination after rarefaction (e.g. all treatment groups individual graphs for each day post Salmonella challenge)
#Fecal
###By days post challenge
tiff(filename="./graphics/I4512ivac_fecal_alltreatmentsbydpc_nmds.tiff", width=8, height=8, units="in", res=600, bg="transparent", family = "sans", pointsize = 12,  compression="lzw")
METf %>% 
  ggplot(aes(x=NMDS1, y=NMDS2, color=Treatment)) +
  geom_point() +
  geom_segment(aes(xend=centnmds1, yend=centnmds2)) +
  ggtitle('NMDS of all treatment groups over time by days post challenge') +
  scale_color_manual(values=c(MVMC='green', MVC='skyblue', EVC='orange')) +
  facet_wrap(~Days_post_challenge, scales = "free") +
  theme_bw()
dev.off()
    ###Excluding MVMC (which is usually not included in comparisons of vaccination/challenged groups)
sub1 <- subset(METf, Treatment %in% c("MVC", "EVC"))
tiff(filename="./graphics/I4512ivac_fecal_alltreatmentsbydpc_nomvmc_nmds.tiff", width=8, height=8, units="in", res=600, bg="transparent", family = "sans", pointsize = 12,  compression="lzw")
sub1 %>% 
  ggplot(aes(x=NMDS1, y=NMDS2, color=Treatment)) +
  geom_point() +
  geom_segment(aes(xend=centnmds1, yend=centnmds2)) +
  ggtitle('NMDS of all treatment groups over time by days post challenge') +
  scale_color_manual(values=c(MVC='skyblue', EVC='orange')) +
  facet_wrap(~Days_post_challenge, scales = "free") +
  theme_bw()
dev.off()

###NMDS with centroids
tiff(filename="./graphics/I4512ivac_fecal_alltreatments_nmds_centroids.tiff", width=8, height=8, units="in", res=600, bg="transparent", family = "sans", pointsize = 12,  compression="lzw")
METf %>% 
  dplyr::select(centnmds1, centnmds2, Treatment, dpc) %>% 
  unique() %>% 
  arrange(dpc) %>% 
  ggplot(aes(x=centnmds1, y=centnmds2, color=Treatment)) +
  geom_point(data=METf, aes(x=NMDS1, y=NMDS2, color=Treatment), alpha=.5) +
  geom_segment(data=METf, aes(x=NMDS1, y=NMDS2, xend=centnmds1, yend=centnmds2, color=Treatment), alpha=.25) +
  geom_path(aes(group=Treatment)) +
  geom_point(aes(shape = Treatment), size=6) +
  geom_text(aes(label=dpc), color='black') +
  ggtitle('NMDS of all treatment groups over time',
          'Fecal samples') +
  xlab('NMDS1') +
  ylab('NMDS2') +
  scale_color_manual(values=c(MVMC='green', MVC='skyblue', EVC='orange')) +
  annotate("text", x=0.5, y=-0.30, label="Stress = 0.144", size = 4) +
  theme_bw() +
  theme(legend.position="bottom")
dev.off()
    ###Excluding MVMC
tiff(filename="./graphics/I4512ivac_fecal_alltreatments_nomvmc_nmds_centroids.tiff", width=8, height=8, units="in", res=600, bg="transparent", family = "sans", pointsize = 12,  compression="lzw")
sub1 %>% 
  dplyr::select(centnmds1, centnmds2, Treatment, dpc) %>% 
  unique() %>% 
  arrange(dpc) %>% 
  ggplot(aes(x=centnmds1, y=centnmds2, color=Treatment)) +
  geom_point(data=sub1, aes(x=NMDS1, y=NMDS2, color=Treatment), alpha=.5) +
  geom_segment(data=sub1, aes(x=NMDS1, y=NMDS2, xend=centnmds1, yend=centnmds2, color=Treatment), alpha=.25) +
  geom_path(aes(group=Treatment)) +
  geom_point(size=6) +
  geom_text(aes(label=dpc), color='black') +
  ggtitle('NMDS of treatment groups over time') +
  xlab('NMDS1') +
  ylab('NMDS2') +
  scale_color_manual(values=c(MVMC='green', MVC='skyblue', EVC='orange')) +
  theme_bw()
dev.off()

#Cecal contents
tiff(filename="./graphics/I4512ivac_cecalcontents_alltreatmentsbydpc_nmds.tiff", width=8, height=8, units="in", res=600, bg="transparent", family = "sans", pointsize = 12,  compression="lzw")
METcc %>% 
  ggplot(aes(x=NMDS1, y=NMDS2, color=Treatment)) +
  geom_point() +
  geom_segment(aes(xend=centnmds1, yend=centnmds2)) +
  ggtitle('NMDS of cecal contents 14 days post challenge') +
  scale_color_manual(values=c(MVMC='green', MVC='skyblue', EVC='orange')) +
  facet_wrap(~Days_post_challenge, scales = "free") +
  theme_bw()
dev.off()

##By days post challenge
    ###No MVMC Fecal
tiff(filename="./graphics/I4512ivac_fecal_alltreatmentsbydpc_nomvmc_nmds2_centroids.tiff", width=8, height=8, units="in", res=600, bg="transparent", family = "sans", pointsize = 12,  compression="lzw")
METf2 %>% 
  ggplot(aes(x=NMDS1, y=NMDS2, color=Treatment)) +
  geom_point() +
  geom_segment(aes(xend=centnmds1, yend=centnmds2)) +
  ggtitle('NMDS of all treatment groups over time by days post challenge') +
  scale_color_manual(values=c(MVC='skyblue', EVC='orange')) +
  facet_wrap(~Days_post_challenge, scales = "free") +
  theme_bw()
dev.off()
    ###NMDS with centroids
tiff(filename="./graphics/I4512ivac_fecal_alltreatments_nomvmc_nmds2_centroids.tiff", width=8, height=8, units="in", res=600, bg="transparent", family = "sans", pointsize = 12,  compression="lzw")
METf2 %>% 
  select(centnmds1, centnmds2, Treatment, dpc) %>% 
  unique() %>% 
  arrange(dpc) %>% 
  ggplot(aes(x=centnmds1, y=centnmds2, color=Treatment)) +
  geom_point(data=METf2, aes(x=NMDS1, y=NMDS2, color=Treatment), alpha=.5) +
  geom_segment(data=METf2, aes(x=NMDS1, y=NMDS2, xend=centnmds1, yend=centnmds2, color=Treatment), alpha=.25) +
  geom_path(aes(group=Treatment)) +
  geom_point(size=6) +
  geom_text(aes(label=dpc), color='black') +
  ggtitle('NMDS of all treatment groups over time') +
  xlab('NMDS1') +
  ylab('NMDS2') +
  scale_color_manual(values=c(MVC='skyblue', EVC='orange')) +
  theme_bw()
dev.off()



# All treatments and days post Salmonella challenge (all groups over time) (centroids)
#Fecal
tiff(filename="./graphics/I4512ivac_fecal_alltreatments_centroids.tiff", width=8, height=8, units="in", res=600, bg="transparent", family = "sans", pointsize = 12,  compression="lzw")
METf %>% 
  select(centnmds1, centnmds2, Treatment, dpc) %>% 
  unique() %>% 
  arrange(dpc) %>% 
  ggplot(aes(x=centnmds1, y=centnmds2, color=Treatment, group=Treatment)) +
  geom_path(aes(group=Treatment)) +
  geom_point(size=6) +
  geom_text(aes(label=dpc), color="black") +
  ggtitle('NMDS of all treatment groups over time (centroids)') +
  xlab('NMDS1') +
  ylab('NMDS2') +
  scale_color_manual(values=c(MVMC='green', MVC='skyblue', EVC='orange')) +
  theme_bw()
dev.off()
#Excluding MVMC 
tiff(filename="./graphics/I4512ivac_fecal_alltreatments_noMVMC_centroids.tiff", width=8, height=8, units="in", res=600, bg="transparent", family = "sans", pointsize = 12,  compression="lzw")
sub1 %>% 
  select(centnmds1, centnmds2, Treatment, dpc) %>% 
  unique() %>% 
  arrange(dpc) %>% 
  ggplot(aes(x=centnmds1, y=centnmds2, color=Treatment, group=Treatment)) +
  geom_path(aes(group=Treatment)) +
  geom_point(size=6) +
  geom_text(aes(label=dpc), color="black") +
  ggtitle('NMDS of treatment groups over time (centroids)') +
  xlab('NMDS1') +
  ylab('NMDS2') +
  scale_color_manual(values=c(MVC='skyblue', EVC='orange')) +
  theme_bw()
dev.off()
#No MVMC 
tiff(filename="./graphics/I4512ivac_fecal_alltreatments_noMVMC_centroids2.tiff", width=8, height=8, units="in", res=600, bg="transparent", family = "sans", pointsize = 12,  compression="lzw")
METf2 %>% 
  select(centnmds1, centnmds2, Treatment, dpc) %>% 
  unique() %>% 
  arrange(dpc) %>% 
  ggplot(aes(x=centnmds1, y=centnmds2, color=Treatment, group=Treatment)) +
  geom_path(aes(group=Treatment)) +
  geom_point(size=6) +
  geom_text(aes(label=dpc), color="black") +
  ggtitle('NMDS of treatment groups over time (centroids)') +
  xlab('NMDS1') +
  ylab('NMDS2') +
  scale_color_manual(values=c(MVC='skyblue', EVC='orange')) +
  theme_bw()
dev.off()

# (A) VACCINATION EFFECT ON THE PORCINE FECAL MICROBIOME (MVC and EVC; D_1 or dpc -1; Fecal)
ggplot(subset(METf, Treatment %in% c("MVC", "EVC") & dpc %in% c("-1")), aes(x=NMDS1, y=NMDS2, color=Vaccine)) +
  geom_point() +
  geom_segment(aes(xend=centnmds1, yend=centnmds2)) +
  scale_color_manual(values=c(Mock='skyblue', Enterisol='orange')) +
  theme_bw() +
  facet_wrap(~Days_post_challenge, scales =
               "free") 
#OR
subA <- subset(METf, Treatment %in% c("MVC", "EVC") & dpc %in% c("-1"))
tiff(filename="./graphics/I4512ivac_fecal_vaceff_nmds_centroids.tiff", width=8, height=8, units="in", res=600, bg="transparent", family = "sans", pointsize = 12,  compression="lzw")
subA %>% 
  select(centnmds1, centnmds2, Vaccine, dpc) %>% 
  unique() %>% 
  arrange(dpc) %>% 
  ggplot(aes(x=centnmds1, y=centnmds2, color=Vaccine)) +
  geom_point(data=subA, aes(x=NMDS1, y=NMDS2, color=Vaccine), alpha=.5) +
  geom_segment(data=subA, aes(x=NMDS1, y=NMDS2, xend=centnmds1, yend=centnmds2, color=Vaccine), alpha=.25) +
  geom_path(aes(group=Vaccine)) +
  geom_point(size=6) +
  geom_text(aes(label=dpc), color='black') +
  ggtitle('NMDS of MVC (Mock_I 4,[5],12:i:-) vs EVC (Enterisol_I 4,[5],12:i:-)', 
          '1 day prior challenge and 27 days after vaccination') +
  xlab('NMDS1') +
  ylab('NMDS2') +
  scale_color_manual(values=c(Mock='skyblue', Enterisol='orange')) +
  theme_bw()
dev.off()
#Including dpc -28 and -1 day prior to vaccine
subAa <- subset(METf, Treatment %in% c("MVC", "EVC") & dpc %in% c("-1", "-28"))
tiff(filename="./graphics/I4512ivac_fecal_vaceff_with-28_nmds_centroids.tiff", width=8, height=8, units="in", res=600, bg="transparent", family = "sans", pointsize = 12,  compression="lzw")
subAa %>% 
  select(centnmds1, centnmds2, Vaccine, dpc) %>% 
  unique() %>% 
  arrange(dpc) %>% 
  ggplot(aes(x=centnmds1, y=centnmds2, color=Vaccine)) +
  geom_point(data=subAa, aes(x=NMDS1, y=NMDS2, color=Vaccine), alpha=.5) +
  geom_segment(data=subAa, aes(x=NMDS1, y=NMDS2, xend=centnmds1, yend=centnmds2, color=Vaccine), alpha=.25) +
  geom_path(aes(group=Vaccine)) +
  geom_point(size=6) +
  geom_text(aes(label=dpc), color='black') +
  ggtitle('NMDS of MVC (Mock_I 4,[5],12:i:-) vs EVC (Enterisol_I 4,[5],12:i:-)', 
          '1 day prior to vaccine (dpc -28) and 27 days after vaccination (dpc -1)') +
  xlab('NMDS1') +
  ylab('NMDS2') +
  scale_color_manual(values=c(Mock='skyblue', Enterisol='orange')) +
  theme_bw()
dev.off()
#No MVMC
subA2 <- subset(METf2, Treatment %in% c("MVC", "EVC") & dpc %in% c("-1"))
tiff(filename="./graphics/I4512ivac_fecal_vaceff_nmds_centroids2.tiff", width=8, height=8, units="in", res=600, bg="transparent", family = "sans", pointsize = 12,  compression="lzw")
subA2 %>% 
  select(centnmds1, centnmds2, Vaccine, dpc) %>% 
  unique() %>% 
  arrange(dpc) %>% 
  ggplot(aes(x=centnmds1, y=centnmds2, color=Vaccine)) +
  geom_point(data=subA2, aes(x=NMDS1, y=NMDS2, color=Vaccine), alpha=.5) +
  geom_segment(data=subA2, aes(x=NMDS1, y=NMDS2, xend=centnmds1, yend=centnmds2, color=Vaccine), alpha=.25) +
  geom_path(aes(group=Vaccine)) +
  geom_point(size=6) +
  geom_text(aes(label=dpc), color='black') +
  ggtitle('NMDS of MVC (Mock_I 4,[5],12:i:-) vs EVC (Enterisol_I 4,[5],12:i:-)', 
          '1 day prior challenge and 27 days after vaccination') +
  xlab('NMDS1') +
  ylab('NMDS2') +
  scale_color_manual(values=c(Mock='skyblue', Enterisol='orange')) +
  theme_bw()
dev.off()
#Including dpc -28 and -1 day prior to vaccine
subAa2 <- subset(METf2, Treatment %in% c("MVC", "EVC") & dpc %in% c("-1", "-28"))
tiff(filename="./graphics/I4512ivac_fecal_vaceff_with-28_nmds_centroids2.tiff", width=8, height=8, units="in", res=600, bg="transparent", family = "sans", pointsize = 12,  compression="lzw")
subAa2 %>% 
  select(centnmds1, centnmds2, Vaccine, dpc) %>% 
  unique() %>% 
  arrange(dpc) %>% 
  ggplot(aes(x=centnmds1, y=centnmds2, color=Vaccine)) +
  geom_point(data=subAa2, aes(x=NMDS1, y=NMDS2, color=Vaccine), alpha=.5) +
  geom_segment(data=subAa2, aes(x=NMDS1, y=NMDS2, xend=centnmds1, yend=centnmds2, color=Vaccine), alpha=.25) +
  geom_path(aes(group=Vaccine)) +
  geom_point(size=6) +
  geom_text(aes(label=dpc), color='black') +
  ggtitle('NMDS of MVC (Mock_I 4,[5],12:i:-) vs EVC (Enterisol_I 4,[5],12:i:-)', 
          '1 day prior to vaccine (dpc -28) and 27 days after vaccination (dpc -1)') +
  xlab('NMDS1') +
  ylab('NMDS2') +
  scale_color_manual(values=c(Mock='skyblue', Enterisol='orange')) +
  theme_bw()
dev.off()

# (B) SALMONELLA CHALLENGE EFFECT ON PORCINE FECAL MICROBIOME (MVMC and MVC; dpc 2, 3, 7, 10, 14; Fecal)
#All dpc
subB <- subset(METf, Treatment %in% c("MVMC", "MVC"))
tiff(filename="./graphics/I4512ivac_fecal_salchaeff_alldpc_nmds_centroids.tiff", width=8, height=8, units="in", res=600, bg="transparent", family = "sans", pointsize = 12,  compression="lzw")
subB %>% 
  select(centnmds1, centnmds2, Challenge, dpc) %>% 
  unique() %>% 
  arrange(dpc) %>% 
  ggplot(aes(x=centnmds1, y=centnmds2, color=Challenge)) +
  geom_point(data=subB, aes(x=NMDS1, y=NMDS2, color=Challenge), alpha=.5) +
  geom_segment(data=subB, aes(x=NMDS1, y=NMDS2, xend=centnmds1, yend=centnmds2, color=Challenge), alpha=.25) +
  geom_path(aes(group=Challenge)) +
  geom_point(size=6) +
  geom_text(aes(label=dpc), color='black') +
  ggtitle('NMDS of MVMC (Mock_Mock) vs MVC (Mock_I 4,[5],12:i:-)') +
  xlab('NMDS1') +
  ylab('NMDS2') +
  scale_color_manual(values=c(Mock='green', I4512i='skyblue')) +
  theme_bw()
dev.off()
#Excluding dpc -28 and dpc -1 (only including days after Salmonella challenge)
subBb <- subset(METf, Treatment %in% c("MVMC", "MVC") & !dpc %in% c("-28", "-1"))
tiff(filename="./graphics/I4512ivac_fecal_salchaeff_nodpc-28-1_nmds_centroids.tiff", width=8, height=8, units="in", res=600, bg="transparent", family = "sans", pointsize = 12,  compression="lzw")
subBb %>% 
  select(centnmds1, centnmds2, Challenge, dpc) %>% 
  unique() %>% 
  arrange(dpc) %>% 
  ggplot(aes(x=centnmds1, y=centnmds2, color=Challenge)) +
  geom_point(data=subBb, aes(x=NMDS1, y=NMDS2, color=Challenge), alpha=.5) +
  geom_segment(data=subBb, aes(x=NMDS1, y=NMDS2, xend=centnmds1, yend=centnmds2, color=Challenge), alpha=.25) +
  geom_path(aes(group=Challenge)) +
  geom_point(size=6) +
  geom_text(aes(label=dpc), color='black') +
  ggtitle('NMDS of MVMC (Mock_Mock) vs MVC (Mock_I 4,[5],12:i:-)', 'At different days after challenge') +
  xlab('NMDS1') +
  ylab('NMDS2') +
  scale_color_manual(values=c(Mock='green', I4512i='skyblue')) +
  theme_bw()
dev.off()

# (C) VACCINATION EFFECT ON SALMONELLA CHALLENGED PORCINE FECAL MICROBIOME (MVC and EVC; dpc 2, 3, 7, 10, 14; Fecal)
#All dpc (including -28)
subC <- subset(METf, Treatment %in% c("MVC", "EVC"))
tiff(filename="./graphics/I4512ivac_fecal_vaceffonsalcha_alldpc_nmds_centroids.tiff", width=8, height=8, units="in", res=600, bg="transparent", family = "sans", pointsize = 12,  compression="lzw")
subC %>% 
  select(centnmds1, centnmds2, Vaccine, dpc) %>% 
  unique() %>% 
  arrange(dpc) %>% 
  ggplot(aes(x=centnmds1, y=centnmds2, color=Vaccine)) +
  geom_point(data=subC, aes(x=NMDS1, y=NMDS2, color=Vaccine), alpha=.5) +
  geom_segment(data=subC, aes(x=NMDS1, y=NMDS2, xend=centnmds1, yend=centnmds2, color=Vaccine), alpha=.25) +
  geom_path(aes(group=Vaccine)) +
  geom_point(size=6) +
  geom_text(aes(label=dpc), color='black') +
  ggtitle('NMDS of MVC (Mock_I 4,[5],12:i:-) vs EVC (Enterisol_I 4,[5],12:i:-)') +
  xlab('NMDS1') +
  ylab('NMDS2') +
  scale_color_manual(values=c(Mock='skyblue', Enterisol='orange')) +
  theme_bw()
dev.off()
#Excluding dpc -28 and dpc -1 (only including days after Salmonella challenge)
subCc <- subset(METf, Treatment %in% c("MVC", "EVC") & !dpc %in% c("-28", "-1"))
tiff(filename="./graphics/I4512ivac_fecal_vaceffonsalcha__nodpc-28-1_nmds_centroids.tiff", width=8, height=8, units="in", res=600, bg="transparent", family = "sans", pointsize = 12,  compression="lzw")
subCc %>% 
  select(centnmds1, centnmds2, Vaccine, dpc) %>% 
  unique() %>% 
  arrange(dpc) %>% 
  ggplot(aes(x=centnmds1, y=centnmds2, color=Vaccine)) +
  geom_point(data=subCc, aes(x=NMDS1, y=NMDS2, color=Vaccine), alpha=.5) +
  geom_segment(data=subCc, aes(x=NMDS1, y=NMDS2, xend=centnmds1, yend=centnmds2, color=Vaccine), alpha=.25) +
  geom_path(aes(group=Vaccine)) +
  geom_point(size=6) +
  geom_text(aes(label=dpc), color='black') +
  ggtitle('NMDS of MVC (Mock_I 4,[5],12:i:-) vs EVC (Enterisol_I 4,[5],12:i:-)') +
  xlab('NMDS1') +
  ylab('NMDS2') +
  scale_color_manual(values=c(Mock='skyblue', Enterisol='orange')) +
  theme_bw()
dev.off()
#No MVMC
#All dpc (including -28)
subC2 <- subset(METf2, Treatment %in% c("MVC", "EVC"))
tiff(filename="./graphics/I4512ivac_fecal_vaceffonsalcha_alldpc_nmds_centroids2.tiff", width=8, height=8, units="in", res=600, bg="transparent", family = "sans", pointsize = 12,  compression="lzw")
subC2 %>% 
  select(centnmds1, centnmds2, Vaccine, dpc) %>% 
  unique() %>% 
  arrange(dpc) %>% 
  ggplot(aes(x=centnmds1, y=centnmds2, color=Vaccine)) +
  geom_point(data=subC2, aes(x=NMDS1, y=NMDS2, color=Vaccine), alpha=.5) +
  geom_segment(data=subC2, aes(x=NMDS1, y=NMDS2, xend=centnmds1, yend=centnmds2, color=Vaccine), alpha=.25) +
  geom_path(aes(group=Vaccine)) +
  geom_point(size=6) +
  geom_text(aes(label=dpc), color='black') +
  ggtitle('NMDS of MVC (Mock_I 4,[5],12:i:-) vs EVC (Enterisol_I 4,[5],12:i:-)') +
  xlab('NMDS1') +
  ylab('NMDS2') +
  scale_color_manual(values=c(Mock='skyblue', Enterisol='orange')) +
  theme_bw()
dev.off()
#Excluding dpc -28 and dpc -1 (only including days after Salmonella challenge)
subCc2 <- subset(METf2, Treatment %in% c("MVC", "EVC") & !dpc %in% c("-28", "-1"))
tiff(filename="./graphics/I4512ivac_fecal_vaceffonsalcha__nodpc-28-1_nmds_centroids2.tiff", width=8, height=8, units="in", res=600, bg="transparent", family = "sans", pointsize = 12,  compression="lzw")
subCc2 %>% 
  select(centnmds1, centnmds2, Vaccine, dpc) %>% 
  unique() %>% 
  arrange(dpc) %>% 
  ggplot(aes(x=centnmds1, y=centnmds2, color=Vaccine)) +
  geom_point(data=subCc2, aes(x=NMDS1, y=NMDS2, color=Vaccine), alpha=.5) +
  geom_segment(data=subCc2, aes(x=NMDS1, y=NMDS2, xend=centnmds1, yend=centnmds2, color=Vaccine), alpha=.25) +
  geom_path(aes(group=Vaccine)) +
  geom_point(size=6) +
  geom_text(aes(label=dpc), color='black') +
  ggtitle('NMDS of MVC (Mock_I 4,[5],12:i:-) vs EVC (Enterisol_I 4,[5],12:i:-)') +
  xlab('NMDS1') +
  ylab('NMDS2') +
  scale_color_manual(values=c(Mock='skyblue', Enterisol='orange')) +
  theme_bw()
dev.off()

# (D) SALMONELLA CHALLENGE EFFECT ON PORCINE CECAL MICROBIOME (MVMC and MVC; D14 or dpc 14; Cecal_contents)
subD <- subset(METcc, Treatment %in% c("MVMC", "MVC"))
tiff(filename="./graphics/I4512ivac_cecalcontents_salchaeff_dpc14_nmds_centroids.tiff", width=8, height=8, units="in", res=600, bg="transparent", family = "sans", pointsize = 12,  compression="lzw")
subD %>% 
  select(centnmds1, centnmds2, Challenge, dpc) %>% 
  unique() %>% 
  arrange(dpc) %>% 
  ggplot(aes(x=centnmds1, y=centnmds2, color=Challenge)) +
  geom_point(data=subD, aes(x=NMDS1, y=NMDS2, color=Challenge), alpha=.5) +
  geom_segment(data=subD, aes(x=NMDS1, y=NMDS2, xend=centnmds1, yend=centnmds2, color=Challenge), alpha=.25) +
  geom_path(aes(group=Challenge)) +
  geom_point(size=6) +
  geom_text(aes(label=dpc), color='black') +
  ggtitle('NMDS of MVMC (Mock_Mock) vs MVC (Mock_I 4,[5],12:i:-)',
          'Cecal contents at days post challenge 14') +
  xlab('NMDS1') +
  ylab('NMDS2') +
  scale_color_manual(values=c(Mock='green', I4512i='skyblue')) +
  theme_bw()
dev.off()

# (E) VACCINATION EFFECT ON SALMONELLA CHALLENGED PORCINE CECAL MICROBIOME (MVC and EVC; D14 or dpc 14; Cecal_contents)
subE <- subset(METcc, Treatment %in% c("MVC", "EVC"))
tiff(filename="./graphics/I4512ivac_cecalcontents_vacchaeff_dpc14_nmds_centroids.tiff", width=8, height=8, units="in", res=600, bg="transparent", family = "sans", pointsize = 12,  compression="lzw")
subE %>% 
  select(centnmds1, centnmds2, Vaccine, dpc) %>% 
  unique() %>% 
  arrange(dpc) %>% 
  ggplot(aes(x=centnmds1, y=centnmds2, color=Vaccine)) +
  geom_point(data=subE, aes(x=NMDS1, y=NMDS2, color=Vaccine), alpha=.5) +
  geom_segment(data=subE, aes(x=NMDS1, y=NMDS2, xend=centnmds1, yend=centnmds2, color=Vaccine), alpha=.25) +
  geom_path(aes(group=Vaccine)) +
  geom_point(size=6) +
  geom_text(aes(label=dpc), color='black') +
  ggtitle('NMDS of MVC (Mock_I 4,[5],12:i:-) vs EVC (Enterisol_I 4,[5],12:i:-)',
          'Cecal contents at days post challenge 14') +
  xlab('NMDS1') +
  ylab('NMDS2') +
  scale_color_manual(values=c(Mock='skyblue', Enterisol='orange')) +
  theme_bw()
dev.off()

# (F) ALL PORCINE CECAL MICROBIOME (MVC, EVC, MVMC; D14 or dpc 14; Cecal_contents)
tiff(filename="./graphics/I4512ivac_cecalcontents_alltreatments_nmds_centroids.tiff", width=8, height=8, units="in", res=600, bg="transparent", family = "sans", pointsize = 12,  compression="lzw")
METcc %>% 
  dplyr::select(centnmds1, centnmds2, Treatment, dpc) %>% 
  unique() %>% 
  arrange(dpc) %>% 
  ggplot(aes(x=centnmds1, y=centnmds2, color=Treatment)) +
  geom_point(data=METcc, aes(x=NMDS1, y=NMDS2, color=Treatment), alpha=.5) +
  geom_segment(data=METcc, aes(x=NMDS1, y=NMDS2, xend=centnmds1, yend=centnmds2, color=Treatment), alpha=.25) +
  geom_path(aes(group=Treatment)) +
  geom_point(aes(shape = Treatment), size=6) +
  geom_text(aes(label=dpc), color='black') +
  ggtitle('NMDS of all treatment groups 14 days post challenge', 
          'Cecal contents') +
  xlab('NMDS1') +
  ylab('NMDS2') +
  scale_color_manual(values=c(MVMC='green', MVC='skyblue', EVC='orange')) +
  annotate("text", x=0.32, y=-0.20, label="Stress = 0.133", size = 4) +
  theme_bw()
dev.off()


### Calculate distances using vegdist (not recommended)
#set.seed(19760620) #NMDS uses random # generator so it's consistent run after run
#dist <- vegdist(OTUmat, method="bray") #creates/calculates lower triangle distance matrix
#nmds <- metaMDS(DIST) #run the ordination - random algorithm
#str(nmds) #list, see points matrix (position of each diff samples)
#scores(nmds) #position of all points on axis, stored as matrix (to plot with ggplot make into df)
#plot(nmds) #ordination in base R graphics
#### Plot ordination data with ggplot2, make NMDS into df
#scores(NMDS) %>% 
#  as_tibble(rownames="Samples") %>% 
#  ggplot(aes(x=NMDS1, y=NMDS2)) +
#  geom_point()

# Other useful commands
NMDSf$points #sample scores, matrix that contains all points
NMDSf$dims #number of NMS axes or dimensions
NMDSf$stress #stress value of final solution
NMDSf$data #what was ordinated, including any transformations
NMDSf$distance #distance metric used
NMDSf$converged #whether solution converged or not (T/F)
NMDSf$tries #number of random initial configurations tried
NMDSf$species #scores of variables (species / taxa in ecology)
NMDSf$call #restates how the function was called
