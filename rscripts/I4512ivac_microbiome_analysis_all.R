# I 4,[5],12:i:- MICROBIOME ANALYSIS
# Selma Burciaga
# January 10, 2023
# Rscript: I4512ivac_veganadonis_commands.R
# R v4.2.2/R v2022.12.0

# PART 1: Using the vegan R package to calculate ecological distances (CC188) - https://riffomonas.org/code_club/2022-02-17-distances
# PART 2: Using the pairwiseAdonis GitHub package to run PERMANOVA pairwise comparisons tests - https://github.com/pmartinezarbizu/pairwiseAdonis



# PART 1 BETA DIVERSITY: NMDS___________________________________________________
## To use package: library(package), packageVersion("package")
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
sample_data(PHYLO)$Vaccine <- factor(sample_data(PHYLO)$Vaccince, levels=c('Mock', 'Enterisol'))
sample_data(PHYLO)$Challenge <- factor(sample_data(PHYLO)$Challenge, levels=c('Mock', 'I4512i'))

## Create subsets (Fecal and Cecal_contents)
PHYLOf <- prune_samples(samples = PHYLO@sam_data$Sample_Type=="Fecal", x = PHYLO)
PHYLOcc <- prune_samples(samples = PHYLO@sam_data$Sample_Type=="Cecal_contents", x = PHYLO)
#No MVMC
PHYLOf2 <- prune_samples(samples = PHYLOf@sam_data$Challenge=="I4512i", x = PHYLOf)

## Create OTU R matrix and metadata R objects 
#Fecal
OTUmatf <- as(otu_table(PHYLOf), "matrix") #creating R matrix from phyloseq object 
METf <- sample_data(PHYLOf) %>% 
  as(Class="data.frame") #creating df from phyloseq object
#Cecal_contents
OTUmatcc <- as(otu_table(PHYLOcc), "matrix")
METcc <- sample_data(PHYLOcc) %>% 
  as(Class="data.frame")
#No MVMC
OTUmatf2 <- as(otu_table(PHYLOf2), "matrix")
METf2 <- sample_data(PHYLOf2) %>% 
  as(Class="data.frame")


## Rarefying distance calculations using avgdist (recommended)
set.seed(19760620) #usually before func call, best for reproducibility

### If file doesn't exist, calculate dist and save, if it does read file and save as DIST (sample= min otu counts)
#Fecal
if (!file.exists("outputs/I4512ivac_fecal_dist.rds")) {
  DISTf <- avgdist(OTUmatf, dmethod="bray", sample=min(rowSums(OTUmatf))) 
  write_rds(DISTf, "outputs/I4512ivac_fecal_dist.rds")
} else {
  DISTf <- readRDS("outputs/I4512ivac_fecal_dist.rds")
}
#Cecal_contents
if (!file.exists("outputs/I4512ivac_cecalcontents_dist.rds")) {
  DISTcc <- avgdist(OTUmatcc, dmethod="bray", sample=min(rowSums(OTUmatcc)))
  write_rds(DISTcc, "outputs/I4512ivac_cecalcontents_dist.rds")
} else {
  DISTcc <- readRDS("outputs/I4512ivac_cecalcontents_dist.rds")
}
#No MVMC
if (!file.exists("outputs/I4512ivac_fecal_dist2.rds")) {
  DISTf2 <- avgdist(OTUmatf2, dmethod="bray", sample=min(rowSums(OTUmatf2))) 
  write_rds(DISTf2, "outputs/I4512ivac_fecal_dist2.rds")
} else {
  DISTf2 <- readRDS("outputs/I4512ivac_fecal_dist2.rds")
}

## Non-metric Multidimensional Scaling - create NMDS points and adhere to MET objects
#Fecal
set.seed(1)
NMDSf <- metaMDS(DISTf)
NMDSscoresf <- scores(NMDSf) %>% 
  as_tibble(rownames="SampleID")
METf <- METf %>% 
  left_join(NMDSscoresf)
#Cecal_contents
set.seed(1)
NMDScc <- metaMDS(DISTcc)
NMDSscorescc <- scores(NMDScc) %>% 
  as_tibble(rownames="SampleID")
METcc <- METcc %>% 
  left_join(NMDSscorescc)
#No MVMC
set.seed(1)
NMDSf2 <- metaMDS(DISTf2)
NMDSscoresf2 <- scores(NMDSf2) %>% 
  as_tibble(rownames="SampleID")
METf2 <- METf2 %>% 
  left_join(NMDSscoresf2)

#Calculation of ordination stress
NMDSf$stress
NMDScc$stress

# View order of variables in METf columns
METf$Days_post_challenge %>% 
  unique()
METf$Treatment %>% 
  unique()

# Group_by dpc and Treatment - describing centroids
#Fecal
METf <- METf %>% 
  group_by(dpc, Treatment) %>% 
  mutate(centnmds1=mean(NMDS1),
         centnmds2=mean(NMDS2)) %>% 
  ungroup()
#Cecal_contents
METcc <- METcc %>% 
  group_by(dpc, Treatment) %>% 
  mutate(centnmds1=mean(NMDS1),
         centnmds2=mean(NMDS2)) %>% 
  ungroup()
#No MVMC
METf2 <- METf2 %>% 
  group_by(dpc, Treatment) %>% 
  mutate(centnmds1=mean(NMDS1),
         centnmds2=mean(NMDS2)) %>% 
  ungroup()


## Bringing metadata into ordination after rarefaction 
#e.g. All treatment groups individual graphs for each day post Salmonella challenge
#Fecal
#By days post challenge
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



    #Excluding MVMC (which is usually not included in comparisons of vaccination/challenged groups)
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
#NMDS with centroids
tiff(filename="./graphics/I4512ivac_fecal_alltreatments_nmds.tiff", width=8, height=8, units="in", res=600, bg="transparent", family = "sans", pointsize = 12,  compression="lzw")
METf %>% 
  select(centnmds1, centnmds2, Treatment, dpc) %>% 
  unique() %>% 
  arrange(dpc) %>% 
  ggplot(aes(x=centnmds1, y=centnmds2, color=Treatment)) +
  geom_point(data=METf, aes(x=NMDS1, y=NMDS2, color=Treatment), alpha=.5) +
  geom_segment(data=METf, aes(x=NMDS1, y=NMDS2, xend=centnmds1, yend=centnmds2, color=Treatment), alpha=.25) +
  geom_path(aes(group=Treatment)) +
  geom_point(size=6) +
  geom_text(aes(label=dpc), color='black') +
  ggtitle('NMDS of all treatment groups over time') +
  xlab('NMDS1') +
  ylab('NMDS2') +
  scale_color_manual(values=c(MVMC='green', MVC='skyblue', EVC='orange')) +
  theme_bw()
dev.off()
    #Excluding MVMC
tiff(filename="./graphics/I4512ivac_fecal_alltreatments_nomvmc_nmds.tiff", width=8, height=8, units="in", res=600, bg="transparent", family = "sans", pointsize = 12,  compression="lzw")
sub1 %>% 
  select(centnmds1, centnmds2, Treatment, dpc) %>% 
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
tiff(filename="./graphics/I4512ivac_cc_alltreatmentsbydpc_nmds.tiff", width=8, height=8, units="in", res=600, bg="transparent", family = "sans", pointsize = 12,  compression="lzw")
METcc %>% 
  ggplot(aes(x=NMDS1, y=NMDS2, color=Treatment)) +
  geom_point() +
  geom_segment(aes(xend=centnmds1, yend=centnmds2)) +
  ggtitle('NMDS of cecal contents 14 days post challenge') +
  scale_color_manual(values=c(MVMC='green', MVC='skyblue', EVC='orange')) +
  facet_wrap(~Days_post_challenge, scales = "free") +
  theme_bw()
dev.off()
#No MVMC Fecal
    #By days post challenge
tiff(filename="./graphics/I4512ivac_fecal_alltreatmentsbydpc_nomvmc_nmds2.tiff", width=8, height=8, units="in", res=600, bg="transparent", family = "sans", pointsize = 12,  compression="lzw")
METf2 %>% 
  ggplot(aes(x=NMDS1, y=NMDS2, color=Treatment)) +
  geom_point() +
  geom_segment(aes(xend=centnmds1, yend=centnmds2)) +
  ggtitle('NMDS of all treatment groups over time by days post challenge') +
  scale_color_manual(values=c(MVC='skyblue', EVC='orange')) +
  facet_wrap(~Days_post_challenge, scales = "free") +
  theme_bw()
dev.off()
  #NMDS with centroids
tiff(filename="./graphics/I4512ivac_fecal_alltreatments_nomvmc_nmds2.tiff", width=8, height=8, units="in", res=600, bg="transparent", family = "sans", pointsize = 12,  compression="lzw")
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
tiff(filename="./graphics/I4512ivac_fecal_vaceff_nmdscentroids.tiff", width=8, height=8, units="in", res=600, bg="transparent", family = "sans", pointsize = 12,  compression="lzw")
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
tiff(filename="./graphics/I4512ivac_fecal_vaceff_with-28_nmdscentroids.tiff", width=8, height=8, units="in", res=600, bg="transparent", family = "sans", pointsize = 12,  compression="lzw")
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
tiff(filename="./graphics/I4512ivac_fecal_vaceff_nmdscentroids2.tiff", width=8, height=8, units="in", res=600, bg="transparent", family = "sans", pointsize = 12,  compression="lzw")
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
tiff(filename="./graphics/I4512ivac_fecal_vaceff_with-28_nmdscentroids2.tiff", width=8, height=8, units="in", res=600, bg="transparent", family = "sans", pointsize = 12,  compression="lzw")
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
tiff(filename="./graphics/I4512ivac_fecal_salchaeff_alldpc_nmdscentroids.tiff", width=8, height=8, units="in", res=600, bg="transparent", family = "sans", pointsize = 12,  compression="lzw")
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
tiff(filename="./graphics/I4512ivac_fecal_salchaeff_nodpc-28-1_nmdscentroids.tiff", width=8, height=8, units="in", res=600, bg="transparent", family = "sans", pointsize = 12,  compression="lzw")
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
tiff(filename="./graphics/I4512ivac_fecal_vaceffonsalcha_alldpc_nmdscentroids.tiff", width=8, height=8, units="in", res=600, bg="transparent", family = "sans", pointsize = 12,  compression="lzw")
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
tiff(filename="./graphics/I4512ivac_fecal_vaceffonsalcha__nodpc-28-1_nmdscentroids.tiff", width=8, height=8, units="in", res=600, bg="transparent", family = "sans", pointsize = 12,  compression="lzw")
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
tiff(filename="./graphics/I4512ivac_fecal_vaceffonsalcha_alldpc_nmdscentroids2.tiff", width=8, height=8, units="in", res=600, bg="transparent", family = "sans", pointsize = 12,  compression="lzw")
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
tiff(filename="./graphics/I4512ivac_fecal_vaceffonsalcha__nodpc-28-1_nmdscentroids2.tiff", width=8, height=8, units="in", res=600, bg="transparent", family = "sans", pointsize = 12,  compression="lzw")
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
tiff(filename="./graphics/I4512ivac_cc_salchaeff_dpc14_nmdscentroids.tiff", width=8, height=8, units="in", res=600, bg="transparent", family = "sans", pointsize = 12,  compression="lzw")
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
tiff(filename="./graphics/I4512ivac_cc_vacchaeff_dpc14_nmdscentroids.tiff", width=8, height=8, units="in", res=600, bg="transparent", family = "sans", pointsize = 12,  compression="lzw")
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



# PART 2 PAIRWISEADONIS_________________________________________________________
## To use package: library(package), packageVersion("package")
library(pairwiseAdonis) #0.4.1
library(RColorBrewer)

## Create new column combining treatment and dpc

#Fecal
METf <- METf %>%  
  mutate(set=paste(Treatment, dpc, sep="_"))
#Cecal_contents
METcc <- METcc %>% 
  mutate(set=paste(Treatment, dpc, sep="_"))
#No MVMC
METf2 <- METf2 %>% 
  mutate(set=paste(Treatment, dpc, sep="_"))

## Confirm order of SampleID matches in DIST and MET objects
#Fecal
all(attributes(DISTf)$Labels == METf$SampleID)
#Cecal_contents
all(attributes(DISTcc)$Labels == METcc$SampleID)
#No MVMC
all(attributes(DISTf2)$Labels == METf2$SampleID)

## Running pairwise adonis
#Fecal
if (!file.exists('./outputs/I4512ivac_fecal_all_pwadon.rds')){
  PWADONf <- pairwise.adonis(x = DISTf, factors = METf$set, perm =9999)
  write_rds(PWADONf, './outputs/I4512ivac_fecal_all_pwadon.rds')
} else {
  PWADONf <- readRDS('./outputs/I4512ivac_fecal_all_pwadon.rds')
}
#Cecal_contents
if (!file.exists('./outputs/I4512ivac_cecalcontents_all_pwadon.rds')){
  PWADONcc <- pairwise.adonis(x = DISTcc, factors = METcc$set, perm =9999)
  write_rds(PWADONcc, './outputs/I4512ivac_cecalcontents_all_pwadon.rds')
} else {
  PWADONcc <- readRDS('./outputs/I4512ivac_cecalcontents_all_pwadon.rds')
}
#No MVMC
if (!file.exists('./outputs/I4512ivac_fecal_all_pwadon2.rds')){
  PWADONf2 <- pairwise.adonis(x = DISTf2, factors = METf2$set, perm =9999)
  write_rds(PWADONf2, './outputs/I4512ivac_fecal_all_pwadon2.rds')
} else {
  PWADONf2 <- readRDS('./outputs/I4512ivac_fecal_all_pwadon2.rds')
}


## Important pairwise comparisons
colnames(PWADONf)[1] = "Pairs"
colnames(PWADONcc)[1] = "Pairs"
colnames(PWADONf2)[1] = "Pairs"
### Treatments within day (Fecal) - Within day comparisons
twdf <- PWADONf[grep("[A-Z]+_(-?[0-9]+) vs [A-Z]+_(\\1)", PWADONf$Pairs),]
twdfl <- grepl("[A-Z]+_(-?[0-9]+) vs [A-Z]+_(\\1)", PWADONf$Pairs)
twdf$Type <- "Treatment_within_day"
### Treatments over time (Fecal) - Within group comparisons
totf <- PWADONf[grep("([A-Z]+)_-?[0-9]+ vs (\\1)_-?[0-9]+", PWADONf$Pairs),]
totfl <- grepl("([A-Z]+)_-?[0-9]+ vs (\\1)_-?[0-9]+", PWADONf$Pairs)
totf$Type <- "Treatment_over_time"
### Cecal contents, dpc 14
ccc <- PWADONcc[grep("[A-Z]+_([14]+) vs [A-Z]+_(\\1)", PWADONcc$Pairs),]
cccl <- grepl("[A-Z]+_([14]+) vs [A-Z]+_(\\1)", PWADONcc$Pairs)
ccc$Type <- "Cecal_contents_comparions"
###No MVMC
    ### Treatments within day (Fecal2)
twdf2 <- PWADONf2[grep("[A-Z]+_(-?[0-9]+) vs [A-Z]+_(\\1)", PWADONf2$Pairs),]
twdfl2 <- grepl("[A-Z]+_(-?[0-9]+) vs [A-Z]+_(\\1)", PWADONf2$Pairs)
twdf2$Type <- "Treatment_within_day"
    ### Treatments over time (Fecal2)
totf2 <- PWADONf2[grep("([A-Z]+)_-?[0-9]+ vs (\\1)_-?[0-9]+", PWADONf2$Pairs),]
totfl2 <- grepl("([A-Z]+)_-?[0-9]+ vs (\\1)_-?[0-9]+", PWADONf2$Pairs)
totf2$Type <- "Treatment_over_time"

## P-value correction (adjusted p-value) using fdr method (less strict; all = Fecal and Cecal contents)
PCOMPS <- bind_rows(twdf, totf, ccc) #all important pairwise comparisons 
PCOMPS$FDR_pval <- p.adjust(PCOMPS$p.value, method = "fdr")
PCOMPS <- select(PCOMPS, -sig, -p.adjusted)
PCOMPS$FDR_sig[PCOMPS$FDR_pval < 0.045] <- "*" #significant pairs after correction
PCOMPS$FDR_sig[PCOMPS$FDR_pval >= 0.045] <- "" 
#No MVMC (except ccc)
PCOMPS2 <- bind_rows(twdf2, totf2, ccc) #all important pairwise comparisons 
PCOMPS2$FDR_pval <- p.adjust(PCOMPS2$p.value, method = "fdr")
PCOMPS2 <- select(PCOMPS2, -sig, -p.adjusted)
PCOMPS2$FDR_sig[PCOMPS2$FDR_pval < 0.045] <- "*" #significant pairs after correction
PCOMPS2$FDR_sig[PCOMPS2$FDR_pval >= 0.045] <- "" 

## Save adjusted pairwise comparisons of twdf, totf, and ccc
write_tsv(PCOMPS, "./outputs/I4512ivac_pairwise_permanova_tests_all.tsv")
#No MVMC
write_tsv(PCOMPS2, "./outputs/I4512ivac_pairwise_permanova_tests_all2.tsv")

## Set reference group - Mock vaccinated, Mock challenged (MVMC, also Mock_Mock)
pwdayMVMC <-
  PWADONf[twdfl,] %>%
  mutate(G1=sub('(.*) vs (.*)','\\1',Pairs),
         G2=sub('(.*) vs (.*)','\\2',Pairs)) %>%
  filter(grepl('MVMC',G1) |grepl('MVMC',G2)) %>%
  mutate(REF=ifelse(grepl('MVMC',G1), G1,G2),
         COMP=ifelse(grepl('MVMC',G1), G2,G1),
         Day=sub('.*_(-?[0-9]+)','\\1',REF) %>% as.numeric(),
         Treatment=sub('(.*)_-?[0-9]+','\\1',COMP),
         FDR=p.adjust(p.value, method = 'fdr'))

## Set reference group - Mock vaccinated, I5412i- challenged (MVC, also Mock_I4512i)
#Fecal
pwdayMVC <-
  PWADONf[twdfl,] %>%
  mutate(G1=sub('(.*) vs (.*)','\\1',Pairs),
         G2=sub('(.*) vs (.*)','\\2',Pairs)) %>%
  filter(grepl('MVC',G1) |grepl('MVC',G2)) %>%
  mutate(REF=ifelse(grepl('MVC',G1), G1,G2),
         COMP=ifelse(grepl('MVC',G1), G2,G1),
         Day=sub('.*_(-?[0-9]+)','\\1',REF) %>% as.numeric(),
         Treatment=sub('(.*)_-?[0-9]+','\\1',COMP),
         FDR=p.adjust(p.value, method = 'fdr'))
#No MVMC
pwdayMVC2 <-
  PWADONf2[twdfl2,] %>%
  mutate(G1=sub('(.*) vs (.*)','\\1',Pairs),
         G2=sub('(.*) vs (.*)','\\2',Pairs)) %>%
  filter(grepl('MVC',G1) |grepl('MVC',G2)) %>%
  mutate(REF=ifelse(grepl('MVC',G1), G1,G2),
         COMP=ifelse(grepl('MVC',G1), G2,G1),
         Day=sub('.*_(-?[0-9]+)','\\1',REF) %>% as.numeric(),
         Treatment=sub('(.*)_-?[0-9]+','\\1',COMP),
         FDR=p.adjust(p.value, method = 'fdr'))

## All treatment groups
#Fecal
pwdayALL <- 
  PWADONf[twdfl,] %>%
  mutate(Pairs=str_replace(Pairs, "EVC_14 vs MVC_14", "MVC_14 vs EVC_14")) %>% 
  mutate(G1=sub('(.*) vs (.*)','\\1',Pairs),
         G2=sub('(.*) vs (.*)','\\2',Pairs)) %>%
  mutate(REF=ifelse(grepl('MVMC',G1), G1,G2),
         COMP=ifelse(grepl('MVMC',G1), G2,G1),
         Day=sub('.*_(-?[0-9]+)','\\1',REF) %>% as.numeric(),
         TreatmentA=sub('(.*)_-?[0-9]+','\\1',REF),
         TreatmentB=sub('(.*)_-?[0-9]+','\\1',COMP),
         FDR=p.adjust(p.value, method = 'fdr'),
         Group=paste(TreatmentA, TreatmentB, sep=' vs '),
         Group=factor(Group, levels = c("MVMC vs MVC",
                                        "MVMC vs EVC",
                                        "EVC vs MVC")))

pwdayALL %>% 
  ggplot(aes(x=F.Model, y=fct_rev(Group))) +
  geom_col() +
  facet_wrap(~Day, nrow = 1)

### PERMANOVA: Treatment group comparisons at each day
#### Treatment groups at each day (Fecal; bar graph)
tiff(filename="./graphics/I4512ivac_permanova_groupsateachday.tiff", width=8, height=8, units="in", res=600, bg="transparent", family = "sans", pointsize = 12,  compression="lzw")
pwdayALL %>%
  ggplot(aes(x=F.Model, y=fct_rev(Group))) +
  geom_col() +
  facet_wrap(~Day, ncol = 1) +
  ggtitle('PERMANOVA tests between treatment groups at each day')
dev.off()
    #With p-values and color by Group pairs
tiff(filename="./graphics/I4512ivac_permanova_groupsateachday_pvalues.tiff", width=8, height=8, units="in", res=600, bg="transparent", family = "sans", pointsize = 12,  compression="lzw")
pwdayALL %>%
  ggplot(aes(x=F.Model, y=fct_rev(Group), fill=Group)) +
  geom_col() +
  geom_text(aes(label=paste0('P=',signif(FDR,2))),hjust=-0.1, 
            position=position_dodge(width=1) ) +
  facet_wrap(~Day, ncol = 1) +
  scale_fill_brewer(palette="Blues", )
  ggtitle('PERMANOVA tests between treatment groups at each day')
  theme_bw()
dev.off()
#### Difference from MVMC (Mock_Mock) at each day
tiff(filename="./graphics/I4512ivac_permanova_vs_MVMC.tiff", width=8, height=8, units="in", res=600, bg="transparent", family = "sans", pointsize = 12,  compression="lzw")
pwdayMVMC %>%
  ggplot(aes(x=Day, y=F.Model, color=Treatment))+
  geom_point() +
  geom_line() +
  ylim(0,25) +
  ggtitle('Difference from MVMC (Mock_I 4,[5],12:i:-) at each day', 'All comparisons within each day are significant')+
  scale_color_manual(values=c(MVC='skyblue', EVC='orange')) +
  theme_bw()
dev.off()
#### Difference from MVC (Mock_I 4,[5],12:i:-) at each day
tiff(filename="./graphics/I4512ivac_permanova_vs_MVC.tiff", width=8, height=8, units="in", res=600, bg="transparent", family = "sans", pointsize = 12,  compression="lzw")
pwdayMVC %>%
  ggplot(aes(x=Day, y=F.Model, color=Treatment))+
  geom_point() +
  geom_line() +
  ylim(0,25) +
  ggtitle('Difference from MVC (Mock_I 4,[5],12:i:-) at each day', 'All comparisons within each day are significant (except EVC day -28)')+
  scale_color_manual(values=c(MVMC='green', EVC='orange')) +
  theme_bw()
dev.off()
    #No MVMC
tiff(filename="./graphics/I4512ivac_permanova_vs_MVC2.tiff", width=8, height=8, units="in", res=600, bg="transparent", family = "sans", pointsize = 12,  compression="lzw")
pwdayMVC2 %>%
  ggplot(aes(x=Day, y=F.Model, color=Treatment))+
  geom_point() +
  geom_line() +
  ylim(0,25) +
  ggtitle('Difference from MVC (Mock_I 4,[5],12:i:-) at each day', 'All comparisons within each day are significant (except EVC day -28)')+
  scale_color_manual(values=c(EVC='orange')) +
  theme_bw()
dev.off()
#### Within group changes over time (comparing back to themselves at day -28)
pw_group <- PWADONf[totfl,] %>% 
  mutate(Pairs=str_replace(Pairs, 
        "MVC_-1 vs MVC_-28", "MVC_-28 vs MVC_-1")) %>%
  mutate(Pairs=str_replace(Pairs, 
        "MVC_-1 vs MVC_14", "MVC_14 vs MVC_-1")) %>%
  mutate(Pairs=str_replace(Pairs, 
        "EVC_-1 vs EVC_-28", "EVC_-28 vs EVC_-1")) %>%
  mutate(Pairs=str_replace(Pairs, 
        "MVMC_-1 vs MVMC_-28", "MVMC_-28 vs MVMC_-1")) %>%
  mutate(G1=sub('(.*) vs (.*)','\\1',Pairs),
         G2=sub('(.*) vs (.*)','\\2',Pairs)) %>%
  filter(grepl('-1',G1) | grepl('-1',G2)) %>%
  mutate(REF=ifelse(grepl('-1',G2), G1,G2),
         COMP=ifelse(grepl('-1',G2), G2,G1),
         Day=sub('\\1','*_(-?[0-9]+).',COMP) %>% as.numeric(),
         Treatment=sub('(.*)_(-?[0-9]+)','\\2',REF),
         fdr_pval=p.adjust(p.value, method = 'fdr')) 
 









# ALPHA DIVERSITY (vegan)
set.seed(19760620)

TAXf <-
  PHYLOf@tax_table %>%
  as(., Class = 'matrix') %>%
  as.data.frame() %>%
  rownames_to_column(var = 'ASV')

TAXcc <-
  PHYLOcc@tax_table %>%
  as(., Class = 'matrix') %>%
  as.data.frame() %>%
  rownames_to_column(var = 'ASV')

EFGf <-
  PHYLOf %>%
  subset_samples(PHYLOf@sam_data$Vaccine %in% c('Mock', 'Enterisol') &
                   PHYLOf@sam_data$Challenge %in% c('Mock', 'I4512i')) %>%
  prune_taxa(taxa = taxa_sums(.) > 0) %>%
  prune_taxa(taxa=colSums(.@otu_table > 0) > 5) # taxa detected in more than 5 samples

EFGf@sam_data$dpc <- EFGf@sam_data$Days_post_challenge
EFGf@sam_data$set <- paste(EFGf@sam_data$Days_post_challenge, EFGf@sam_data$Vaccine,EFGf@sam_data$Challenge, sep = '_')
EFGf@sam_data$Group <- paste(EFGf@sam_data$Vaccine,EFGf@sam_data$Challenge, sep = '_')

EFGcc <-
  PHYLOcc %>%
  subset_samples(PHYLOcc@sam_data$Vaccine %in% c('Mock', 'Enterisol') &
                   PHYLOcc@sam_data$Challenge %in% c('Mock', 'I4512i')) %>%
  prune_taxa(taxa = taxa_sums(.) > 0) %>%
  prune_taxa(taxa=colSums(.@otu_table > 0) > 5) # taxa detected in more than 5 samples

EFGcc@sam_data$dpc <- EFGcc@sam_data$Days_post_challenge
EFGcc@sam_data$set <- paste(EFGcc@sam_data$Days_post_challenge, EFGcc@sam_data$Vaccine,EFGcc@sam_data$Challenge, sep = '_')
EFGcc@sam_data$Group <- paste(EFGcc@sam_data$Vaccine,EFGcc@sam_data$Challenge, sep = '_')


# Plot_richness (phyloseq) - total # OTUs/taxa in sample
plot_richness(PHYLOf, x="Treatment", color="Treatment", measures=c("Shannon", "Simpson", "Chao1"))
plot_richness(PHYLOcc, x="Treatment", color="Treatment", measures=c("Shannon", "Simpson", "Chao1"))

## see how spread out your data actually is overall (all points not constrained by 25%/75% percentile)
#Fecal - EDIT TO UPDATE COLORS - not geomboxplot
tiff(filename="./graphics/I4512ivac_fecal_plotrichness_groups.tiff", width=8, height=8, units="in", res=600, bg="transparent", family = "sans", pointsize = 12,  compression="lzw")
richf <- plot_richness(PHYLOf, x="Treatment", color="Treatment", measures=c("Shannon")) %>%
  ggplot(aes(x=Group, y=Shannon, group=set, fill=Group)) +
  geom_point() +
  scale_fill_manual(values=c(MVMC='green', MVC='skyblue', EVC='orange')) +
  theme_bw()
richf + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
dev.off()

#Cecal_contents


##Rarefy_even_depth
EFGf_rare <- rarefy_even_depth(EFGf)
EFGcc_rare <- rarefy_even_depth(EFGcc)

### Shannon's Diversity
#Fecal
EFGf_rare@sam_data$Shannon <- diversity(EFGf_rare@otu_table)

tiff(filename="./graphics/I4512ivac_fecal_shannon_groups.tiff", width=8, height=8, units="in", res=600, bg="transparent", family = "sans", pointsize = 12,  compression="lzw")
shanf <- EFGf_rare@sam_data %>%
  ggplot(aes(x=Group, y=Shannon, group=set, fill=Group)) +
  geom_boxplot() +
  facet_wrap(~Days_post_challenge, nrow=1) +
  scale_fill_manual(values=c(Mock_Mock='green', Mock_I4512i='skyblue', Enterisol_I4512i='orange')) +
  theme_bw()
shanf + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
dev.off()

tiff(filename="./graphics/I4512ivac_fecal_shannon_groups2.tiff", width=8, height=8, units="in", res=600, bg="transparent", family = "sans", pointsize = 12,  compression="lzw")
shan1f <- EFGf_rare@sam_data %>%
  ggplot(aes(x=Treatment, y=Shannon, group=set, fill=Treatment)) +
  geom_boxplot() +
  facet_wrap(~Days_post_challenge, nrow=1) +
  scale_fill_manual(values=c(MVMC='green', MVC='skyblue', EVC='orange')) +
  theme_bw()
shan1f + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
dev.off()

#Cecal_contents
EFGcc_rare@sam_data$Shannon <- diversity(EFGcc_rare@otu_table)

tiff(filename="./graphics/I4512ivac_cc_shannon_groups.tiff", width=8, height=8, units="in", res=600, bg="transparent", family = "sans", pointsize = 12,  compression="lzw")
shancc <- EFGcc_rare@sam_data %>%
  ggplot(aes(x=Group, y=Shannon, group=set, fill=Group)) +
  geom_boxplot() +
  facet_wrap(~Days_post_challenge, nrow=1) +
  scale_fill_manual(values=c(Mock_Mock='green', Mock_I4512i='skyblue', Enterisol_I4512i='orange')) +
  theme_bw()
shancc + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
dev.off()

tiff(filename="./graphics/I4512ivac_cc_shannon_groups2.tiff", width=8, height=8, units="in", res=600, bg="transparent", family = "sans", pointsize = 12,  compression="lzw")
shan1cc <- EFGcc_rare@sam_data %>%
  ggplot(aes(x=Treatment, y=Shannon, group=set, fill=Treatment)) +
  geom_boxplot() +
  facet_wrap(~Days_post_challenge, nrow=1) +
  scale_fill_manual(values=c(MVMC='green', MVC='skyblue', EVC='orange')) +
  theme_bw()
shan1cc + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
dev.off()

### Simpson's Diversity
#Fecal
EFGf_rare@sam_data$Simpson <- diversity(EFGf_rare@otu_table)

tiff(filename="./graphics/I4512ivac_fecal_simpson_groups.tiff", width=8, height=8, units="in", res=600, bg="transparent", family = "sans", pointsize = 12,  compression="lzw")
simpf <- EFGf_rare@sam_data %>%
  ggplot(aes(x=Group, y=Simpson, group=set, fill=Group)) +
  geom_boxplot() +
  facet_wrap(~Days_post_challenge, nrow=1) +
  scale_fill_manual(values=c(Mock_Mock='green', Mock_I4512i='skyblue', Enterisol_I4512i='orange')) +
  theme_bw()
simpf + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
dev.off()

tiff(filename="./graphics/I4512ivac_fecal_simpson_groups2.tiff", width=8, height=8, units="in", res=600, bg="transparent", family = "sans", pointsize = 12,  compression="lzw")
simp1f <- EFGf_rare@sam_data %>%
  ggplot(aes(x=Treatment, y=Simpson, group=set, fill=Treatment)) +
  geom_boxplot() +
  facet_wrap(~Days_post_challenge, nrow=1) +
  scale_fill_manual(values=c(MVMC='green', MVC='skyblue', EVC='orange')) +
  theme_bw()
simp1f + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
dev.off()

#Cecal_contents
EFGcc_rare@sam_data$Simpson <- diversity(EFGcc_rare@otu_table)

tiff(filename="./graphics/I4512ivac_cc_simpson_groups.tiff", width=8, height=8, units="in", res=600, bg="transparent", family = "sans", pointsize = 12,  compression="lzw")
simpcc <- EFGcc_rare@sam_data %>%
  ggplot(aes(x=Group, y=Simpson, group=set, fill=Group)) +
  geom_boxplot() +
  facet_wrap(~Days_post_challenge, nrow=1) +
  scale_fill_manual(values=c(Mock_Mock='green', Mock_I4512i='skyblue', Enterisol_I4512i='orange')) +
  theme_bw()
simpcc + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
dev.off()

tiff(filename="./graphics/I4512ivac_cc_simpson_groups2.tiff", width=8, height=8, units="in", res=600, bg="transparent", family = "sans", pointsize = 12,  compression="lzw")
simp1cc <- EFGcc_rare@sam_data %>%
  ggplot(aes(x=Treatment, y=Simpson, group=set, fill=Treatment)) +
  geom_boxplot() +
  facet_wrap(~Days_post_challenge, nrow=1) +
  scale_fill_manual(values=c(MVMC='green', MVC='skyblue', EVC='orange')) +
  theme_bw()
simp1cc + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
dev.off()

### Chao1's Diversity
#Fecal
EFGf_rare@sam_data$Chao1 <- diversity(EFGf_rare@otu_table)

tiff(filename="./graphics/I4512ivac_fecal_chao1_groups.tiff", width=8, height=8, units="in", res=600, bg="transparent", family = "sans", pointsize = 12,  compression="lzw")
chaof <- EFGf_rare@sam_data %>%
  ggplot(aes(x=Group, y=Chao1, group=set, fill=Group)) +
  geom_boxplot() +
  facet_wrap(~Days_post_challenge, nrow=1) +
  scale_fill_manual(values=c(Mock_Mock='green', Mock_I4512i='skyblue', Enterisol_I4512i='orange')) +
  theme_bw()
chaof + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
dev.off()

tiff(filename="./graphics/I4512ivac_fecal_chao1_groups2.tiff", width=8, height=8, units="in", res=600, bg="transparent", family = "sans", pointsize = 12,  compression="lzw")
chao1f <- EFGf_rare@sam_data %>%
  ggplot(aes(x=Treatment, y=Chao1, group=set, fill=Treatment)) +
  geom_boxplot() +
  facet_wrap(~Days_post_challenge, nrow=1) +
  scale_fill_manual(values=c(MVMC='green', MVC='skyblue', EVC='orange')) +
  theme_bw()
chao1f + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
dev.off()

#Cecal_contents
EFGcc_rare@sam_data$Chao1 <- diversity(EFGcc_rare@otu_table)

tiff(filename="./graphics/I4512ivac_cc_shannon_groups.tiff", width=8, height=8, units="in", res=600, bg="transparent", family = "sans", pointsize = 12,  compression="lzw")
chaocc <- EFGcc_rare@sam_data %>%
  ggplot(aes(x=Group, y=Chao1, group=set, fill=Group)) +
  geom_boxplot() +
  facet_wrap(~Days_post_challenge, nrow=1) +
  scale_fill_manual(values=c(Mock_Mock='green', Mock_I4512i='skyblue', Enterisol_I4512i='orange')) +
  theme_bw()
chaocc + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
dev.off()

tiff(filename="./graphics/I4512ivac_cc_shannon_groups2.tiff", width=8, height=8, units="in", res=600, bg="transparent", family = "sans", pointsize = 12,  compression="lzw")
chao1cc <- EFGcc_rare@sam_data %>%
  ggplot(aes(x=Treatment, y=Chao1, group=set, fill=Treatment)) +
  geom_boxplot() +
  facet_wrap(~Days_post_challenge, nrow=1) +
  scale_fill_manual(values=c(MVMC='green', MVC='skyblue', EVC='orange')) +
  theme_bw()
chao1cc + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
dev.off()





# Linear model
#Fecal
mod_datf <-
  EFGf_rare@sam_data %>%
  data.frame() %>%
  mutate(Group=factor(Group, levels=c('Enterisol_I4512i', 'Mock_I4512i', 'Mock_Mock')),
         Days_post_challenge = factor(Days_post_challenge, levels = c('-28', '-1', '2', '3', '7','10', '14')))
colnames(EFGf_rare@sam_data)

shannon_modf <- lm(data = mod_datf, formula= Shannon ~ Days_post_challenge * Group)
shannon_modf %>% summary()

#Cecal_contents
mod_datcc <-
  EFGcc_rare@sam_data %>%
  data.frame() %>%
  mutate(Group=factor(Group, levels=c('Enterisol_I4512i', 'Mock_I4512i', 'Mock_Mock')),
         Days_post_challenge = factor(Days_post_challenge, levels = c('14')))
colnames(EFGcc_rare@sam_data)

shannon_modcc <- lm(data = mod_datcc, formula= Shannon ~ Group)
shannon_modcc %>% summary()
simpson_modf <- lm(data = mod_datf, formula= Simpson ~ Group)


# EMMEANS - estimated marginal means (aka least squares means) for summarizing a fitted linear model that include factors
#Fecal
library(emmeans)
EMMEANSf <- emmeans(shannon_modf, ~ Group | Days_post_challenge)
shan_contrastsf <- emmeans::contrast(EMMEANSf, method='pairwise', adjust='fdr')
contrasts_tablef <- shan_contrastsf %>% as.data.frame()
plot_datf <- EMMEANSf %>% as.data.frame()

tiff(filename="./graphics/I4512ivac_fecal_alpha_div_over_time.tiff", width=8, height=8, units="in", res=600, bg="transparent", family = "sans", pointsize = 12,  compression="lzw")
plot_datf %>%
  ggplot(aes(x=Days_post_challenge,
             y=emmean,ymin=lower.CL,
             ymax=upper.CL,
             group=Group, color=Group)) +
  geom_pointrange(position = position_dodge(width = .2)) +
  geom_line() + ggtitle('Alpha diversity over time') +
  ylab('Shannon') +
  scale_color_manual(values=c(Mock_Mock = 'green',
                              Mock_I4512i = 'skyblue',
                              Enterisol_I4512i = 'orange')) +
  theme_bw()
dev.off()

#Cecal_contents
EMMEANScc <- emmeans(shannon_modcc, ~ Group)
shan_contrastscc <- emmeans::contrast(EMMEANScc, method='pairwise', adjust='fdr')
contrasts_tablecc <- shan_contrastscc %>% as.data.frame()
plot_datcc <- EMMEANScc %>% as.data.frame()


# Alpha contrasts table
#Fecal
alpha_contrasts_tablef <-
  contrasts_tablef %>%
  mutate(across(where(is.numeric), signif, digits=2),
         sig=ifelse(p.value < 0.05, TRUE, FALSE))
write_tsv(alpha_contrasts_tablef, "./outputs/I4512ivac_fecal_alpha_tests.tsv")

#Cecal_contents
alpha_contrasts_tablecc <-
  contrasts_tablecc %>%
  mutate(across(where(is.numeric), signif, digits=2),
         sig=ifelse(p.value < 0.05, TRUE, FALSE))
write_tsv(alpha_contrasts_tablecc, "./outputs/I4512ivac_cecalcontents_alpha_tests.tsv")


# Beta dispersion
BETA_DISPERf <- betadisper(DISTf, group = sample_data(EFGf)$set)
EFGf_rare@sam_data$Distance_from_centroid <- BETA_DISPERf$distances
EFGf_rare@sam_data %>% ggplot(aes(x=Group, y=Distance_from_centroid)) + geom_boxplot() +
  facet_wrap(~Days_post_challenge, ncol=1)

BETA_DISPERcc <- betadisper(DISTcc, group = sample_data(EFGcc)$set)
EFGcc_rare@sam_data$Distance_from_centroid <- BETA_DISPERcc$distances
EFGcc_rare@sam_data %>% ggplot(aes(x=Group, y=Distance_from_centroid)) + geom_boxplot() +
  facet_wrap(~Days_post_challenge, ncol=1)


EMMEANSf <- emmeans(shannon_modf, ~ Group | Days_post_challenge)
shan_contrastsf <- emmeans::contrast(EMMEANSf, method='pairwise', adjust='fdr')
contrasts_tablef <- shan_contrastsf %>% as.data.frame()
plot_datf <- EMMEANSf %>% as.data.frame()



tiff(filename="./graphics/I4512ivac_fecal_alpha_div_over_time.tiff", width=8, height=8, units="in", res=600, bg="transparent", family = "sans", pointsize = 12,  compression="lzw")
EFGf_rare@sam_data %>%
  ggplot(aes(x=Days_post_challenge,
             y=Distance_from_centroid,
             group=Group, color=Group)) +
  geom_line() + ggtitle('Beta dispersion over time') +
  scale_color_manual(values=c(Mock_Mock = 'green',
                              Mock_I4512i = 'skyblue',
                              Enterisol_I4512i = 'orange')) +
  theme_bw()
dev.off()




#RELATIVE ABUNDANCE
table(phyloseq::tax_table(PHYLOf)[, "Phylum"]) #count of phyla
PHYLOf_rel_abund = phyloseq::transform_sample_counts(PHYLOf, function(x){x / sum(x)})

phyloseq::plot_bar(PHYLOf_rel_abund, fill = "Phylum") +
  geom_bar(aes(color = Phylum, fill = Phylum), stat = "identity", position = "stack") +
  labs(x = "", y = "Relative Abundance\n") +
  facet_wrap(~ Treatment, scales = "free") +
  theme(panel.background = element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())

relat_abund <- PHYLOf %>% tax_glom(taxrank="Phylum") %>% psmelt()
head(relat_abund)
unique(relat_abund$Phylum)


# Set colors for plotting
my_colors <- c(
  '#a6cee3','#1f78b4','#b2df8a','#33a02c','#fb9a99','#e31a1c',
           '#fdbf6f','#ff7f00','#cab2d6','#6a3d9a','#ffff99','#b15928', 
           "#CBD588", "#5F7FC7", "orange","#DA5724", "#508578", "#CD9BCD",
           "#AD6F3B", "#673770","#D14285", "#652926", "#C84248", 
           "#8569D5", "#5E738F","#D1A33D", "#8A7C64", "#599861", "black"
)

#If you want different taxonomic level, find and replace the taxonomic level listed here
#Fecal
otu.summaryf <- prop.table(as.matrix(OTUmatf), 1) 
otu_abundf <- colSums(otu.summaryf)
otu.summaryf <- rbind(otu_abundf, otu.summaryf)
otu.summary_sortedf <- otu.summaryf[,order(otu.summaryf[1,], decreasing = TRUE)]
#Cecal_contents
otu.summarycc <- prop.table(as.matrix(OTUmatcc), 1) 
otu_abundcc <- colSums(otu.summarycc)
otu.summarycc <- rbind(otu_abundcc, otu.summarycc)
otu.summary_sortedcc <- otu.summarycc[,order(otu.summarycc[1,], decreasing = TRUE)]


#top 16 most abundant genera
num_genera <- 16 # enter the number of genera you want

#melts large dataset into two columns
#Fecal
melt_otuf <- melt(otu.summary_sortedf[,c(1:num_genera)])
colnames(melt_otuf) <- c("SampleID", "OTU", "Abundance")
tail(melt_otuf)
#Cecal_contents
melt_otucc <- melt(otu.summary_sortedcc[,c(1:num_genera)])
colnames(melt_otucc) <- c("SampleID", "OTU", "Abundance")
tail(melt_otucc)

#Putting it all together: merge melt_otu, metadata, taxonomy tables
#Fecal
meta_otuf <- merge(METf, melt_otuf, by.x = "SampleID", by.y = "SampleID") #merge by SampleID
if (!file.exists("outputs/I4512ivac_fecal_metaotutax.rds")) {
  meta_otu_taxf <- merge(meta_otuf, PHYLOf@tax_table) #merge meta & otu abudunces %>% 
    meta_otu_taxf <- meta_otu_taxf[order(meta_otu_taxf$NMDS1),]
  write_rds(meta_otu_taxf, "outputs/I4512ivac_fecal_metaotutax.rds")
} else {
  meta_otu_taxf <- readRDS("outputs/I4512ivac_fecal_metaotutax.rds")
}
str(meta_otu_taxf)
summary(meta_otu_taxf$Treatment)
#Cecal_contents
meta_otucc <- merge(METcc, melt_otucc, by.x = "SampleID", by.y = "SampleID") #merge by SampleID
if (!file.exists("outputs/I4512ivac_cecalcontents_metaotutax.rds")) {
  meta_otu_taxcc <- merge(meta_otucc, PHYLOcc@tax_table) #merge meta & otu abudunces %>% 
    meta_otu_taxcc <- meta_otu_taxcc[order(meta_otu_taxcc$NMDS1),] 
  write_rds(meta_otu_taxcc, "outputs/I4512ivac_cecalcontents_metaotutax.rds")
} else {
  meta_otu_taxcc <- readRDS("outputs/I4512ivac_cecalcontents_metaotutax.rds")
}
str(meta_otu_taxcc)
summary(meta_otu_taxcc$Treatment)

#sorting based on MDS1 from negative to positive (NMDS axis 1) and save as rds
#Fecal
meta_otu_taxf <- meta_otu_taxf[order(meta_otu_taxf$NMDS1),]
if (!file.exists("outputs/I4512ivac_fecal_metaotutax.rds")) {
  meta_otu_taxf <- meta_otu_taxf[order(meta_otu_taxf$NMDS1),] 
  write_rds(meta_otu_taxf, "outputs/I4512ivac_fecal_metaotutax.rds")
} else {
  meta_otu_taxf <- readRDS("outputs/I4512ivac_fecal_metaotutax.rds")
}
#Cecal_contents
meta_otu_taxcc <- meta_otu_taxcc[order(meta_otu_taxcc$NMDS1),]
if (!file.exists("outputs/I4512ivac_cecalcontents_metaotutax.rds")) {
  meta_otu_taxcc <- meta_otu_taxcc[order(meta_otu_taxcc$NMDS1),] 
  write_rds(meta_otu_taxcc, "outputs/I4512ivac_cecalcontents_metaotutax.rds")
} else {
  meta_otu_taxcc <- readRDS("outputs/I4512ivac_cecalcontents_metaotutax.rds")
}
#ordering samples based on NMDS axis 1
meta_otu_tax$Treatment <- factor(meta_otu_tax$Treatment, levels=unique(as.character(meta_otu_tax$Treatment)) )

#MAKE A GRAPH! Plot individuals not group means
ggplot(meta_otu_tax, aes(x = Treatment, y = Abundance, fill = genus)) + 
  geom_bar(stat = "identity") +
  scale_fill_manual(values = my_colors) +
  # Remove x axis title
  theme(axis.title.x = element_blank()) + 
  ylim(c(0,1)) +
  guides(fill = guide_legend(reverse = F, keywidth = .5, keyheight = .5, ncol = 1)) +
  theme(legend.text=element_text(size=8)) +
  #theme(legend.position="bottom") +
  theme(axis.text.x = element_blank()) +
  ylab(paste0("Relative Abundance (top ", num_genera, " genera)")) +
  ggtitle("Genus Composition by Bull sorted by NMDS Axis 1") 
ggsave("GenusBarPlotNMDS1_Phylo_pwpt.png", width = 7, height = 3)







tmp <- NMDS_from_phyloseq(EFGf, 'set')

tmp[[2]] <-
  tmp[[2]] %>%
  mutate(Days_post_challenge=sub('(-?[0-9]+)_([A-Za-z0-9]+)_([A-Za-z]+)','\\1',group),
         Vaccine=sub('(-?[0-9]+)_([A-Za-z0-9]+)_([A-Za-z]+)','\\2',group),
         Challenge=sub('(-?[0-9]+)_([A-Za-z0-9]+)_([A-Za-z]+)','\\3',group),
         GROUP=paste(Vaccine, Challenge, sep='_'))







