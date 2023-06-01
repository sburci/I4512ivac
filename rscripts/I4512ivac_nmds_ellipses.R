# I 4,[5],12:i:- MICROBIOME ANALYSIS
# Selma Burciaga
# March 20, 2023
# Rscript: I4512ivac_nmds.R
# R v4.2.2/R v2022.12.0

# PART 1: Using the vegan R package to calculate ecological distances (CC188) - https://riffomonas.org/code_club/2022-02-17-distances
## anything ending with f is for Fecal samples
## anything ending with f2 is for Fecal samples without MVMC treatment
## anything ending with cc is for Cecal Content samples


# BETA DIVERSITY: NMDS WITH ELLIPSES___________________________________________________
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
PHYLOcc <- prune_samples(samples = PHYLO@sam_data$Sample_Type=="Cecal_contents", x = PHYLO)

## Create OTU R matrix and dataframe with metadata from phyloseq object
OTUmatf <- as(otu_table(PHYLOf), "matrix")
METf <- sample_data(PHYLOf) %>% 
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
METcc <- METcc %>% 
  group_by(dpc, Treatment) %>% 
  dplyr::mutate(centnmds1=mean(NMDS1),
         centnmds2=mean(NMDS2)) %>% 
  ungroup()


###By days post challenge
#Fecal
tiff(filename="./graphics/I4512ivac_fecal_alltreatmentsbydpc_nmds_ell.tiff", width=8, height=8, units="in", res=600, bg="transparent", family = "sans", pointsize = 12,  compression="lzw")
METf %>% 
  ggplot(aes(x=NMDS1, y=NMDS2, color=Treatment)) +
  geom_point() +
  ggtitle('NMDS of all treatment groups over time by days post challenge') +
  scale_color_manual(values=c(MVMC='green', MVC='skyblue', EVC='orange')) +
  facet_wrap(~Days_post_challenge, scales = "free") +
  stat_ellipse() + #Plots the ellipses
  theme_bw()
dev.off()

tiff(filename="./graphics/I4512ivac_fecal_alltreatments_nmds_centroids.tiff", width=8, height=8, units="in", res=600, bg="transparent", family = "sans", pointsize = 12,  compression="lzw")
METf %>% 
  unique() %>% 
  arrange(dpc) %>% 
  ggplot(aes(x=centnmds1, y=centnmds2, color=Treatment)) +
  geom_point(data=METf, aes(x=NMDS1, y=NMDS2, color=Treatment), alpha=.5) +
  stat_ellipse() +
  geom_path(aes(group=Treatment)) +
  geom_point(size=6) +
  geom_text(aes(label=dpc), color='black') +
  ggtitle('NMDS of all treatment groups over time') +
  xlab('NMDS1') +
  ylab('NMDS2') +
  scale_color_manual(values=c(MVMC='green', MVC='skyblue', EVC='orange')) +
  theme_bw()
dev.off()

#Cecal_contents
tiff(filename="./graphics/I4512ivac_cecalcontents_alltreatments_nmds_centroids.tiff", width=8, height=8, units="in", res=600, bg="transparent", family = "sans", pointsize = 12,  compression="lzw")
METcc %>% 
  dplyr::select(centnmds1, centnmds2, Treatment, dpc) %>% 
  unique() %>% 
  arrange(dpc) %>% 
  ggplot(aes(x=centnmds1, y=centnmds2, color=Treatment)) +
  geom_point(data=METcc, aes(x=NMDS1, y=NMDS2, color=Treatment), alpha=.5) +
  geom_segment(data=METcc, aes(x=NMDS1, y=NMDS2, xend=centnmds1, yend=centnmds2, color=Treatment), alpha=.25) +
  geom_path(aes(group=Treatment)) +
  geom_point(size=6) +
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
