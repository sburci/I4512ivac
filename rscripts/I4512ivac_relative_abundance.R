# I 4,[5],12:i:- MICROBIOME ANALYSIS
# Selma Burciaga
# March 20, 2023
# Rscript: I4512ivac_relative_abundance.R
# R v4.2.2/R v2022.12.0

# PART 4: 
## anything ending with f is for Fecal samples
## anything ending with f2 is for Fecal samples without MVMC treatment
## anything ending with cc is for Cecal Content samples



# RELATIVE ABUNDANCE: Richness & Evenness___________________________________________________
library(phyloseq) #v1.42.0
library(tidyverse) #v1.3.2
library(vegan) #v2.6.4
library(ggplot2) #3.4.0
library(dplyr)
library(readr)
library(reshape2)

## Open I4512ivac_phylo_all.rds and create phylo subset objects
PHYLO <- read_rds("./outputs/I4512ivac_phylo_all.rds")
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

# Relative abundance: PHYLUM
table(phyloseq::tax_table(PHYLOf)[, "Phylum"]) #count of phyla
PHYLOf_rel_abund = phyloseq::transform_sample_counts(PHYLOf, function(x){x / sum(x)})

phyloseq::plot_bar(PHYLOf_rel_abund, fill = "Phylum") +
  geom_bar(aes(color = Phylum, fill = Phylum), stat = "identity", position = "stack") +
  labs(x = "", y = "Relative Abundance\n") +
  facet_wrap(~ Treatment, scales = "free") +
  theme(panel.background = element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())

tiff(filename="./graphics/I4512ivac_fecal_relative_abundance.tiff", width=8, height=8, units="in", res=600, bg="transparent", family = "sans", pointsize = 12,  compression="lzw")
PHYLOf_rel_abund_graph <- 
  phyloseq::plot_bar(PHYLOf_rel_abund, fill = "Phylum") +
  geom_bar(aes(color = Phylum, fill = Phylum), stat = "identity", position = "stack") +
  labs(x = "", y = "Relative Abundance\n") +
  facet_wrap(~ Treatment, scales = "free") +
  theme(panel.background = element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())
dev.off()

relat_abundf <- PHYLOf %>% tax_glom(taxrank="Phylum") %>% psmelt()
head(relat_abundf)
unique(relat_abundf$Phylum)


table(phyloseq::tax_table(PHYLOcc)[, "Phylum"]) #count of phyla
PHYLOcc_rel_abund = phyloseq::transform_sample_counts(PHYLOcc, function(x){x / sum(x)})

phyloseq::plot_bar(PHYLOcc_rel_abund, fill = "Phylum") +
  geom_bar(aes(color = Phylum, fill = Phylum), stat = "identity", position = "stack") +
  labs(x = "", y = "Relative Abundance\n") +
  facet_wrap(~ Treatment, scales = "free") +
  theme(panel.background = element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())

relat_abundcc <- PHYLOcc %>% tax_glom(taxrank="Phylum") %>% psmelt()
head(relat_abundcc)
unique(relat_abundcc$Phylum)

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
