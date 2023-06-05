# I 4,[5],12:i:- MICROBIOME ANALYSIS
# Selma Burciaga
# May 25, 2023
# Rscript: I4512ivac_relative_abundance.R
# R v4.2.2/R v2022.12.0

## PART 4: 
## anything ending with f is for Fecal samples
## anything ending with cc is for Cecal Content samples
## code in this script was created with the help of David J. Bradshaw

## Note: Save as graphics as PDF - for presentation (4x6) or manuscript
##       Center all plot titles - theme_update(plot.title = element_text(hjust = 0.5))


## RELATIVE ABUNDANCE: % OTU detected___________________________________________________
library(phyloseq)
library(tidyverse)
library(vegan)
library(ggplot2)
library(dplyr)
library(readr)
library(reshape2)
library(plyr)
library(Rmisc)
library(RColorBrewer)

## Open I4512ivac_phylo_all.rds and create phylo subset objects
PHYLO <- readRDS("./outputs/I4512ivac_phylo_all.rds")
PHYLOf <- prune_samples(samples = PHYLO@sam_data$Sample_Type=="Fecal", x = PHYLO)
PHYLOcc <- prune_samples(samples = PHYLO@sam_data$Sample_Type=="Cecal_contents", x = PHYLO)

## View metadata for reference
MET <- sample_data(PHYLO) %>% 
  as(Class="data.frame")

## Read the taxonomy file into R
TAX <- read.table("./mothuro/I4512ivac.otu.taxonomy", sep = "\t", header = T)
TAX <- separate(data = TAX, col = Taxonomy, into = c("Kingdom", "Phylum", "Class", "Family", "Order", "Genus", "Species"), sep = ";")
str(TAX)

## Filter taxanomy - total # seq associated with taxa (x) must be over 10
#data = filter_taxa(data, function(x) sum(x) >10, TRUE)
#ntaxa(data)

#_____________________________________________________________________________________________________________________________________________________________

## RELATIVE ABUNDANCE (ra)
raf = phyloseq::transform_sample_counts(PHYLOf, function(x)100* x / sum(x))
racc = phyloseq::transform_sample_counts(PHYLOcc, function(x)100* x / sum(x))

#GENRA
#Fecal
genraf <- tax_glom(raf, taxrank="Genus") #agglomerate taxa
dataf <- psmelt(genraf) #create df from phyloseq object
#dataf <- unite(dataf, Genus, Family:Genus, sep=";") #optional: combine family and genus columns
dataf$Genus <- as.character(dataf$Genus) #convert genus to charactor vector from factor (b/c R)
meansf <- ddply(dataf, ~Genus, function(x) c(Mean=mean(x$Abundance))) #group df by genus, calculate mean ra
Other <- meansf[meansf$Mean <=1,]$Genus #find phyla whose ra is < than 1%
dataf[dataf$Genus %in% Other,]$Genus <- 'Other Prokaryotes' #change their name to "Other Prokaryotes"
dataf <- dataf[!dataf$Genus=="Other Prokaryotes",] #remove all genera labeled "Other Prokaryotes"
dataf <- subset(dataf, select=c(SampleID, Treatment, dpc, tdpc, Abundance, Genus)) #remove unnecessary columns
dataf <- summarySE(data = dataf, measurevar = "Abundance", groupvars = c("dpc", "Treatment", "Genus"), na.rm = TRUE) #summarize based on target parameter
dataf <- arrange(dataf, Abundance) #arrange by abundance
dataf$dpc_Treatment <- paste(dataf$dpc, dataf$Treatment, sep = "_") #add column combining target columns (x and faceting)
Abundance <- ddply(dataf, ~dpc_Treatment, function(x) c(Abundance=100-sum(x$Abundance))) #create table (Abundance) that is the leftover Prokaryotes
Abundance$Genus <- "Other Prokaryotes" #add column to Abundance labeling leftover Prokaryotes
Abundance <- separate(data = Abundance, col = dpc_Treatment, into = c("dpc", "Treatment"), sep = "_", remove=FALSE) #split the combined column in Abundance into separate columns but keep the original
dataf <- subset(dataf, select=c(dpc, dpc_Treatment, Abundance, Treatment, Genus)) #remove unnecessary columns from dataf
colnames(dataf) #check column names of both dataf and Abundance & change order of Abundance to match dataf
colnames(Abundance) #[1] "dpc" "dpc_Treatment" "Abundance" "Treatment" "Genus"
Abundance <- Abundance[,c(2,1,4,3,5)] #[1] "dpc_Treatment" "dpc" "Treatment""Abundance" "Genus"
colnames(Abundance) #confirm column names match order of dataf
dataf <- rbind(dataf, Abundance) #combine with original table

#Save relative abundance (%)
write_tsv(dataf, "./outputs/I4512ivac_rel_abund_percentages.tsv")

#Create relative abundance stacked bar graphs
tiff('./graphics/I4512ivac_fecal_relative_abundance_percent_bytreatment.tiff', units="in", width=10, height=6, res=600, family = "sans", pointsize = 12,  compression="lzw")
ggplot(dataf, aes(x = dpc, y = Abundance, fill = Genus)) +
  geom_bar(stat = "identity", color = "black") +
  facet_grid(~Treatment) + 
  labs(x = "Days post challenge", y = "Relative abundance (%)") +
  ggtitle("Relative abundance of fecal samples by treatment") +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    panel.background = element_rect(fill = "transparent",colour = NA),
    plot.background = element_rect(fill = "transparent",colour = NA)) +
  scale_fill_manual("Genus", values=c("#B2182B", "#D6604D", "#F4A582", "#FDDBC7",
                                      "#C51B7D", "#DE77AE", "#F1B6DA", "#FDE0EF", 
                                      "#B35806", "#E08214", "#FDB863", "#FEE0B6",
                                      "#4D9221", "#7FBC41", "#B8E186", "#E6F5D0", 
                                      "#2166AC", "#4393C3", "#92C5DE", "#D1E5F0", 
                                      "#542788", "#8073AC", "#B2ABD2", "#D8DAEB"))
dev.off()

tiff('./graphics/I4512ivac_fecal_relative_abundance_percent_bydpi.tiff', units="in", width=10, height=6, res=600, family = "sans", pointsize = 12,  compression="lzw")
ggplot(dataf, aes(x = Treatment, y = Abundance, fill = Genus)) +
  geom_bar(stat = "identity", color = "black") +
  facet_grid(~dpc) + 
  labs(x = "Treatment", y = "Relative abundance (%)") +
  ggtitle("Relative abundance of fecal samples by days post challenge") +
  theme(plot.title = element_text(hjust = 0.5)) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  theme(
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    panel.background = element_rect(fill = "transparent",colour = NA),
    plot.background = element_rect(fill = "transparent",colour = NA)) +
  scale_fill_manual("Genus", values=c("#B2182B", "#D6604D", "#F4A582", "#FDDBC7",
                                      "#C51B7D", "#DE77AE", "#F1B6DA", "#FDE0EF", 
                                      "#B35806", "#E08214", "#FDB863", "#FEE0B6",
                                      "#4D9221", "#7FBC41", "#B8E186", "#E6F5D0", 
                                      "#2166AC", "#4393C3", "#92C5DE", "#D1E5F0", 
                                      "#542788", "#8073AC", "#B2ABD2", "#D8DAEB"))
dev.off()

#Cecal_contents
genracc <- tax_glom(racc, taxrank="Genus") #agglomerate taxa
datacc <- psmelt(genracc) #create df from phyloseq object
#datacc <- unite(datacc, Genus, Family:Genus, sep=";") #optional: combine family and genus columns
datacc$Genus <- as.character(datacc$Genus) #convert genus to charactor vector from factor (b/c R)
meanscc <- ddply(datacc, ~Genus, function(x) c(Mean=mean(x$Abundance))) #group df by genus, calculate mean ra
Other <- meanscc[meanscc$Mean <=1,]$Genus #find phyla whose ra is < than 1%
datacc[datacc$Genus %in% Other,]$Genus <- 'Other Prokaryotes' #change their name to "Other Prokaryotes"
datacc <- datacc[!datacc$Genus=="Other Prokaryotes",] #remove all genera labeled "Other Prokaryotes"
datacc <- subset(datacc, select=c(SampleID, Treatment, dpc, tdpc, Abundance, Genus)) #remove unnecessary columns
datacc <- summarySE(data = datacc, measurevar = "Abundance", groupvars = c("dpc", "Treatment", "Genus"), na.rm = TRUE) #summarize based on target parameter
datacc <- arrange(datacc, Abundance) #arrange by abundance
datacc$dpc_Treatment <- paste(datacc$dpc, datacc$Treatment, sep = "_") #add column combining target columns (x and faceting)
Abundance <- ddply(datacc, ~dpc_Treatment, function(x) c(Abundance=100-sum(x$Abundance))) #create table (Abundance) that is the leftover Prokaryotes
Abundance$Genus <- "Other Prokaryotes" #add column to Abundance labeling leftover Prokaryotes
Abundance <- separate(data = Abundance, col = dpc_Treatment, into = c("dpc", "Treatment"), sep = "_", remove=FALSE) #split the combined column in Abundance into separate columns but keep the original
datacc <- subset(datacc, select=c(dpc, dpc_Treatment, Abundance, Treatment, Genus)) #remove unnecessary columns from dataf
colnames(datacc) #check column names of both dataf and Abundance & change order of Abundance to match dataf
colnames(Abundance) #[1] "dpc" "dpc_Treatment" "Abundance" "Treatment" "Genus"
Abundance <- Abundance[,c(2,1,4,3,5)] #[1] "dpc_Treatment" "dpc" "Treatment""Abundance" "Genus"
colnames(Abundance) #confirm column names match order of dataf
datacc <- rbind(datacc, Abundance) #combine with original table

tiff('./graphics/I4512ivac_cecalcontents_relative_abundance_percent_bytreatment.tiff', units="in", width=10, height=6, res=600, family = "sans", pointsize = 12,  compression="lzw")
ggplot(datacc, aes(x = dpc, y = Abundance, fill = Genus)) +
  geom_bar(stat = "identity", color = "black", width = 0.20) +
  facet_grid(~Treatment) + 
  labs(x = "Days post challenge", y = "Relative abundance (%)") +
  ggtitle("Relative abundance of cecal contents by treatment") +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    panel.background = element_rect(fill = "transparent",colour = NA),
    plot.background = element_rect(fill = "transparent",colour = NA)) +
  scale_fill_manual("Genus", values=c("#e31a1c", "#B2182B", "#D6604D", "#F4A582", 
                                      "#C51B7D", "#FDE0EF", "#B35806", "#FEE0B6",
                                      "#4D9221", "#7FBC41", "#B8E186", "#E6F5D0",
                                      "#2166AC", "#4393C3", "#FFFF00", "#FFFACD",
                                      "#92C5DE", "#D1E5F0",
                                      "#542788", "#8073AC", "#B2ABD2", "#D8DAEB"))
dev.off()

tiff('./graphics/I4512ivac_cecalcontents_relative_abundance_percent_bydpi.tiff', units="in", width=10, height=6, res=600, family = "sans", pointsize = 12,  compression="lzw")
ggplot(datacc, aes(x = Treatment, y = Abundance, fill = Genus)) +
  geom_bar(stat = "identity", color = "black", width = 0.15, position = "fill") +
  facet_grid(~dpc) + 
  labs(x = "Treatment", y = "Relative abundance (%)") +
  ggtitle("Relative abundance of cecal contents by days post challenge") +
  theme(plot.title = element_text(hjust = 0.5)) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  theme(
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    panel.background = element_rect(fill = "transparent",colour = NA),
    plot.background = element_rect(fill = "transparent",colour = NA)) +
  scale_fill_manual("Genus", values=c("#e31a1c", "#B2182B", "#D6604D", "#F4A582", 
                                      "#C51B7D", "#FDE0EF", "#B35806", "#FEE0B6",
                                      "#4D9221", "#7FBC41", "#B8E186", "#E6F5D0",
                                      "#2166AC", "#4393C3", "#FFFF00", "#FFFACD",
                                      "#92C5DE", "#D1E5F0",
                                      "#542788", "#8073AC", "#B2ABD2", "#D8DAEB"))
dev.off()











#_____________________________________________________________________________________________________________________________________________________________


##Other commands

# Top 10 genra
genra10 <- names(sort(taxa_sums(PHYLOf_genra), TRUE)[1:10])
PHYLOf_genra10 <- prune_taxa(genra10, PHYLOf_genra)

dffra <- psmelt(PHYLOf_genra10) %>%
  group_by(Treatment, dpc, tdpc)

ggplot(dffra, aes(x = dpc, y = Abundance, fill = Genus)) +
  geom_bar(stat = "identity") +
  facet_grid(~Treatment) + 
  labs(x = "Days post challenge", y = "Relative abundance")

genus10 <- names(sort(taxa_sums(PHYLOf_genus), TRUE)[1:10])
PHYLOf_genus10 <- prune_taxa(genus10, PHYLOf_genus)
PHYLOf_genus10_rel_abund = phyloseq::transform_sample_counts(PHYLOf_genus10, function(x){x / sum(x)})
PHYLOf_genus10_rel_abund2 = phyloseq::transform_sample_counts(PHYLOf_genus10, function(x)100* x / sum(x))

df <- psmelt(PHYLOf_genra10) %>%
  filter(!is.na(Genus)) %>%
  group_by(Treatment, dpc, tdpc, Genus) %>%
  # Add abundances within each genus and sample
  summarize_at("Abundance", sum)

df4 <- psmelt(PHYLOf_genus10) %>%
  filter(!is.na(Genus)) %>%
  group_by(Treatment, dpc, tdpc, Genus) %>%
  # Add abundances within each genus and sample
  summarize_at("Abundance", sum) %>% 
  mutate(Percent = 100* Abundance/sum(Abundance))

ggplot(df, aes(x = dpc, y = Abundance, fill = Genus)) +
  geom_bar(stat = "identity") +
  facet_grid(~Treatment) + 
  labs(x = "Days post challenge", y = "Relative abundance")

ggplot(df, aes(x = Treatment, y = Abundance, fill = Genus)) +
  geom_bar(stat = "identity") +
  facet_grid(~dpc) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

ggplot(df4, aes(x = dpc, y = Percent, fill = Genus)) +
  geom_bar(stat = "identity") +
  facet_grid(~Treatment) + 
  labs(x = "Days post challenge")

ggplot(df4, aes(x = Treatment, y = Percent, fill = Genus)) +
  geom_bar(stat = "identity") +
  facet_grid(~dpc) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))


#Cecal_contents
genus10 <- names(sort(taxa_sums(PHYLOcc_genus), TRUE)[1:10])
PHYLOcc_genus10 <- prune_taxa(genus10, PHYLOcc_genus)
PHYLOcc_genus10_rel_abund = phyloseq::transform_sample_counts(PHYLOcc_genus10, function(x){x / sum(x)})

datacc2 <- psmelt(PHYLOcc_genus10) %>%
  filter(!is.na(Genus)) %>%
  group_by(Treatment, dpc, tdpc, Genus) %>%
  # Add abundances within each genus and sample
  summarize_at("Abundance", sum)

ggplot(datacc2, aes(x = Treatment, y = Abundance, fill = Genus)) +
  geom_bar(stat = "identity") +
  facet_grid(~dpc) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))































# Phylum
##Fecal
table(phyloseq::tax_table(PHYLOf)[, "Phylum"]) #count of phyla
relat_abund_phylum_f <- PHYLOf %>% tax_glom(taxrank="Phylum") %>% psmelt()
head(relat_abund_phylum_f)
unique(relat_abund_phylum_f$Phylum)

tiff(filename="./graphics/I4512ivac_fecal_relative_abundance_phylum")
phyloseq::plot_bar(PHYLOf_rel_abund, fill = "Phylum") +
  geom_bar(aes(color = Phylum, fill = Phylum), stat = "identity", position = "stack") +
  labs(x = "", y = "Relative Abundance\n") +
  facet_wrap(~ Treatment, scales = "free") +
  theme(panel.background = element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())
dev.off()
##Cecal_contents
table(phyloseq::tax_table(PHYLOcc)[, "Phylum"]) #count of phyla
relat_abund_phylum_cc <- PHYLOcc %>% tax_glom(taxrank="Phylum") %>% psmelt()
head(relat_abund_phylum_cc)
unique(relat_abund_phylum_cc$Phylum)

tiff(filename="./graphics/I4512ivac_cecalcontents_relative_abundance_phylum")
phyloseq::plot_bar(PHYLOcc_rel_abund, fill = "Phylum") +
  geom_bar(aes(color = Phylum, fill = Phylum), stat = "identity", position = "stack") +
  labs(x = "", y = "Relative Abundance\n") +
  facet_wrap(~ Treatment, scales = "free") +
  theme(panel.background = element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())
dev.off()


# Genera
#Fecal
table(phyloseq::tax_table(PHYLOf)[, "Genus"]) #count of genera
relat_abund_genus_f <- PHYLOf %>% tax_glom(taxrank="Genus") %>% psmelt() #all
head(relat_abund_genus_f)
unique(relat_abund_genus_f$Genus)

pdf(file="./graphics/I4512ivac_fecal_relative_abundance_percentages_genra.pdf", 
    width = 6, height = 4)
phyloseq::plot_bar(PHYLOf_genus10_rel_abund, fill = "Genus") +
  geom_bar(aes(color = Genus, fill = Genus), stat = "identity", position = "stack") +
  labs(x = "", y = "Relative Abundance\n") +
  facet_wrap(~ Treatment, scales = "free") +
  ggtitle('Top 10 genera overall of fecal contents by treatment') +
  theme(panel.background = element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        strip.background = element_blank())
dev.off()

##Cecal_contents
table(phyloseq::tax_table(PHYLOcc)[, "Genus"]) #count of phyla
relat_abund_genus_cc <- PHYLOcc %>% tax_glom(taxrank="Genus") %>% psmelt()
head(relat_abund_genus_cc)
unique(relat_abund_genus_cc$Genus)

pdf(file="./graphics/I4512ivac_cecalcontents_relative_abundance_percentages_genra.pdf", 
    width = 6, height = 4)
phyloseq::plot_bar(PHYLOcc_genus10_rel_abund, fill = "Genus") +
  geom_bar(aes(color = Genus, fill = Genus), stat = "identity", position = "stack") +
  labs(x = "", y = "Relative Abundance\n") +
  facet_wrap(~ Treatment, scales = "free") +
  ggtitle('Top 10 genera of cecal contents by treatment \n 14 dpc') +
  theme(panel.background = element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        strip.background = element_blank())
dev.off()
