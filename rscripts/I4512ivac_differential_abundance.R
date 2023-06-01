# I 4,[5],12:i:- MICROBIOME ANALYSIS
# Selma Burciaga
# May 25, 2023
# Rscript: I4512ivac_differential_abundance.R
# R v4.3.0

## PART 4:
## anything ending with f is for Fecal samples
## anything ending with cc is for Cecal Content samples
## code in this script was created with the help of Carmen Wickware

##before starting script, run following commands on mothur
# count.groups(shared=I4512ivac.otu.shared) #size of smallest group: 3661.
# sub.sample(shared=I4512ivac.otu.shared, size=3661) #subsample for future analyses; output file I4512ivac.otu.0.03.subsample.shared
# count.groups(shared=I4512ivac.otu.0.03.subsample.shared) #confrim all are same size
# get.groups(shared=I4512ivac.otu.0.03.subsample.shared, design=mouse.time.design, sets=Fecal)
# lefse(shared=I4512ivac.otu.0.03.subsample.shared, design=/project/fsepru113/sburciaga/I4512ivac/mothuro/I4512ivac_design.txt)

# DIFFERENTIAL ABUNDANCE: STATS USING LEFSE___________________________________________________

#LefSe

library(MASS)
library(tidyr)
library(ggplot2)
library(dplyr)

#mothur output files into dataframe (.lefse_summary & .taxonomy)
summary_lefse <- "./mothuro/X.subsample.1.lefse_summary"
taxfile <- "./mothuro/I4512ivac.otu.taxonomy"

#read in lefse data
lefse_out <- read.table(file=summary_lefse, sep="\t", header=T)

#filter out unnwanted data - water samples, those without LDA, and p-values > 0.0001
lefse_noH2O <- lefse_out[lefse_out$Class != "water",]
lefse_NMDS <- lefse_noH2O[lefse_noH2O$Class != "-",]
lefse_NMDS_filt <- lefse_NMDS[which(lefse_NMDS$pValue <= 0.0001),]

#read in the taxonomy file and separate the taxa into columns
tax <- read.table(file=taxfile, sep="\t", header=T)
tax <- separate(data = tax, col = Taxonomy, into = c("Kingdom", "Phylum", "Class", "Family", "Order", "Genus", "Species"), sep = ";")

#merge the lefse and taxonomy by OTU
lefse_tax <- merge(lefse_NMDS_filt, tax, by.x="OTU", by.y="OTU")

#order the data.frame by logLDA
lefse_tax$genus<- reorder(lefse_tax$genus, lefse_tax$LDA)
str(lefse_tax)
lefse_tax$Class <- factor(lefse_tax$Class, levels = c("LD", "HD"))
diversity <- c("Low Diversity", "High Diversity")
names(diversity) <- c("LD", "HD")


#plot
#pdf(file="Bull_Lefse.pdf", width = 6, height = 3)
ggplot(lefse_tax, aes(x = genus, y = LDA, fill = Class)) +
  geom_bar(stat = "identity") +
  facet_grid(Class~., labeller = labeller(Class=diversity)) +
  #scale_fill_manual(values = my_colors) +
  # Remove x axis title
  #theme(axis.title.x = element_blank()) +
  ylab(paste0("Linear Discriminant Analysis Effect Size Score"))+
  #ylim(c(0,1)) +
  #guides(fill = guide_legend(reverse = F, keywidth = .5, keyheight = .5)) +
  #theme(legend.text=element_blank()) +
  theme(legend.position="none") +
  theme(axis.text.x = element_text(hjust = 1, vjust = 0.5, size = 10)) +
  coord_flip()
ggsave(file="lefse_.png", width=6, height=6)
#dev.off()







#ANCOMBC2

knitr::opts_chunk$set(message = FALSE, warning = FALSE, comment = NA, 
                      fig.width = 6.25, fig.height = 5)
library(ANCOMBC)
library(tidyverse)
library(caret)
library(DT)

## Open I4512ivac_phylo_all.rds and create phylo subset objects
PHYLO <- readRDS("./outputs/I4512ivac_phylo_all.rds")
PHYLOf <- prune_samples(samples = PHYLO@sam_data$Sample_Type=="Fecal", x = PHYLO)
PHYLOcc <- prune_samples(samples = PHYLO@sam_data$Sample_Type=="Cecal_contents", x = PHYLO)

## View metadata for reference
MET <- sample_data(PHYLO) %>% 
  as(Class="data.frame")

set.seed(123)

output = ancombc2(data = PHYLOf, tax_level = "Genus",
                  fix_formula = "dpc + Treatment",
                  rand_formula = "(dpc | SampleID)",
                  group = "Treatment",
                  dunnet = FALSE)

res_prim = output$res %>%
  mutate_if(is.numeric, function(x) round(x, 2))
res_prim %>%
  datatable(caption = "ANCOM-BC2 Primary Results")