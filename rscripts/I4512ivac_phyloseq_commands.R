# I 4,[5],12:i:- MICROBIOME ANALYSIS
# Selma Burciaga
# December 29, 2022
# Rscript: I4512ivac_phyloseq_commands.R
# R v4.2.2/RStudio v1.4.1717

# Preparing matrix of sample by taxa counts by working with mothur output files

## To use package: library(package), packageVersion("package")
library(tidyverse) #v1.3.2
library(phyloseq) #v1.42.0
library(vegan) #v2.6.4
library(funfuns) #v0.1.2
library(broom) #v1.0.1
library(DESeq2) #v1.38.1
library(cowplot) #v1.1.1
library(forcats) #v0.5.2
library(ggplot2) #v3.4.0

## Build phyloseq objects
OTU <- phyloseq::import_mothur(mothur_shared_file = './mothuro/I4512ivac.otu.shared') %>% t()
TAX <- phyloseq::import_mothur(mothur_constaxonomy_file = './mothuro/I4512ivac.otu.taxonomy')
colnames(TAX) <- c('Domain', 'Phylum', 'Class', 'Order', 'Family' , 'Genus' )

MET <- read_csv('./I4512ivac_metadata.csv') %>% 
  mutate(ID=SampleID, 
         dpc=factor(Days_post_challenge, levels=c("-28", "-1", "2", "3", "7", "10", "14"))) %>% #creating categorical variable dpc and setting order of days
  select(ID, everything()) %>% 
  column_to_rownames(var = "ID") %>% 
  sample_data()

PHYLO <- phyloseq(MET, TAX, OTU)

## Distribution of sampling sequencing depth (to remove samples with low counts)
#min(sample_sums(PHYLO))
#tiff(filename = "I4512ivac_sample_sums_histo", width =4, height=4, unit = 'in', bg="transparent", res = 600,  family = "sans", pointsize = 12,  compression="lzw") 
#hist(sample_sums(PHYLO), breaks=100)
#dev.off()

## Determine counts in each otu, if low samples check if all from same group
PHYLO@sam_data$num_counts <- PHYLO %>% sample_sums()
PHYLO@sam_data %>% ggplot(aes(x=num_counts, fill=Vaccine)) + geom_histogram()

## Determine counts in each taxonomy group 
#taxa_sums(PHYLO)
#Keep OTU with >5 counts = TRUE, if FALSE remove, command to remove taxa with low OTU count:
#PHYLO <- prune_taxa(taxa_sums(PHYLO) > 5, PHYLO)

## Save phyloseq object containing all info (metadata, taxonomy, otus) for use in vegan
write_rds(PHYLO, "outputs/I4512ivac_phylo_all.rds") #create rds file (as single object)


### PROCEED TO ALPHA (VEGAN/PAIRWISEADONIS) AND BETA DIVERSITY ANALYSES 



## Other useful commands
rank_names(TAX)
sample_variable(MET)
otu_table(OTU)[1:5, 1:5]
tax_table(TAX)[1:5, 1:4]
taxa_names(TAX)[1:10]
MET$Days_post_challenge %>% 
  unique() #prints values in column of df
