# I 4,[5],12:i:- MICROBIOME ANALYSIS
# Selma Burciaga
# March 20, 2023
# Rscript: I4512ivac_permanova.R
# R v4.2.2/R v2022.12.0

# PART 2: Using the pairwiseAdonis GitHub package to run PERMANOVA pairwise comparisons tests - https://github.com/pmartinezarbizu/pairwiseAdonis
## anything ending with f is for Fecal samples
## anything ending with f2 is for Fecal samples without MVMC treatment
## anything ending with cc is for Cecal Content samples


# BETA DIVERSITY: PAIRWISEADONIS_________________________________________________________
library(pairwiseAdonis) #0.4.1
library(RColorBrewer)
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

## Create new column combining treatment and dpc
METf <- METf %>%  
  mutate(set=paste(Treatment, dpc, sep="_"))
METf2 <- METf2 %>% 
  mutate(set=paste(Treatment, dpc, sep="_"))
METcc <- METcc %>% 
  mutate(set=paste(Treatment, dpc, sep="_"))

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

## Confirm SampleID order matches in DIST and MET objects
all(attributes(DISTf)$Labels == METf$SampleID)
all(attributes(DISTcc)$Labels == METcc$SampleID)
all(attributes(DISTf2)$Labels == METf2$SampleID)

## Running pairwiseAdonis
if (!file.exists('./outputs/I4512ivac_fecal_all_pwadon.rds')){
  PWADONf <- pairwise.adonis(x = DISTf, factors = METf$set, perm =9999)
  write_rds(PWADONf, './outputs/I4512ivac_fecal_all_pwadon.rds')
} else {
  PWADONf <- readRDS('./outputs/I4512ivac_fecal_all_pwadon.rds')
}

if (!file.exists('./outputs/I4512ivac_fecal_all_pwadon2.rds')){
  PWADONf2 <- pairwise.adonis(x = DISTf2, factors = METf2$set, perm =9999)
  write_rds(PWADONf2, './outputs/I4512ivac_fecal_all_pwadon2.rds')
} else {
  PWADONf2 <- readRDS('./outputs/I4512ivac_fecal_all_pwadon2.rds')
}

if (!file.exists('./outputs/I4512ivac_cecalcontents_all_pwadon.rds')){
  PWADONcc <- pairwise.adonis(x = DISTcc, factors = METcc$set, perm =9999)
  write_rds(PWADONcc, './outputs/I4512ivac_cecalcontents_all_pwadon.rds')
} else {
  PWADONcc <- readRDS('./outputs/I4512ivac_cecalcontents_all_pwadon.rds')
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

PCOMPS2 <- bind_rows(twdf2, totf2, ccc) #all important pairwise comparisons 
PCOMPS2$FDR_pval <- p.adjust(PCOMPS2$p.value, method = "fdr")
PCOMPS2 <- select(PCOMPS2, -sig, -p.adjusted)
PCOMPS2$FDR_sig[PCOMPS2$FDR_pval < 0.045] <- "*" #significant pairs after correction
PCOMPS2$FDR_sig[PCOMPS2$FDR_pval >= 0.045] <- "" 

## Save adjusted pairwise comparisons of twdf, totf, and ccc
write_tsv(PCOMPS, "./outputs/I4512ivac_pairwise_permanova_tests_all.tsv")
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

## PERMANOVA bar graphs 
pwdayALL %>% 
  ggplot(aes(x=F.Model, y=fct_rev(Group))) +
  geom_col() +
  facet_wrap(~Day, nrow = 1)
  ###Treatment group comparisons at each day
tiff(filename="./graphics/I4512ivac_permanova_groupsateachday.tiff", width=8, height=8, units="in", res=600, bg="transparent", family = "sans", pointsize = 12,  compression="lzw")
pwdayALL %>%
  ggplot(aes(x=F.Model, y=fct_rev(Group))) +
  geom_col() +
  facet_wrap(~Day, ncol = 1) +
  ggtitle('PERMANOVA tests between treatment groups at each day')
dev.off()
  ###With p-values and color by Group pairs
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
    ### Difference from MVMC (Mock_Mock) at each day
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
    ### Difference from MVC (Mock_I 4,[5],12:i:-) at each day
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
    ###No MVMC
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

## Within group changes over time (comparing back to themselves at day -28)
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