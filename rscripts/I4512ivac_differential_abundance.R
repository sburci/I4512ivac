# I 4,[5],12:i:- MICROBIOME ANALYSIS
# Selma Burciaga
# June 23, 2023
# Rscript: I4512ivac_differential_abundance.R
# R v4.3.0

## PART 5:
## anything ending with f is for Fecal samples
## anything ending with cc is for Cecal Content samples
## code in this script was created with the help of Julian Trachsel



# DIFFERENTIAL ABUNDANCE: STATS USING DESeq2___________________________________________________
# DESeq2 normalizes data, fits model, and runs statistical test to identify log2FoldChange in OTUs
library(tidyverse)
library(phyloseq)
library(cowplot)
library(DESeq2)
#library(ANCOMBC)
#library(funfuns)


#Open I4512ivac_phylo_all.rds and create phylo subset objects
PHYLO <- readRDS("./outputs/I4512ivac_phylo_all.rds")
PHYLOf <- prune_samples(samples = PHYLO@sam_data$Sample_Type=="Fecal", x = PHYLO)
PHYLOcc <- prune_samples(samples = PHYLO@sam_data$Sample_Type=="Cecal_contents", x = PHYLO)

#Convert treatment and dpc to factor
PHYLOf@sam_data$dpc <- factor(PHYLOf@sam_data$Days_post_challenge, levels=c('-28','-1','2', '3', '7', '10', '14'))
PHYLOf@sam_data$set <- paste(PHYLOf@sam_data$dpc, PHYLOf@sam_data$Treatment, sep = '_')
  #Confirm
PHYLOf@sam_data$Treatment
PHYLOf@sam_data$treat
PHYLOf@sam_data$Days_post_challenge
PHYLOf@sam_data$dpc
PHYLOf@sam_data$set

#Prune_taxa
PHYLOf %>% 
  prune_taxa(taxa = taxa_sums(.) > 0) %>%  
  prune_taxa(taxa=colSums(.@otu_table >0) > 5)

#Create DESeq2 object (dds)
dds <- phyloseq_to_deseq2(physeq = PHYLOf, design = ~ Treatment + dpc + Treatment:dpc)





#COMPARISONS BY DAY: MVMC VS MVC & MVMC VS EVC_________________________________________________________________________________________
PHYLOf@sam_data$treat <- factor(PHYLOf@sam_data$Treatment, levels=c('MVMC','MVC','EVC')) #MVMC as control

#Days post challenge -28
PHYLOf_d28 <- prune_samples(x = PHYLOf, samples = PHYLOf@sam_data$dpc == -28)
PHYLOf_d28 <- prune_taxa(taxa_sums(PHYLOf_d28) > 1, PHYLOf_d28)
dds_d28 <- phyloseq_to_deseq2(physeq = PHYLOf_d28, design = ~ treat)
dds_d28 <- DESeq(dds_d28, test = 'Wald', fitType = 'parametric') 
resultsNames(dds_d28) #results group comparisons and effect size estimation with apeglm for GLM

  #Comparison: mock vaccinated challenge (MVC) vs mock vaccinated mock challenged (MVMC)
  #res_d28 <- lfcShrink(dds_d28, coef = resultsNames(dds_d28)[2], type = 'normal') #apeglm is less bias than type = 'normal'
  res_d28_MVC <- lfcShrink(dds_d28, coef = 'treat_MVC_vs_MVMC', type = 'apeglm') 
  #create column with significant pvalues <0.05
  sig_d28_MVC = res_d28_MVC[which(res_d28_MVC$padj < 0.05), ]
  sig_d28_MVC = cbind(as(sig_d28_MVC, "data.frame"), as(tax_table(PHYLOf_d28)[rownames(sig_d28_MVC), ], "matrix"))
  sig_d28_MVC$Enrichedin <- ifelse(sig_d28_MVC$log2FoldChange <=0, 'MVMC', 'MVC') #negative log2Fold change = enriched in MVMC
  sig_d28_MVC$dpc <- '-28' #add column of dpc
  sig_d28_MVC$vs <- 'MVMCvsMVC' #add column of comparison
  #create log2FoldChange graph
  ggplot(sig_d28_MVC, aes(x = log2FoldChange, y = Genus)) +
    geom_point()

  #Comparison: Enterisol vaccinated challenge (EVC) vs mock vaccinated mock challenged (MVMC)
  res_d28_EVC <- lfcShrink(dds_d28, coef = 'treat_EVC_vs_MVMC', type = 'apeglm')
  #create column with significant pvalues <0.05
  sig_d28_EVC = res_d28_EVC[which(res_d28_EVC$padj < 0.05), ]
  sig_d28_EVC = cbind(as(sig_d28_EVC, "data.frame"), as(tax_table(PHYLOf_d28)[rownames(sig_d28_EVC), ], "matrix"))
  sig_d28_EVC$Enrichedin <- ifelse(sig_d28_EVC$log2FoldChange <=0, 'MVMC', 'EVC') #negative log2Fold change = enriched in MVMC
  sig_d28_EVC$dpc <- '-28' #add column of dpc
  sig_d28_EVC$vs <- 'MVMCvsEVC' #add column of comparison
  #create log2FoldChange graph
  ggplot(sig_d28_EVC, aes(x = log2FoldChange, y = Genus)) +
    geom_point() 

#Days post challenge -1
PHYLOf_d1 <- prune_samples(x = PHYLOf, samples = PHYLOf@sam_data$dpc == -1)
PHYLOf_d1 <- prune_taxa(taxa_sums(PHYLOf_d1) > 1, PHYLOf_d1)
dds_d1 <- phyloseq_to_deseq2(physeq = PHYLOf_d1, design = ~ treat)
dds_d1 <- DESeq(dds_d1, test = 'Wald', fitType = 'parametric') 
resultsNames(dds_d1) #results group comparisons and effect size estimation with apeglm for GLM

  #Comparison: mock vaccinated challenge (MVC) vs mock vaccinated mock challenged (MVMC)
  res_d1_MVC <- lfcShrink(dds_d1, coef = 'treat_MVC_vs_MVMC', type = 'apeglm') 
  #create column with significant pvalues <0.05
  sig_d1_MVC = res_d1_MVC[which(res_d1_MVC$padj < 0.05), ]
  sig_d1_MVC = cbind(as(sig_d1_MVC, "data.frame"), as(tax_table(PHYLOf_d1)[rownames(sig_d1_MVC), ], "matrix"))
  sig_d1_MVC$Enrichedin <- ifelse(sig_d1_MVC$log2FoldChange <=0, 'MVMC', 'MVC') #negative log2Fold change = enriched in MVMC
  sig_d1_MVC$dpc <- '-1' #add column of dpc
  sig_d1_MVC$vs <- 'MVMCvsMVC' #add column of comparison
  #create log2FoldChange graph
  ggplot(sig_d1_MVC, aes(x = log2FoldChange, y = Genus)) +
    geom_point()
  
  #Comparison: Enterisol vaccinated challenge (EVC) vs mock vaccinated mock challenged (MVMC)
  res_d1_EVC <- lfcShrink(dds_d1, coef = 'treat_EVC_vs_MVMC', type = 'apeglm')
  #create column with significant pvalues <0.05
  sig_d1_EVC = res_d1_EVC[which(res_d1_EVC$padj < 0.05), ]
  sig_d1_EVC = cbind(as(sig_d1_EVC, "data.frame"), as(tax_table(PHYLOf_d1)[rownames(sig_d1_EVC), ], "matrix"))
  sig_d1_EVC$Enrichedin <- ifelse(sig_d1_EVC$log2FoldChange <=0, 'MVMC', 'EVC') #negative log2Fold change = enriched in MVMC
  sig_d1_EVC$dpc <- '-1' #add column of dpc
  sig_d1_EVC$vs <- 'MVMCvsEVC' #add column of comparison
  #create log2FoldChange graph
  ggplot(sig_d1_EVC, aes(x = log2FoldChange, y = Genus)) +
    geom_point() 
  
#Days post challenge 2
PHYLOf_d2 <- prune_samples(x = PHYLOf, samples = PHYLOf@sam_data$dpc == 2)
PHYLOf_d2 <- prune_taxa(taxa_sums(PHYLOf_d2) > 1, PHYLOf_d2)
dds_d2 <- phyloseq_to_deseq2(physeq = PHYLOf_d2, design = ~ treat)
dds_d2 <- DESeq(dds_d2, test = 'Wald', fitType = 'parametric') 
resultsNames(dds_d2) #results group comparisons and effect size estimation with apeglm for GLM

  #Comparison: mock vaccinated challenge (MVC) vs mock vaccinated mock challenged (MVMC)
  res_d2_MVC <- lfcShrink(dds_d2, coef = 'treat_MVC_vs_MVMC', type = 'apeglm') 
  #create column with significant pvalues <0.05
  sig_d2_MVC = res_d2_MVC[which(res_d2_MVC$padj < 0.05), ]
  sig_d2_MVC = cbind(as(sig_d2_MVC, "data.frame"), as(tax_table(PHYLOf_d2)[rownames(sig_d2_MVC), ], "matrix"))
  sig_d2_MVC$Enrichedin <- ifelse(sig_d2_MVC$log2FoldChange <=0, 'MVMC', 'MVC') #negative log2Fold change = enriched in MVMC
  sig_d2_MVC$dpc <- '2' #add column of dpc
  sig_d2_MVC$vs <- 'MVMCvsMVC' #add column of comparison
  #create log2FoldChange graph
  ggplot(sig_d2_MVC, aes(x = log2FoldChange, y = Genus)) +
    geom_point()

  #Comparison: Enterisol vaccinated challenge (EVC) vs mock vaccinated mock challenged (MVMC)
  res_d2_EVC <- lfcShrink(dds_d2, coef = 'treat_EVC_vs_MVMC', type = 'apeglm')
  #create column with significant pvalues <0.05
  sig_d2_EVC = res_d2_EVC[which(res_d2_EVC$padj < 0.05), ]
  sig_d2_EVC = cbind(as(sig_d2_EVC, "data.frame"), as(tax_table(PHYLOf_d2)[rownames(sig_d2_EVC), ], "matrix"))
  sig_d2_EVC$Enrichedin <- ifelse(sig_d2_EVC$log2FoldChange <=0, 'MVMC', 'EVC') #negative log2Fold change = enriched in MVMC
  sig_d2_EVC$dpc <- '2' #add column of dpc
  sig_d2_EVC$vs <- 'MVMCvsEVC' #add column of comparison
  #create log2FoldChange graph
  ggplot(sig_d2_EVC, aes(x = log2FoldChange, y = Genus)) +
    geom_point() 

#Days post challenge 3
PHYLOf_d3 <- prune_samples(x = PHYLOf, samples = PHYLOf@sam_data$dpc == 3)
PHYLOf_d3 <- prune_taxa(taxa_sums(PHYLOf_d3) > 1, PHYLOf_d3)
dds_d3 <- phyloseq_to_deseq2(physeq = PHYLOf_d3, design = ~ treat)
dds_d3 <- DESeq(dds_d3, test = 'Wald', fitType = 'parametric') 
resultsNames(dds_d3) #results group comparisons and effect size estimation with apeglm for GLM
  
  #Comparison: mock vaccinated challenge (MVC) vs mock vaccinated mock challenged (MVMC)
  res_d3_MVC <- lfcShrink(dds_d3, coef = 'treat_MVC_vs_MVMC', type = 'apeglm') 
  #create column with significant pvalues <0.05
  sig_d3_MVC = res_d3_MVC[which(res_d3_MVC$padj < 0.05), ]
  sig_d3_MVC = cbind(as(sig_d3_MVC, "data.frame"), as(tax_table(PHYLOf_d3)[rownames(sig_d3_MVC), ], "matrix"))
  sig_d3_MVC$Enrichedin <- ifelse(sig_d3_MVC$log2FoldChange <=0, 'MVMC', 'MVC') #negative log2Fold change = enriched in MVMC
  sig_d3_MVC$dpc <- '3' #add column of dpc
  sig_d3_MVC$vs <- 'MVMCvsMVC' #add column of comparison
  #create log2FoldChange graph
  ggplot(sig_d3_MVC, aes(x = log2FoldChange, y = Genus)) +
    geom_point()
  
  #Comparison: Enterisol vaccinated challenge (EVC) vs mock vaccinated mock challenged (MVMC)
  res_d3_EVC <- lfcShrink(dds_d3, coef = 'treat_EVC_vs_MVMC', type = 'apeglm')
  #create column with significant pvalues <0.05
  sig_d3_EVC = res_d3_EVC[which(res_d3_EVC$padj < 0.05), ]
  sig_d3_EVC = cbind(as(sig_d3_EVC, "data.frame"), as(tax_table(PHYLOf_d3)[rownames(sig_d3_EVC), ], "matrix"))
  sig_d3_EVC$Enrichedin <- ifelse(sig_d3_EVC$log2FoldChange <=0, 'MVMC', 'EVC') #negative log2Fold change = enriched in MVMC
  sig_d3_EVC$dpc <- '3' #add column of dpc
  sig_d3_EVC$vs <- 'MVMCvsEVC' #add column of comparison
  #create log2FoldChange graph
  ggplot(sig_d3_EVC, aes(x = log2FoldChange, y = Genus)) +
    geom_point()   
    
#Days post challenge 7
PHYLOf_d7 <- prune_samples(x = PHYLOf, samples = PHYLOf@sam_data$dpc == 7)
PHYLOf_d7 <- prune_taxa(taxa_sums(PHYLOf_d7) > 1, PHYLOf_d7)
dds_d7 <- phyloseq_to_deseq2(physeq = PHYLOf_d7, design = ~ treat)
dds_d7 <- DESeq(dds_d7, test = 'Wald', fitType = 'parametric') 
resultsNames(dds_d7) #results group comparisons and effect size estimation with apeglm for GLM
  
  #Comparison: mock vaccinated challenge (MVC) vs mock vaccinated mock challenged (MVMC)
  res_d7_MVC <- lfcShrink(dds_d7, coef = 'treat_MVC_vs_MVMC', type = 'apeglm') 
  #create column with significant pvalues <0.05
  sig_d7_MVC = res_d7_MVC[which(res_d7_MVC$padj < 0.05), ]
  sig_d7_MVC = cbind(as(sig_d7_MVC, "data.frame"), as(tax_table(PHYLOf_d7)[rownames(sig_d7_MVC), ], "matrix"))
  sig_d7_MVC$Enrichedin <- ifelse(sig_d7_MVC$log2FoldChange <=0, 'MVMC', 'MVC') #negative log2Fold change = enriched in MVMC
  sig_d7_MVC$dpc <- '7' #add column of dpc
  sig_d7_MVC$vs <- 'MVMCvsMVC' #add column of comparison
  #create log2FoldChange graph
  ggplot(sig_d7_MVC, aes(x = log2FoldChange, y = Genus)) +
    geom_point()
  
  #Comparison: Enterisol vaccinated challenge (EVC) vs mock vaccinated mock challenged (MVMC)
  res_d7_EVC <- lfcShrink(dds_d7, coef = 'treat_EVC_vs_MVMC', type = 'apeglm')
  #create column with significant pvalues <0.05
  sig_d7_EVC = res_d7_EVC[which(res_d7_EVC$padj < 0.05), ]
  sig_d7_EVC = cbind(as(sig_d7_EVC, "data.frame"), as(tax_table(PHYLOf_d7)[rownames(sig_d7_EVC), ], "matrix"))
  sig_d7_EVC$Enrichedin <- ifelse(sig_d7_EVC$log2FoldChange <=0, 'MVMC', 'EVC') #negative log2Fold change = enriched in MVMC
  sig_d7_EVC$dpc <- '7' #add column of dpc
  sig_d7_EVC$vs <- 'MVMCvsEVC' #add column of comparison
  #create log2FoldChange graph
  ggplot(sig_d7_EVC, aes(x = log2FoldChange, y = Genus)) +
    geom_point()

#Days post challenge 10
PHYLOf_d10 <- prune_samples(x = PHYLOf, samples = PHYLOf@sam_data$dpc == 10)
PHYLOf_d10 <- prune_taxa(taxa_sums(PHYLOf_d10) > 1, PHYLOf_d10)
dds_d10 <- phyloseq_to_deseq2(physeq = PHYLOf_d10, design = ~ treat)
dds_d10 <- DESeq(dds_d10, test = 'Wald', fitType = 'parametric') 
resultsNames(dds_d10) #results group comparisons and effect size estimation with apeglm for GLM
  
  #Comparison: mock vaccinated challenge (MVC) vs mock vaccinated mock challenged (MVMC)
  res_d10_MVC <- lfcShrink(dds_d10, coef = 'treat_MVC_vs_MVMC', type = 'apeglm') 
  #create column with significant pvalues <0.05
  sig_d10_MVC = res_d10_MVC[which(res_d10_MVC$padj < 0.05), ]
  sig_d10_MVC = cbind(as(sig_d10_MVC, "data.frame"), as(tax_table(PHYLOf_d10)[rownames(sig_d10_MVC), ], "matrix"))
  sig_d10_MVC$Enrichedin <- ifelse(sig_d10_MVC$log2FoldChange <=0, 'MVMC', 'MVC') #negative log2Fold change = enriched in MVMC
  sig_d10_MVC$dpc <- '10' #add column of dpc
  sig_d10_MVC$vs <- 'MVMCvsMVC' #add column of comparison
  #create log2FoldChange graph
  ggplot(sig_d10_MVC, aes(x = log2FoldChange, y = Genus)) +
    geom_point()
  
  #Comparison: Enterisol vaccinated challenge (EVC) vs mock vaccinated mock challenged (MVMC)
  res_d10_EVC <- lfcShrink(dds_d10, coef = 'treat_EVC_vs_MVMC', type = 'apeglm')
  #create column with significant pvalues <0.05
  sig_d10_EVC = res_d10_EVC[which(res_d10_EVC$padj < 0.05), ]
  sig_d10_EVC = cbind(as(sig_d10_EVC, "data.frame"), as(tax_table(PHYLOf_d10)[rownames(sig_d10_EVC), ], "matrix"))
  sig_d10_EVC$Enrichedin <- ifelse(sig_d10_EVC$log2FoldChange <=0, 'MVMC', 'EVC') #negative log2Fold change = enriched in MVMC
  sig_d10_EVC$dpc <- '10' #add column of dpc
  sig_d10_EVC$vs <- 'MVMCvsEVC' #add column of comparison
  #create log2FoldChange graph
  ggplot(sig_d10_EVC, aes(x = log2FoldChange, y = Genus)) +
    geom_point()  
  
#Days post challenge 14
PHYLOf_d14 <- prune_samples(x = PHYLOf, samples = PHYLOf@sam_data$dpc == 14)
PHYLOf_d14 <- prune_taxa(taxa_sums(PHYLOf_d14) > 1, PHYLOf_d14)
dds_d14 <- phyloseq_to_deseq2(physeq = PHYLOf_d14, design = ~ treat)
dds_d14 <- DESeq(dds_d14, test = 'Wald', fitType = 'parametric') 
resultsNames(dds_d14) #results group comparisons and effect size estimation with apeglm for GLM
  
  #Comparison: mock vaccinated challenge (MVC) vs mock vaccinated mock challenged (MVMC)
  res_d14_MVC <- lfcShrink(dds_d14, coef = 'treat_MVC_vs_MVMC', type = 'apeglm') 
  #create column with significant pvalues <0.05
  sig_d14_MVC = res_d14_MVC[which(res_d14_MVC$padj < 0.05), ]
  sig_d14_MVC = cbind(as(sig_d14_MVC, "data.frame"), as(tax_table(PHYLOf_d14)[rownames(sig_d14_MVC), ], "matrix"))
  sig_d14_MVC$Enrichedin <- ifelse(sig_d14_MVC$log2FoldChange <=0, 'MVMC', 'MVC') #negative log2Fold change = enriched in MVMC
  sig_d14_MVC$dpc <- '14' #add column of dpc
  sig_d14_MVC$vs <- 'MVMCvsMVC' #add column of comparison
  #create log2FoldChange graph
  ggplot(sig_d14_MVC, aes(x = log2FoldChange, y = Genus)) +
    geom_point()
  
  #Comparison: Enterisol vaccinated challenge (EVC) vs mock vaccinated mock challenged (MVMC)
  res_d14_EVC <- lfcShrink(dds_d14, coef = 'treat_EVC_vs_MVMC', type = 'apeglm')
  #create column with significant pvalues <0.05
  sig_d14_EVC = res_d14_EVC[which(res_d14_EVC$padj < 0.05), ]
  sig_d14_EVC = cbind(as(sig_d14_EVC, "data.frame"), as(tax_table(PHYLOf_d14)[rownames(sig_d14_EVC), ], "matrix"))
  sig_d14_EVC$Enrichedin <- ifelse(sig_d14_EVC$log2FoldChange <=0, 'MVMC', 'EVC') #negative log2Fold change = enriched in MVMC
  sig_d14_EVC$dpc <- '14' #add column of dpc
  sig_d14_EVC$vs <- 'MVMCvsEVC' #add column of comparison
  #create log2FoldChange graph
  ggplot(sig_d14_EVC, aes(x = log2FoldChange, y = Genus)) +
    geom_point()  
  
  
  
  
  
#COMPARISONS BY DAY: EVC VS MVC - both challenged groups_________________________________________________________________________________________
#change order of treatments before lfcShrink (comparison always goes back to the first listed)
PHYLOf <- prune_samples(PHYLOf@sam_data$Treatment %in% c('MVC', 'EVC'), PHYLOf) #prune_taxa to exclude MVMC
PHYLOf@sam_data$treat <- factor(PHYLOf@sam_data$Treatment, levels=c('MVC','EVC')) 

#Days post challenge -28
PHYLOf_d28_both <- prune_samples(x = PHYLOf, samples = PHYLOf@sam_data$dpc == -28)
PHYLOf_d28_both <- prune_taxa(taxa_sums(PHYLOf_d28_both) > 1, PHYLOf_d28_both)
dds_d28_both <- phyloseq_to_deseq2(physeq = PHYLOf_d28_both, design = ~ treat)
dds_d28_both <- DESeq(dds_d28_both, test = 'Wald', fitType = 'parametric')   
resultsNames(dds_d28_both)

  #Comparison: Enterisol vaccinated challenge (EVC) vs mock vaccinated challenged (MVC)
  res_d28_both <- lfcShrink(dds_d28_both, coef = 'treat_EVC_vs_MVC', type = 'apeglm')
  #create column with significant pvalues <0.05
  sig_d28_both = res_d28_both[which(res_d28_both$padj < 0.05), ]
  sig_d28_both = cbind(as(sig_d28_both, "data.frame"), as(tax_table(PHYLOf_d28)[rownames(sig_d28_both), ], "matrix"))
  sig_d28_both $Enrichedin <- ifelse(sig_d28_both$log2FoldChange <=0, 'MVC', 'EVC') #negative log2Fold change = enriched in MVC
  sig_d28_both$dpc <- '-28'
  sig_d28_both$vs <- 'MVCvsEVC' #add column of comparison
  #create log2FoldChange graph
  ggplot(sig_d28_both, aes(x = log2FoldChange, y = Genus)) +
    geom_point()

#Days post challenge -1
PHYLOf_d1_both <- prune_samples(x = PHYLOf, samples = PHYLOf@sam_data$dpc == -1)
PHYLOf_d1_both <- prune_taxa(taxa_sums(PHYLOf_d1_both) > 1, PHYLOf_d1_both)
dds_d1_both <- phyloseq_to_deseq2(physeq = PHYLOf_d1_both, design = ~ treat)
dds_d1_both <- DESeq(dds_d1_both, test = 'Wald', fitType = 'parametric')   
resultsNames(dds_d1_both)

  #Comparison: Enterisol vaccinated challenge (EVC) vs mock vaccinated challenged (MVC)
  res_d1_both <- lfcShrink(dds_d1_both, coef = 'treat_EVC_vs_MVC', type = 'apeglm')
  #create column with significant pvalues <0.05
  sig_d1_both = res_d1_both[which(res_d1_both$padj < 0.05), ]
  sig_d1_both = cbind(as(sig_d1_both, "data.frame"), as(tax_table(PHYLOf_d1)[rownames(sig_d1_both), ], "matrix"))
  sig_d1_both $Enrichedin <- ifelse(sig_d1_both$log2FoldChange <=0, 'MVC', 'EVC') #negative log2Fold change = enriched in MVC
  sig_d1_both$dpc <- '-1'
  sig_d1_both$vs <- 'MVCvsEVC' #add column of comparison
  #create log2FoldChange graph
  ggplot(sig_d1_both, aes(x = log2FoldChange, y = Genus)) +
    geom_point()

#Days post challenge 2
PHYLOf_d2_both <- prune_samples(x = PHYLOf, samples = PHYLOf@sam_data$dpc == 2)
PHYLOf_d2_both <- prune_taxa(taxa_sums(PHYLOf_d2_both) > 1, PHYLOf_d2_both)
dds_d2_both <- phyloseq_to_deseq2(physeq = PHYLOf_d2_both, design = ~ treat)
dds_d2_both <- DESeq(dds_d2_both, test = 'Wald', fitType = 'parametric')   
resultsNames(dds_d2_both)

  #Comparison: mock vaccinated challenged (MVC) vs Enterisol vaccinated challenge (EVC)  
  res_d2_both <- lfcShrink(dds_d2_both, coef = 'treat_EVC_vs_MVC', type = 'apeglm')
  #create column with significant pvalues <0.05
  sig_d2_both = res_d2_both[which(res_d2_both$padj < 0.05), ]
  sig_d2_both = cbind(as(sig_d2_both, "data.frame"), as(tax_table(PHYLOf_d2)[rownames(sig_d2_both), ], "matrix"))
  sig_d2_both $Enrichedin <- ifelse(sig_d2_both$log2FoldChange <=0, 'MVC', 'EVC') #negative log2Fold change = enriched in MVC
  sig_d2_both$dpc <- '2'
  sig_d2_both$vs <- 'MVCvsEVC' #add column of comparison
  #create log2FoldChange graph
  ggplot(sig_d2_both, aes(x = log2FoldChange, y = Genus)) +
    geom_point()
  
#Days post challenge 3
PHYLOf_d3_both <- prune_samples(x = PHYLOf, samples = PHYLOf@sam_data$dpc == 3)
PHYLOf_d3_both <- prune_taxa(taxa_sums(PHYLOf_d3_both) > 1, PHYLOf_d3_both)
dds_d3_both <- phyloseq_to_deseq2(physeq = PHYLOf_d3_both, design = ~ treat)
dds_d3_both <- DESeq(dds_d3_both, test = 'Wald', fitType = 'parametric')   
resultsNames(dds_d3_both)

  #Comparison: mock vaccinated challenged (MVC) vs Enterisol vaccinated challenge (EVC)  
  res_d3_both <- lfcShrink(dds_d3_both, coef = 'treat_EVC_vs_MVC', type = 'apeglm')
  #create column with significant pvalues <0.05
  sig_d3_both = res_d3_both[which(res_d3_both$padj < 0.05), ]
  sig_d3_both = cbind(as(sig_d3_both, "data.frame"), as(tax_table(PHYLOf_d3)[rownames(sig_d3_both), ], "matrix"))
  sig_d3_both $Enrichedin <- ifelse(sig_d3_both$log2FoldChange <=0, 'MVC', 'EVC') #negative log2Fold change = enriched in MVC
  sig_d3_both$dpc <- '3'
  sig_d3_both$vs <- 'MVCvsEVC' #add column of comparison
  #create log2FoldChange graph
  ggplot(sig_d2_both, aes(x = log2FoldChange, y = Genus)) +
    geom_point()

#Days post challenge 7
PHYLOf_d7_both <- prune_samples(x = PHYLOf, samples = PHYLOf@sam_data$dpc == 7)
PHYLOf_d7_both <- prune_taxa(taxa_sums(PHYLOf_d7_both) > 1, PHYLOf_d7_both)
dds_d7_both <- phyloseq_to_deseq2(physeq = PHYLOf_d7_both, design = ~ treat)
dds_d7_both <- DESeq(dds_d7_both, test = 'Wald', fitType = 'parametric')   
resultsNames(dds_d7_both)
  
  #Comparison: Enterisol vaccinated challenge (EVC) vs mock vaccinated challenged (MVC)
  res_d7_both <- lfcShrink(dds_d7_both, 'treat_EVC_vs_MVC', type = 'apeglm')
  #create column with significant pvalues <0.05
  sig_d7_both = res_d7_both[which(res_d7_both$padj < 0.05), ]
  sig_d7_both = cbind(as(sig_d7_both, "data.frame"), as(tax_table(PHYLOf_d7)[rownames(sig_d7_both), ], "matrix"))
  sig_d7_both $Enrichedin <- ifelse(sig_d7_both$log2FoldChange <=0, 'MVC', 'EVC') #negative log2Fold change = enriched in MVC
  sig_d7_both$dpc <- '7'
  sig_d7_both$vs <- 'MVCvsEVC' #add column of comparison
  #create log2FoldChange graph
  ggplot(sig_d7_both, aes(x = log2FoldChange, y = Genus)) +
    geom_point()

#Days post challenge 10
PHYLOf_d10_both <- prune_samples(x = PHYLOf, samples = PHYLOf@sam_data$dpc == 10)
PHYLOf_d10_both <- prune_taxa(taxa_sums(PHYLOf_d10_both) > 1, PHYLOf_d10_both)
dds_d10_both <- phyloseq_to_deseq2(physeq = PHYLOf_d10_both, design = ~ treat)
dds_d10_both <- DESeq(dds_d10_both, test = 'Wald', fitType = 'parametric')   
resultsNames(dds_d10_both) 

dds_d10_both <- DESeq(dds_d10_both, test = 'Wald', fitType = 'parametric', minReplicatesForReplace=Inf)   


  #Comparison: Enterisol vaccinated challenge (EVC) vs mock vaccinated challenged (MVC)
  res_d10_both <- lfcShrink(dds_d10_both, coef = 'treat_EVC_vs_MVC', type = 'apeglm')
  #create column with significant pvalues <0.05
  sig_d10_both = res_d1_both[which(res_d10_both$padj < 0.05), ] # NA pvalues (https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#pvaluesNA)
  sig_d10_both = cbind(as(sig_d10_both, "data.frame"), as(tax_table(PHYLOf_d10)[rownames(sig_d10_both), ], "matrix"))
  sig_d10_both $Enrichedin <- ifelse(sig_d10_both$log2FoldChange <=0, 'MVC', 'EVC') #negative log2Fold change = enriched in MVC
  sig_d10_both$dpc <- '10'
  sig_d10_both$vs <- 'MVCvsEVC' #add column of comparison
  #create log2FoldChange graph
  ggplot(sig_d10_both, aes(x = log2FoldChange, y = Genus)) +
    geom_point()
    
#Days post challenge 14
PHYLOf_d14_both <- prune_samples(x = PHYLOf, samples = PHYLOf@sam_data$dpc == 14)
PHYLOf_d14_both <- prune_taxa(taxa_sums(PHYLOf_d14_both) > 1, PHYLOf_d14_both)
dds_d14_both <- phyloseq_to_deseq2(physeq = PHYLOf_d14_both, design = ~ treat)
dds_d14_both <- DESeq(dds_d14_both, test = 'Wald', fitType = 'parametric')   
resultsNames(dds_d14_both)
  
  #Comparison: Enterisol vaccinated challenge (EVC) vs mock vaccinated challenged (MVC)
  res_d14_both <- lfcShrink(dds_d14_both, 'treat_EVC_vs_MVC', type = 'apeglm')
  #create column with significant pvalues <0.05
  sig_d14_both = res_d14_both[which(res_d14_both$padj < 0.05), ]
  sig_d14_both = cbind(as(sig_d14_both, "data.frame"), as(tax_table(PHYLOf_d14)[rownames(sig_d14_both), ], "matrix"))
  sig_d14_both $Enrichedin <- ifelse(sig_d14_both$log2FoldChange <=0, 'EVC', 'MVC') #negative log2Fold change = enriched in MVC
  sig_d14_both$dpc <- '14'
  sig_d14_both$vs <- 'MVCvsEVC' #add column of comparison
  #create log2FoldChange graph
  ggplot(sig_d14_both, aes(x = log2FoldChange, y = Genus)) +
    geom_point()

  
  
  
#Days post challenge 14
PHYLOf_d14_both <- prune_samples(x = PHYLOf, samples = PHYLOf@sam_data$dpc == 14)
PHYLOf_d14_both <- prune_taxa(taxa_sums(PHYLOf_d14_both) > 1, PHYLOf_d14_both)
dds_d14_both <- phyloseq_to_deseq2(physeq = PHYLOf_d14_both, design = ~ treat)
dds_d14_both <- DESeq(dds_d14_both, test = 'Wald', fitType = 'parametric')   
resultsNames(dds_d14_both)
  
  #Comparison: Enterisol vaccinated challenge (EVC) vs mock vaccinated challenged (MVC)
  res_d14_both <- lfcShrink(dds_d14_both, 'treat_EVC_vs_MVC', type = 'apeglm')
  #create column with significant pvalues <0.05
  sig_d14_both = res_d14_both[which(res_d14_both$padj < 0.05), ]
  sig_d14_both = cbind(as(sig_d14_both, "data.frame"), as(tax_table(PHYLOf_d14)[rownames(sig_d14_both), ], "matrix"))
  sig_d14_both $Enrichedin <- ifelse(sig_d14_both$log2FoldChange <=0, 'EVC', 'MVC') #negative log2Fold change = enriched in MVC
  sig_d14_both$dpc <- '14'
  sig_d14_both$vs <- 'MVCvsEVC' #add column of comparison
  #create log2FoldChange graph
  ggplot(sig_d14_both, aes(x = log2FoldChange, y = Genus)) +
    geom_point()

  
  


#GENERATE LOG2FOLDCHANGE GRAPHS_________________________________________________________________________________________
#Join multiple data frames (missing sig_d10_both)
sig_all <- rbind(sig_d28_both, sig_d1_both, sig_d2_both, sig_d3_both, sig_d7_both, sig_d14_both, 
                 sig_d28_MVC, sig_d1_MVC, sig_d2_MVC, sig_d3_MVC, sig_d7_MVC, sig_d10_MVC, sig_d14_MVC, 
                 sig_d28_EVC, sig_d1_EVC, sig_d2_EVC, sig_d3_MVC, sig_d7_EVC, sig_d10_EVC, sig_d14_EVC)
sig_all$dpc <- factor(sig_all$dpc, levels=c('-28','-1','2', '3', '7', '10', '14'))
sig_all$dpc



#create log2FoldChange graphs comparing each treatment
MVMCvsEVC <- 
  sig_all %>% 
  filter(abs(log2FoldChange) > .25) %>% 
  filter(vs %in% "MVMCvsEVC") %>% 
  ggplot(aes(x = -log2FoldChange, y = Genus, color = dpc)) +
  geom_point(aes(shape = dpc), size = 2) +
  ggtitle('MVMC vs EVC') +
  guides(color=guide_legend("Days post challenge"), shape=guide_legend("Days post challenge")) +
  annotate("text", x=-6.3, y=2, label="Enriched in MVMC", size = 2.5) +
  annotate("text", x=12.6, y=2, label="Enriched in EVC", size = 2.5) +
  theme(legend.key=element_blank(),legend.background=element_blank()) +
  theme_bw() + geom_vline(xintercept=0) +
  theme(legend.position = 'bottom')
tiff('./graphics/I4512ivac_fecal_differential_abundance_MVMCvsEVC.tiff', units="in", width=8, height=6, res=600, family = "sans", pointsize = 12,  compression="lzw")
MVMCvsEVC
dev.off()

MVMCvsMVC <- 
  sig_all %>% 
  filter(abs(log2FoldChange) > .25) %>% 
  filter(vs %in% "MVMCvsMVC") %>% 
  ggplot(aes(x = log2FoldChange, y = Genus, color = dpc)) +
  geom_point(aes(shape = dpc), size = 2) +
  ggtitle('MVMC vs MVC') +
  guides(color=guide_legend("Days post challenge"), shape=guide_legend("Days post challenge")) +
  annotate("text", x=-2.9, y=2, label="Enriched in MVMC", size = 2.5) +
  annotate("text", x=11.2, y=2, label="Enriched in MVC", size = 2.5) +
  theme(legend.key=element_blank(),legend.background=element_blank()) +
  theme_bw() + geom_vline(xintercept=0) +
  theme(legend.position = 'bottom')
tiff('./graphics/I4512ivac_fecal_differential_abundance_MVMCvsMVC.tiff', units="in", width=8, height=6, res=600, family = "sans", pointsize = 12,  compression="lzw")
MVMCvsMVC 
dev.off()

EVCvsMVC <- 
  sig_all %>% 
  filter(abs(log2FoldChange) > .25) %>% 
  filter(vs %in% "EVCvsMVC") %>% 
  ggplot(aes(x = log2FoldChange, y = Genus, color = dpc)) +
  geom_point(aes(shape = dpc), size = 2) +
  ggtitle('EVC vs MVC') +
  guides(color=guide_legend("Days post challenge"), shape=guide_legend("Days post challenge")) +
  annotate("text", x=-2.5, y=2, label="Enriched in EVC", size = 2.5) +
  annotate("text", x=9.3, y=2, label="Enriched in MVC", size = 2.5) +
  theme(legend.key=element_blank(),legend.background=element_blank()) +
  theme_bw() + geom_vline(xintercept=0) +
  theme(legend.position = 'bottom')
tiff('./graphics/I4512ivac_fecal_differential_abundance_EVCvsMVC.tiff', units="in", width=8, height=6, res=600, family = "sans", pointsize = 12,  compression="lzw")
EVCvsMVC 
dev.off()
  
EVC2vsMVMC <- 
  sig_all %>% 
  filter(abs(log2FoldChange) > .25) %>% 
  filter(vs %in% "EVC2vsMVMC") %>% 
  ggplot(aes(x = log2FoldChange, y = Genus, color = dpc)) +
  geom_point(aes(shape = dpc), size = 2) +
  ggtitle('EVC vs MVMC') +
  guides(color=guide_legend("Days post challenge"), shape=guide_legend("Days post challenge")) +
  annotate("text", x=-2.5, y=2, label="Enriched in EVC", size = 2.5) +
  annotate("text", x=9.3, y=2, label="Enriched in MVMC", size = 2.5) +
  theme(legend.key=element_blank(),legend.background=element_blank()) +
  theme_bw() + geom_vline(xintercept=0) +
  theme(legend.position = 'bottom')
tiff('./graphics/I4512ivac_fecal_differential_abundance_EVC2vsMVMC.tiff', units="in", width=8, height=6, res=600, family = "sans", pointsize = 12,  compression="lzw")
EVC2vsMVMC 
dev.off()
  
EVC2vsMVMC 
MVMCvsEVC
  
  
  
  
  
  
  
  
  
  






  
  
  










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






#LefSe

##before starting script, run following commands on mothur
# count.groups(shared=I4512ivac.otu.shared) #size of smallest group: 3661.
# sub.sample(shared=I4512ivac.otu.shared, size=3661) #subsample for future analyses; output file I4512ivac.otu.0.03.subsample.shared
# count.groups(shared=I4512ivac.otu.0.03.subsample.shared) #confrim all are same size
# get.groups(shared=I4512ivac.otu.0.03.subsample.shared, design=mouse.time.design, sets=Fecal)
# lefse(shared=I4512ivac.otu.0.03.subsample.shared, design=/project/fsepru113/sburciaga/I4512ivac/mothuro/I4512ivac_design.txt)

# Carmen Wickware's code
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


#REFERENCES
#using 'apeglm' for LFC shrinkage.Cite:
# Zhu, A., Ibrahim, J.G., Love, M.I. (2018) Heavy-tailed prior distributions for sequence count data: removing the noise and preserving large differences. Bioinformatics. https://doi.org/10.1093/bioinformatics/bty895
